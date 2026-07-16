#=
Main workflow:
- Generate PFIM reference networks
- Assign body- and biomass
- Burn-in networks
- Topological extinctions
- Dynamic extinctions
=#

# --- 1. General Set-up ---

using CSV
using DataFrames
using DifferentialEquations
using Distributions
using EcologicalNetworksDynamics
using Extinctions
using FoodWebTools
using JLD2
using pfim
using ProgressMeter
using SpeciesInteractionNetworks
using Statistics

include("src/internals.jl")

import Random
Random.seed!(66)

# --- storage ---
rows = Dict[]
topo_curve_store = DataFrame()
dyn_curve_store = DataFrame()
species_metadata_store = DataFrame()

# --- data ---
traits = CSV.read("data/community.csv", DataFrame)
feeding_rules = CSV.read("data/feeding_rules.csv", DataFrame)
plankton = CSV.read("data/community_plankton.csv", DataFrame)

# --- global params ---
n_networks = 20
t = 50
survival_threshold = 1e-12
C_min = 0.05
C_max = 0.15

# --- Distributions ---
C_dist = truncated(Normal(0.15, 0.05), C_min, C_max)
global_dist = LogNormal(log(30), 1.5)

size_bounds = Dict(
    "primary" => (0.01, 0.1),
    "tiny" => (0.1, 10.0),
    "small" => (10.0, 50.0),
    "medium" => (50.0, 100.0),
    "large" => (100.0, 300.0),
    "very_large" => (300.0, 500.0),
    "gigantic" => (500.0, Inf)
)

# MAIN LOOP

for i in 1:n_networks
    println(">>> Processing run $i of $n_networks...")

    # --- 1. Sample parameters & Body sizes ---
    C_targ = rand(C_dist)
    y = collect(String, traits.size)

    bodysize = [
        begin
            lo, hi = size_bounds[s]
            rand(truncated(global_dist, lo, hi))
        end
        for s in y
    ]

    traits[!, :bodymass] = bodysize
    traits.biomass = fill(missing, nrow(traits))
    df = vcat(traits, plankton)

    # --- 2. Biomass estimates ---
    known = .!ismissing.(df.biomass)
    b = -3/4
    a = exp(mean(log.(df.biomass[known]) .- b .* log.(df.bodymass[known])))
    predicted = df.bodymass .^ b
    df.biomass[.!known] .= predicted[.!known]
    biomass = float.(df.biomass)

    # --- 3. Base Networks Generation (Creation) ---
    mass_rule = (res, con) -> con >= 0.5 * res ? 1 : 0

    pfim_cont = PFIM(df, feeding_rules; return_type=:matrix)
    pfim_size = PFIM(df, feeding_rules; return_type=:matrix, size_col=:bodymass, num_size_rule=mass_rule)

    # Construct Foodweb objects directly
    pfim_down = Foodweb(Matrix(transpose(downsample_network(pfim_cont, 2.5; target_co=C_targ, max_iter=100))))
    pfim_down_size = Foodweb(Matrix(transpose(downsample_network(pfim_size, 2.5; target_co=C_targ, max_iter=100))))
    niche_fw = Foodweb(:niche; S=size(pfim_cont, 1), C=C_targ)

    prods = map(==("primary"), string.(df.tiering))
    atn_fw = Foodweb(Matrix(transpose(lmatrix(df.species, df.bodymass, prods))))

    # Consolidate all initial networks into one dictionary
    initial_networks = Dict(
        "down" => pfim_down,
        "down_size" => pfim_down_size,
        "niche" => niche_fw,
        "atn" => atn_fw
    )

    # Dictionary to keep initial S and C values for summary reference
    creation_metrics = Dict()

    for (net_name, fw) in initial_networks
        S_init = size(fw.A, 1)
        C_init = sum(fw.A) / (S_init^2)
        creation_metrics[net_name] = (S=S_init, C=C_init)

        # Extract Species Metadata at creation
        record_species_stage!(species_metadata_store, i, net_name, "creation", fw.A, nothing)
    end

    # --- 4. Burn-In & Realisation ---
    realised_networks = Dict()
    for (net_name, fw) in initial_networks
        realised = realise_network(
            fw;
            t=t,
            threshold=survival_threshold
        )
        if realised !== nothing
            realised_networks[net_name] = realised
        end
    end

    # Skip this replicate iteration if no networks successfully survived burn-in
    if isempty(realised_networks)
        @warn "Iteration $i: no realised networks. Skipping."
        continue
    end

    # --- 5. Run Simulations ---
    for (net_name, realised) in realised_networks
        A_realised = realised.A
        params = realised.params
        final_biomasses = realised.biomasses
        survivors = realised.survivors

        S_realised = length(survivors)
        C_realised = sum(A_realised) / (S_realised^2)

        # Stage 2: Extract Species Metadata Post Burn-in (Realised)
        record_species_stage!(species_metadata_store, i, net_name, "post_burn_in", A_realised, params, survivors)

        # Topological extinctions
        N = build_network(A_realised)
        topo_results = run_topological_extinctions(N, params)
        R_topo = compute_robustness(topo_results)

        topo_curves = Dict(k => extinction_breakdown(v) for (k, v) in topo_results)
        topo_df = export_curves(topo_curves, "topo_$net_name", i)
        append!(topo_curve_store, topo_df)

        # Dynamic extinctions
        dyn_results = run_dynamic_extinctions(params, final_biomasses; t=t)
        R_dyn = compute_robustness(dyn_results)

        dyn_curves = Dict(k => extinction_breakdown(v) for (k, v) in dyn_results)
        dyn_df = export_curves(dyn_curves, "dyn_$net_name", i)
        append!(dyn_curve_store, dyn_df)

        # Compile Summary row
        row = Dict(
            :net_id => i,
            :net_type => net_name,
            :S_creation => creation_metrics[net_name].S,
            :C_creation => creation_metrics[net_name].C,
            :S_realised => S_realised,
            :C_realised => C_realised,
        )

        for (k, v) in R_topo
            row[Symbol("topo_" * k)] = v
        end

        for (k, v) in R_dyn
            row[Symbol("dyn_" * k)] = v
        end

        push!(rows, row)
    end
end

# --- 6. Save outputs ---

results_df = DataFrame(rows)
all_curve_df = vcat(topo_curve_store, dyn_curve_store)

CSV.write("outputs/paleo_robustness_summaries.csv", results_df)
CSV.write("outputs/paleo_extinction_curves.csv", all_curve_df)
CSV.write("outputs/paleo_species_metadata.csv", species_metadata_store)