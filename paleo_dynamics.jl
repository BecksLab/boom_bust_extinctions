#=
Main workflow:
- Generate PFIM reference networks
- Assign body- and biomass
- Burn-in PFIM networks
- Generate ONE niche model matching the realised PFIM structure
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
species_store = DataFrame()

# --- data ---
traits = CSV.read("data/community.csv", DataFrame)
feeding_rules = CSV.read("data/feeding_rules.csv", DataFrame)
# import plankton df (note this also has ecogenie biomass)
plankton = CSV.read("data/community_plankton.csv", DataFrame)

# --- global params ---
n_networks = 20
t = 5000
survival_threshold = 1e-12
C_min = 0.05
C_max = 0.15

# --- Co distribution ---
C_dist = truncated(Normal(0.15, 0.05), C_min, C_max);

# --- body size distribution ---
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

    # --- 0. Sample parameters ---
    C_targ = rand(C_dist)

    # --- 1. Body sizes ---

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

    select!(df, Not(:biomass, :bodymass))

    # --- 3. PFIM metawebs ---

    mass_rule = (res, con) -> con >= 0.5 * res ? 1 : 0

    # create metaweb - categorical
    pfim_cont = PFIM(
        df,
        feeding_rules;
        return_type=:matrix
    )

    # create metaweb - cont size
    pfim_size = PFIM(
        traits,
        feeding_rules;
        return_type=:matrix,
        size_col=:bodymass,
        num_size_rule=mass_rule
    )

    # --- 3. PFIM downsample ---

    pfim_down = downsample_network(pfim_cont, 2.5;
        target_co=C_targ,
        max_iter=100)

    pfim_down_size = downsample_network(pfim_size, 2.5;
        target_co=C_targ,
        max_iter=100)

    # --- 4. create niche web ---

    # use target Co and richness of pfim web
    niche_fw = Foodweb(
        :niche;
        S=size(pfim_cont, 1),
        C=C_targ
    )

    # --- 5. Realised networks (burn-in) ---

    realised_networks = Dict()

    pfim_down_realised = realise_network(
        pfim_down;
        t=t,
        threshold=survival_threshold
    )

    pfim_down_size_realised = realise_network(
        pfim_down_size;
        #bodymasses = bodysize,
        t=t,
        threshold=survival_threshold
    )

    niche_realised = realise_network(
        # this is ugly and janky AF!!
        Matrix(transpose(Matrix(niche_fw.A)));
        t=t,
        threshold=survival_threshold
    )

    if niche_realised !== nothing
        realised_networks["niche"] = niche_realised
    end

    if pfim_down_realised !== nothing
        realised_networks["down"] = pfim_down_realised
    end

    if pfim_down_size_realised !== nothing
        realised_networks["down_size"] = pfim_down_size_realised
    end

    # Skip this replicate if no networks successfully realised
    if isempty(realised_networks)
        @warn "Iteration $i: no realised networks. Skipping."
        continue
    end

    # build mini df that has the initial richness and Co values (pre burn in)
    pre_networks = Dict()

    pre_networks["niche"] = (
        S = size(pfim_cont, 1),
        C = C_targ
    )
    pre_networks["down"] = (
        S = size(pfim_down, 1),
        C = sum(pfim_down) / (size(pfim_down, 1)^2)
    )
    pre_networks["down_size"] = (
        S = size(pfim_down_size, 1),
        C = sum(pfim_down_size) / (size(pfim_down_size, 1)^2)
    )

    # --- 6. Run simulations ---

    for (net_name, realised) in realised_networks

        A = realised.A
        params = realised.params
        final_biomasses = realised.biomasses
        survivors = realised.survivors

        N = build_network(A)

        # species data

        BM = params.M[survivors]
        TL = params.trophic.levels[survivors]
        MC = params.metabolic_class[survivors]

        species_df = DataFrame(
            net_id=fill(i, length(survivors)),
            net_type=fill(net_name, length(survivors)),
            species_id=1:length(survivors),
            original_id=survivors,
            body_mass=BM,
            trophic_level=TL,
            metabolic_class=MC
        )

        append!(species_store, species_df)

        # topological extinctions

        topo_results = run_topological_extinctions(N, params)
        R_topo = compute_robustness(topo_results)

        topo_curves = Dict(
            k => extinction_breakdown(v)
            for (k, v) in topo_results
        )

        topo_df = export_curves(topo_curves, "topo_$net_name", i)
        append!(topo_curve_store, topo_df)

        # dynamic extinctions

        dyn_results = run_dynamic_extinctions(params, final_biomasses)
        R_dyn = compute_robustness(dyn_results)

        dyn_curves = Dict(
            k => extinction_breakdown(v)
            for (k, v) in dyn_results
        )

        dyn_df = export_curves(dyn_curves, "dyn_$net_name", i)
        append!(dyn_curve_store, dyn_df)

        # summary row

        row = Dict(
            :net_id => i,
            :net_type => net_name,
            :S_final => length(survivors),
            :C_final => sum(A) / (length(survivors)^2),
            :C_initial => pre_networks[net_name].C,
            :S_initial => pre_networks[net_name].S
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

# --- 7. Save outputs ---

results_df = DataFrame(rows)
all_curve_df = vcat(topo_curve_store, dyn_curve_store)

CSV.write("outputs/paleo_robustness_summaries.csv", results_df)
CSV.write("outputs/paleo_extinction_curves.csv", all_curve_df)
CSV.write("outputs/paleo_species_metadata.csv", species_store)