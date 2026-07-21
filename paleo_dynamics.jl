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
include("src/niche_downsample.jl")
include("src/random_downsample.jl")
include("src/link_downsample.jl")

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
survival_threshold = 1e-12 
C_min = 0.05
C_max = 0.15

# Define the t-values we want to test
t_values = [5000]

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

# 1. Outer Loop: Iterate over different values of t
for t in t_values
    println("\n==========================================")
    println(">>> Running extinctions for t = $t ...")
    println("==========================================\n")

    # 2. Inner Loop: Process networks
    for i in 1:n_networks
        println(">>> Processing run $i of $n_networks (t = $t)...")

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

        #pfim_cont = PFIM(df, feeding_rules; return_type=:matrix)
        pfim_size = PFIM(df, feeding_rules; return_type=:matrix, size_col=:bodymass, num_size_rule=mass_rule)

        # Construct Foodweb objects directly
        #pfim_down = Foodweb(Matrix(downsample_network(pfim_cont, 2.5; target_co=C_targ, max_iter=100)))
        pfim_down_size = Foodweb(Matrix(downsample_network(pfim_size, 2.5; target_co=C_targ, max_iter=100)))

        # Create the Niche-Downsampled Network from pfim_cont
        pfim_niche_down = Foodweb(Matrix(downsample_niche_network(pfim_size, 1.0; target_co=C_targ, max_iter=100)))

        # Create the degree-Downsampled Network from pfim_cont
        pfim_link_down = Foodweb(Matrix(downsample_degree_product_network(pfim_size; target_co=C_targ, max_iter=100)))

        # Create the random-Downsampled Network from pfim_cont
        pfim_rand_down = Foodweb(Matrix(downsample_random_network(pfim_size; target_co=C_targ, max_iter=100)))

        niche_fw = Foodweb(:niche; S=size(pfim_size, 1), C=C_targ)

        prods = map(==("primary"), string.(df.tiering))
        atn_fw = Foodweb(lmatrix(df.species, df.bodymass, prods))

        # Consolidate all initial networks into one dictionary
        initial_networks = Dict(
            "down_link" => pfim_link_down,
            "down_power" => pfim_down_size,
            "down_niche" => pfim_niche_down,
            "down_rand" =>  pfim_rand_down,
            "niche" => niche_fw,
            "atn" => atn_fw
        )

        # Dictionary to keep initial S and C values for summary reference
        creation_metrics = Dict()

        for (net_name, fw) in initial_networks
            S_init = size(fw.A, 1)
            C_init = sum(fw.A) / (S_init^2)
            Int_init = intervality(fw.A)
            creation_metrics[net_name] = (S=S_init, C=C_init, I=Int_init)

            # Extract Species Metadata at creation
            # Create a temporary container for this specific record
            temp_metadata = DataFrame()

            # Let the function record metadata into our temporary DataFrame
            record_species_stage!(temp_metadata, i, net_name, "creation", fw.A, nothing, biomass)

            # Inject the current t-value into the temporary DataFrame
            temp_metadata[!, :t_val] .= t

            # Append the completed DataFrame to our global store
            append!(species_metadata_store, temp_metadata, cols=:union)
        end

        # --- 4. Burn-In & Realisation ---
        realised_networks = Dict()
        for (net_name, fw) in initial_networks
            if net_name ∈ ["niche", "niche_down"]
                realised = realise_network(
                    fw;
                    t=t,
                    threshold=survival_threshold
                )
            else
                realised = realise_network(
                    fw;
                    bodymasses=df.bodymass,
                    t=t,
                    threshold=survival_threshold
                )

            end
            if realised !== nothing
                realised_networks[net_name] = realised
            end
        end

        # Skip this replicate iteration if no networks successfully survived burn-in
        if isempty(realised_networks)
            @warn "Iteration $i (t = $t): no realised networks. Skipping."
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
            Int_realised = intervality(A_realised)

            # Stage 2: Extract Species Metadata Post Burn-in (Realised)
            # Create a temporary container for this specific record
            temp_metadata = DataFrame()

            # Let the function record metadata into our temporary DataFrame
            record_species_stage!(temp_metadata, i, net_name, "post_burn_in", A_realised, params, survivors, final_biomasses)

            # Inject the current t-value into the temporary DataFrame
            temp_metadata[!, :t_val] .= t

            # Append the completed DataFrame to our global store
            append!(species_metadata_store, temp_metadata, cols=:union)

            # Topological extinctions
            N = build_network(A_realised)
            topo_results = run_topological_extinctions(N, params)
            R_topo = compute_robustness(topo_results)

            topo_curves = Dict(k => extinction_breakdown(v) for (k, v) in topo_results)
            topo_df = export_curves(topo_curves, "topo_$net_name", i)
            topo_df[!, :t_val] .= t # Record t in curve outputs
            append!(topo_curve_store, topo_df)

            # Dynamic extinctions (Using dynamic t from outer loop)
            dyn_results = run_dynamic_extinctions(params, final_biomasses; t=t)
            R_dyn = compute_robustness(dyn_results)

            dyn_curves = Dict(k => extinction_breakdown(v) for (k, v) in dyn_results)
            dyn_df = export_curves(dyn_curves, "dyn_$net_name", i)
            dyn_df[!, :t_val] .= t # Record t in curve outputs
            append!(dyn_curve_store, dyn_df)

            # Compile Summary row
            row = Dict(
                :t_val => t, # Record t in the summaries
                :net_id => i,
                :net_type => net_name,
                :S_creation => creation_metrics[net_name].S,
                :C_creation => creation_metrics[net_name].C,
                :S_realised => S_realised,
                :C_realised => C_realised,
                :I_creation => creation_metrics[net_name].I,
                :I_realised => Int_realised
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
end

# --- 6. Save outputs ---

results_df = DataFrame(rows)
all_curve_df = vcat(topo_curve_store, dyn_curve_store)

# Ensure the t_val column exists on metadata before writing, handling creation step defaults
if !("t_val" in names(species_metadata_store))
    species_metadata_store[!, :t_val] .= missing
end

CSV.write("outputs/paleo_robustness_summaries_bodysize.csv", results_df)
CSV.write("outputs/paleo_extinction_curves_bodysize.csv", all_curve_df)
CSV.write("outputs/paleo_species_metadata_bodysize.csv", species_metadata_store)