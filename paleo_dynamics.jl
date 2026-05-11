#=
Main workflow:
- Generate networks; 
- Assign body- and bio-mass; 
- First burn-in; 
- Topological Extinctions;  
- Dynamic Extinctions; 
=#

# --- 1. General Set-up ---

# Load Dependencies
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

# call internal functions
include("src/internals.jl")
include("src/adbm.jl")

# set seed
import Random
Random.seed!(66);

# output dict
rows = Dict[]
topo_curve_store = DataFrame()
dyn_curve_store  = DataFrame()
species_store = DataFrame()

# import dataframes
traits = CSV.read("data/community.csv", DataFrame)
feeding_rules = CSV.read("data/feeding_rules.csv", DataFrame)

# --- Global Params ---
n_networks = 2                  # number of networks to make
t = 5000                           # relaxation time after perturbation
survival_threshold = 1e-12         # extinction threshold

# --- Pruning Params ---
n_candidates = 5
S_tolerance = 4
C_tolerance = 0.1

for i in 1:n_networks

    # --- 1. Assign body sizes ---
    y = collect(String, traits.size)

    # Global body-size distribution
    global_dist = LogNormal(log(30), 1.5)
    # median ≈ 30
    # sigma controls spread/right-tail

    # Category bounds
    size_bounds = Dict(
        "tiny"       => (0.1, 10.0),
        "small"      => (10.0, 50.0),
        "medium"     => (50.0, 100.0),
        "large"      => (100.0, 300.0),
        "very_large" => (300.0, 500.0),
        "gigantic"   => (500.0, Inf)
    )

    # Sample from truncated log-normal
    bodysize = [
        begin
            lo, hi = size_bounds[s]
            rand(truncated(global_dist, lo, hi))
        end
        for s in y
    ]

    traits[!, :bodymass] = bodysize

    # Estimate biomass using Metabolic Theory scaling (M^-3/4)
    biomass = bodysize .^ (-3 / 4)

    # --- 2. Generate PFIM reference networks ---

    mass_rule = (res, con) -> con >= 0.5 * res ? 1 : 0

    pfim_downsample = PFIM(
        traits,
        feeding_rules;
        return_type = :matrix,
        y = 30.0,
        downsample = true
    )

    pfim_downsample_contsize = PFIM(
        traits,
        feeding_rules;
        return_type = :matrix,
        size_col = :bodymass,
        num_size_rule = mass_rule,
        y = 30.0,
        downsample = true
    )

    # --- 3. Realise PFIM networks ---

    realised_networks = Dict()

    # PFIM without explicit body sizes
    pfim_realised = realise_network(
        pfim_downsample;
        t = t,
        threshold = survival_threshold
    )

    if pfim_realised !== nothing
        realised_networks["down"] = pfim_realised
    end

    # PFIM with continuous body sizes
    pfim_size_realised = realise_network(
        pfim_downsample_contsize;
        bodymasses = bodysize,
        t = t,
        threshold = survival_threshold
    )

    if pfim_size_realised !== nothing
        realised_networks["down_size"] = pfim_size_realised
    end

    # --- 4. Select ecological reference network ---

    # Use size-structured PFIM as ecological target

    reference_name = "down_size"
    reference_net = realised_networks[reference_name]

    target_S = reference_net.S
    target_C = reference_net.C

    println("\n=================================================")
    println("Ecological reference: $reference_name")
    println("=================================================")

    println("Target realised richness: $target_S")
    println("Target realised connectance: $target_C")

    # Estimate generation parameters

    # Burn-in generally reduces richness/connectance

    S_generate = target_S + 5
    C_generate = target_C * 1.2

    # Generate matched niche network

    println("\nSearching for matching niche network...")

    best_niche = find_matching_network(

        () -> begin

            niche_fw = Foodweb(
                :niche;
                S = S_generate,
                C = C_generate
            )

            niche_A = Matrix(niche_fw.A)

            realise_network(
                niche_A;
                t = t,
                threshold = survival_threshold
            )

        end,

        target_S,
        target_C;

        max_attempts = n_candidates,
        S_tol = S_tolerance,
        C_tol = C_tolerance
    )

    # Generate matched ADBM network

    println("\nSearching for matching ADBM network...")

    best_adbm = find_matching_network(

        () -> begin

            parameters = adbm_parameters(
                traits,
                bodysize
            )

            N = adbmmodel(
                traits,
                parameters,
                biomass
            )

            adbm_A = Matrix(N.edges.edges)

            realise_network(
                adbm_A;
                t = t,
                threshold = survival_threshold
            )

        end,

        target_S,
        target_C;

        max_attempts = n_candidates,
        S_tol = S_tolerance,
        C_tol = C_tolerance
    )

    # Final network collection

    matched_networks = Dict(

        "down" => realised_networks["down"],
        "down_size" => realised_networks["down_size"],
        "niche" => best_niche,
        "ADBM" => best_adbm
    )

    # --- 5. Run sims per a network type ---
    for (net_name, realised) in matched_networks

        A = realised.A
        params = realised.params
        final_biomasses = realised.biomasses
        survivors = realised.survivors

        final_adj_matrix = realised.A

        N = build_network(final_adj_matrix)

        # --- 6. species stats ---
        BM = params.M[survivors]
        TL = params.trophic.levels[survivors]
        MC = params.metabolic_class[survivors]

        species_df = DataFrame(
            net_id = fill(i, length(survivors)),
            net_type = fill(net_name, length(survivors)),
            species_id = 1:length(survivors),
            original_id = survivors,
            body_mass = BM,
            trophic_level = TL,
            metabolic_class = MC,
        )

        append!(species_store, species_df)

        # --- 7. Topological extinctions ---
        topo_results = run_topological_extinctions(N, params)
        R_topo = compute_robustness(topo_results)

        topo_curves = Dict()
        for (k, v) in topo_results
            topo_curves[k] = extinction_breakdown(v)
        end

        topo_df = export_curves(topo_curves, "topo_$net_name", i)
        append!(topo_curve_store, topo_df)

        # --- 8. dynamic extinctions ---
        dynamic_results = run_dynamic_extinctions(params, final_biomasses)
        R_dyn = compute_robustness(dynamic_results)

        dyn_curves = Dict()
        for (k, v) in dynamic_results
            dyn_curves[k] = extinction_breakdown(v)
        end

        dyn_df = export_curves(dyn_curves, "dyn_$net_name", i)
        append!(dyn_curve_store, dyn_df)

        # --- 9. summary row ---
        row = Dict(
            :net_id => i,
            :net_type => net_name,
            :S => length(survivors),
            :C => sum(final_adj_matrix) / (length(survivors)^2),
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

# create data frame object
results_df = DataFrame(rows)

# extinction curves
all_curve_df = vcat(topo_curve_store, dyn_curve_store)

# Write files
CSV.write("outputs/paleo_robustness_summaries.csv", results_df)
CSV.write("outputs/paleo_extinction_curves.csv", all_curve_df)
CSV.write("outputs/paleo_species_metadata.csv", species_store)