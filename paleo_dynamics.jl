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
n_networks = 10                  # number of networks to make
t = 5000                           # relaxation time after perturbation
survival_threshold = 1e-12         # extinction threshold

for i in 1:n_networks

    # --- 1. Assign body sizes ---
    y = collect(String, traits.size)

    bodysize = (
        y ->
            y == "tiny" ? rand(Uniform(0.1, 10.0)) :
            y == "small" ? rand(Uniform(10.0, 50.0)) :
            y == "medium" ? rand(Uniform(50.0, 100.0)) :
            y == "large" ? rand(Uniform(100.0, 300.0)) :
            y == "very_large" ? rand(Uniform(300.0, 500.0)) :
            y == "gigantic" ? rand(Uniform(500.0, 700.0)) : y
    ).(y)

    traits[!, :bodymass] = bodysize

    # --- 2. Build the 4 PFIM networks ---
    mass_rule = (res, con) -> con >= 0.5 * res ? 1 : 0

    pfim_meta = PFIM(traits, feeding_rules; return_type = :matrix)
    pfim_downsample = PFIM(traits, feeding_rules; 
                           return_type = :matrix, y = 30.0, downsample = true)

    pfim_meta_contsize = PFIM(traits, feeding_rules; 
                              return_type = :matrix,
                              size_col = :bodymass,
                              num_size_rule = mass_rule)

    pfim_downsample_contsize = PFIM(traits, feeding_rules; 
                                    return_type = :matrix,
                                    size_col = :bodymass,
                                    num_size_rule = mass_rule, 
                                    y = 30.0, downsample = true)

    # Put them in a structured container
    networks = Dict(
        "meta" => pfim_meta,
        "down" => pfim_downsample,
        "meta_size" => pfim_meta_contsize,
        "down_size" => pfim_downsample_contsize
    )

    # --- 3. Run sims per a network type ---
    for (net_name, pfim_net) in networks

        A = pfim_net
        fw = Foodweb(A)
        S = size(A, 1)

        # body sizes aligned with network ordering
        size_by_index = zeros(Float64, S)

        for (sp, idx) in sp_index
            match = findfirst(==(sp), traits.species)
            if match !== nothing
                size_by_index[idx] = traits.bodymass[match]
            end
        end

        # model params
        params = default_model(
            fw,
            BodyMass(size_by_index),
            ClassicResponse(; h = 2.0),
        )

        A = params.A

        # initial biomass
        B0 = rand(Uniform(0.1, 1), S)

        # --- 4. Burn-in ---
        sol = simulate(params, B0, t;
            show_degenerated = false,
            callback = CallbackSet(
                extinction_callback(params, survival_threshold)
            )
        )

        final_biomasses = sol.u[end]

        # --- 5. prune network ---
        survivors, final_adj_matrix =
            extract_valid_network(A, final_biomasses, survival_threshold)

        if survivors === nothing
            @warn "All extinct: network $i ($net_name)"
            continue
        end

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

        N = build_network(final_adj_matrix)

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