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

    # --- 2. Build the 4 PFIM (+ Niche) networks ---
    mass_rule = (res, con) -> con >= 0.5 * res ? 1 : 0

    pfim_downsample = PFIM(traits, feeding_rules; 
                           return_type = :matrix, y = 15.0, downsample = true)

    pfim_downsample_contsize = PFIM(traits, feeding_rules; 
                                    return_type = :matrix,
                                    size_col = :bodymass,
                                    num_size_rule = mass_rule, 
                                    y = 15.0, downsample = true)

    parameters = adbm_parameters(traits, bodysize)
    N = adbmmodel(traits, parameters, biomass)
    adbm = Matrix(N.edges.edges)

    # for fun we can build a niche model (lets use pfim metweb as params)
    S = size(pfim_downsample, 1)
    L = sum(pfim_downsample)
    C = L / S^2
    niche_fw = Foodweb(:niche; S, C)
    niche = Matrix(niche_fw.A)

    # Put them in a structured container
    networks = Dict(
        "down" => pfim_downsample,
        "down_size" => pfim_downsample_contsize,
        "niche" => niche,
        "ADBM" => adbm
    )

    # --- 3. Run sims per a network type ---
    for (net_name, pfim_net) in networks

        A = pfim_net
        fw = Foodweb(A)
        S = size(A, 1)

        if net_name ∈ ["niche",
                       #"meta", 
                       "down", 
                       "ADBM"]

            # model params - use standard bodysize scaling
            params = default_model(
                fw,
                BodyMass(; Z = 10),
                ClassicResponse(; h = 2.0),
            )

        else

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
            
        end

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