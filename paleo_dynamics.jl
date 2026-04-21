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
traits = CSV.read("data/traits.csv", DataFrame)
edges = CSV.read("data/edgelist.csv", DataFrame)

# --- Global Params ---
n_networks = 10                  # number of networks to make
t = 5000                           # relaxation time after perturbation
survival_threshold = 1e-12         # extinction threshold

# --- Build food web from edgelist ---

species = unique(vcat(edges.taxon_resource, edges.taxon_consumer))
S = length(species)

sp_index = Dict(sp => i for (i, sp) in enumerate(species))

A = zeros(Int64, S, S)

for row in eachrow(edges)
    i = sp_index[row.taxon_resource]
    j = sp_index[row.taxon_consumer]
    A[i, j] = 1.0
end

# --- Define size classes ---
size_map = Dict(row.species => row.size for row in eachrow(traits))

size_ranges = Dict(
    "tiny" => (-12.0, -8.0),
    "small" => (-8.0, -4.0),
    "medium" => (-4.0, 0.0),
    "large" => (0.0, 4.0)
)

for i in 1:n_networks

    # --- 1. Build Network ---
    fw = Foodweb(A)

    S = size(A, 1)

    # --- 2. Simulate body sizes ---
    TL = compute_trophic_levels(A)

    logM = zeros(S)

    β = 1.0          # trophic scaling (~10x per level)
    σ = 0.5          # noise

    for (sp, i) in sp_index
        TL_i = TL[i]
        cls = size_map[sp]

        lo, hi = size_ranges[cls]
        μ = (lo + hi) / 2

        logM[i] = clamp(
                        μ + β * (TL_i - 1) + rand(Normal(0, σ)),
                        lo, hi
                    )
    end

    M = 10 .^ logM

    # feed into the params df
    params = default_model(
        fw,
        BodyMass(M),
        ClassicResponse(;h = 2.0),
    )

    A = params.A

    # simulate initial biomass
    B0 = rand(Uniform(0.1, 1), S)

    # --- 3. Burn-in ---
    sol = simulate(params, B0, t;
        show_degenerated = false,
        callback = CallbackSet(
            extinction_callback(params, survival_threshold)
        )
    )

    final_biomasses = sol.u[end]

    # --- 4. Extract valid network ---
    survivors, final_adj_matrix = extract_valid_network(A, final_biomasses, survival_threshold)

    if survivors === nothing
        @warn "All species extinct after burn-in, skipping network $i"
        continue
    end

    # --- 5. Get species/network stats ---
    # subset vectors
    BM = params.M[survivors]
    TL = params.trophic.levels[survivors]
    MC = params.metabolic_class[survivors]

    # build species-level dataframe
    species_df = DataFrame(
        net_id = fill(i, length(survivors)),
        species_id = 1:length(survivors),   # reindexed within network
        original_id = survivors,
        body_mass = BM,
        trophic_level = TL,
        metabolic_class = MC,
    )

    append!(species_store, species_df)

    N = build_network(final_adj_matrix)

    # --- 6. Topological extinctions ---
    topo_results = run_topological_extinctions(N, params)
    # r50
    R_topo = compute_robustness(topo_results)
    # get extinction curves
    topo_curves = Dict()
    for (k, v) in topo_results
        topo_curves[k] = extinction_breakdown(v)

        if isempty(v.networks)
            @warn "Empty sequence for scenario $k in network $i"
        end
    end
    topo_df = export_curves(topo_curves, "topo", i)
    append!(topo_curve_store, topo_df)

    # --- 7. Dynamic extinctions ---
    dynamic_results = run_dynamic_extinctions(params, final_biomasses)
    # r50
    R_dyn = compute_robustness(dynamic_results)
    # get extinction curves
    dyn_curves = Dict()
    for (k, v) in dynamic_results
        dyn_curves[k] = extinction_breakdown(v)

        if isempty(v.networks)
            @warn "Empty sequence for scenario $k in network $i"
        end
    end
    dyn_df = export_curves(dyn_curves, "dyn", i)
    append!(dyn_curve_store, dyn_df)

    # --- 8. Build row ---
    row = Dict(
        :net_id => i,
        :S => S,
        :C => sum(A) / (S^2),
    )

    # --- Add topological results ---
    for (k, v) in R_topo
        row[Symbol("topo_" * k)] = v
    end

    # --- Add dynamic results ---
    for (k, v) in R_dyn
        row[Symbol("dyn_" * k)] = v
    end

    # --- Append to row dict ---
    push!(rows, row)

end

# create data frame object
results_df = DataFrame(rows)

# extinction curves
all_curve_df = vcat(topo_curve_store, dyn_curve_store)

# Write files
CSV.write("outputs/paleo_robustness_summaries.csv", results_df)
CSV.write("outputs/paleo_extinction_curves.csv", all_curve_df)
CSV.write("outputs/paleo_species_metadata.csv", species_store)