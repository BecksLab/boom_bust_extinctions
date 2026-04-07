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

# --- Global Params ---
t = 5                              # relaxation time after perturbation
survival_threshold = 1e-3          # extinction threshold
S = 10                             # number of species
C = 0.1                            # connectance

# --- 2. Build Network ---

fw = Foodweb(:niche; S, C)

params = default_model(
    fw,
    BodyMass(; Z = 10),
    ClassicResponse(; h = 2),
)

A = params.A

# Initial biomasses
B0 = rand(Uniform(0.1, 1), S)

# --- 3. Initial Burn-in ---

sol = simulate(params, B0, t;
    callback = CallbackSet(
        extinction_callback(params, survival_threshold)
    )
)

final_biomasses = sol.u[end]

# --- 4. Extract topologically relevant network ---
# alive and connected species

survivors, final_adj_matrix = extract_valid_network(A, final_biomasses, survival_threshold)

# Stop if system collapsed during burn-in
if survivors === nothing
    error("All species extinct after burn-in")
end

# Build network object
N = build_network(final_adj_matrix)

# --- 5. Topological Extinctions ---

topo_results = run_topological_extinctions(N, params)

# --- 6. Dynamic Extinctions ---

dynamic_results = run_dynamic_extinctions(params, final_biomasses)

# --- 7. Get Robustness (R50) ---

R_topo = compute_robustness(topo_results)
R_dyn  = compute_robustness(dynamic_results)
