#=
Main workflow:
- Generate networks; 
- Assign body- and bio-mass; 
- First burn-in; 
- Topological Extinctions;  
- Dynamic Extinctions; 
=#


# --- 1. Load Dependencies ---
using CSV
using EcologicalNetworksDynamics
using JLD2
using Statistics

# set seed
import Random
Random.seed!(66)

# --- 2. Build Networks ---

fw1 = Foodweb(:niche; S = 5, C = 0.2)