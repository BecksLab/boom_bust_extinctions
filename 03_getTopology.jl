
# --- 1. General Set-up ---

using CSV
using DataFrames
using Extinctions
using FoodWebTools
using JLD2
using pfim
using ProgressMeter
using SpeciesInteractionNetworks
using Statistics

include("src/topology.jl")
include("src/internals.jl")

import Random
Random.seed!(66)

# import adjacency matrices

networks = load_object("outputs/adjacency_matrices.jld2")

# data frame

rows = Dict[]

for i in 1:nrow(networks)
    
    N = build_network(Matrix(networks.adjacency[i]))

    # push network summary
    d = _network_summary(N)
    d[:net_id] = networks.net_id[i]
    d[:net_type] = networks.net_type[i]
    d[:community] = networks.community[i]
    d[:stage] = networks.stage[i]

    # send to df
    push!(rows, d)

end

topology_df = DataFrame(rows)

CSV.write("outputs/paleo_topology.csv", topology_df)

# get gen and vul of each species from each network

deg_rows = DataFrame(
        spp_id = Any[],
        gen = Any[],
        vul = Any[],
        net_id = Any[],
        net_type = Any[],
        community = Any[],
        stage = Any[])

for i in 1:nrow(networks)
    
    N = build_network(Matrix(networks.adjacency[i]))

    gen = SpeciesInteractionNetworks.generality(N)
    vul = SpeciesInteractionNetworks.vulnerability(N)

    d = DataFrame(
        spp_id = collect(keys(gen)),
        gen = collect(values(gen)),
        vul = collect(values(vul)),
        net_id = fill(networks.net_id[i], length(gen)),
        net_type = fill(networks.net_type[i], length(gen)),
        stage = fill(networks.stage[i], length(gen)),
        community = fill(networks.community[i], length(gen))
    )

    # send to df
    append!(deg_rows, d)

end

CSV.write("outputs/paleo_spp_degrees.csv", deg_rows)