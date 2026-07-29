#=
################################################################################
BUILD REALISED NETWORKS
################################################################################

Purpose
-------
Generate ensembles of ecological networks, assign body sizes, realise each
network through burn-in dynamics, and archive the realised networks for later
analysis.

This script performs all computationally expensive operations exactly once.

Pipeline
--------
1. Load species traits and feeding rules
2. Sample body sizes
3. Generate network topologies
4. Burn-in each network (realisation)
5. Record metadata
6. Save a complete network archive (.jld2)

Outputs
-------
realised_network_archive.jld2

The archive contains

• realised networks
• creation networks
• species metadata
• adjacency matrices
• generation parameters
• summary statistics

The archive can subsequently be loaded by
02_extinctionSims.jl without regenerating any networks.

################################################################################
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

network_archive = Dict()

species_metadata_store = DataFrame()

adjacency_store = DataFrame(
    t=Int[],
    net_id=Int[],
    community=String[],
    net_type=String[],
    stage=String[],
    S=Int[],
    C=Float64[],
    adjacency=Any[]
)

summary_rows = Dict[]

# --- data ---
communities = Dict("dolomites" => Dict(
        "traits" => CSV.read("data/dolomites_community.csv", DataFrame),
        "plankton" => CSV.read("data/dolomites_plankton.csv", DataFrame)
    ), "russia" => Dict(
        "traits" => CSV.read("data/russia_community.csv", DataFrame),
        "plankton" => CSV.read("data/russia_plankton.csv", DataFrame)
    ))
feeding_rules = CSV.read("data/feeding_rules.csv", DataFrame)

# --- global params ---
n_networks = 30  # number of runs per network
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


# --- 2. Build and realise networks ---

for t in t_values

    # 2. Inner Loop: Process networks
    for i in 1:n_networks

        C_targ = rand(C_dist)

        for (community_name, community_data) in communities

            traits = deepcopy(community_data["traits"])
            plankton = deepcopy(community_data["plankton"])

            println(">>> Processing run $i of $n_networks (t = $t) for $community_name...")

            # --- 1. Sample parameters & Body sizes ---
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
            #known = .!ismissing.(df.biomass)
            #b = -3/4
            #a = exp(mean(log.(df.biomass[known]) .- b .* log.(df.bodymass[known])))
            #predicted = df.bodymass .^ b
            #df.biomass[.!known] .= predicted[.!known]
            #biomass = float.(df.biomass)

            # create initial biomass from uniform rand dist
            B_init = rand(nrow(df))

            # --- 3. Base Networks Generation (Creation) ---
            mass_rule = (res, con) -> con >= 0.5 * res ? 1 : 0

            pfim_meta = PFIM(df, feeding_rules; return_type=:matrix)
            #pfim_meta = PFIM(df, feeding_rules; return_type=:matrix, size_col=:bodymass, num_size_rule=mass_rule)

            # Construct Foodweb objects directly
            #pfim_down = Foodweb(Matrix(downsample_network(pfim_cont, 2.5; target_co=C_targ, max_iter=100)))
            pfim_power_down = Foodweb(Matrix(downsample(pfim_meta, :powerlaw; y=2.5, target_co=C_targ, max_iter=300)))

            # Create the Niche-Downsampled Network from pfim_cont
            pfim_niche_down = Foodweb(Matrix(downsample(pfim_meta, :niche; sigma_scale=1.0, target_co=C_targ, max_iter=300)))

            # Create the degree-Downsampled Network from pfim_cont
            pfim_link_down = Foodweb(Matrix(downsample(pfim_meta, :degree; target_co=C_targ, max_iter=300)))

            # Create the random-Downsampled Network from pfim_cont
            pfim_rand_down = Foodweb(Matrix(downsample(pfim_meta, :random; target_co=C_targ, max_iter=300)))

            niche_fw = Foodweb(:niche; S=size(pfim_meta, 1), C=C_targ)

            prods = map(==("primary"), string.(df.tiering))
            atn_fw = Foodweb(lmatrix(df.species, df.bodymass, prods; threshold=0.06))

            # Consolidate all initial networks into one dictionary
            initial_networks = Dict(
                "down_link" => pfim_link_down,
                "down_power" => pfim_power_down,
                "down_niche" => pfim_niche_down,
                "down_rand" => pfim_rand_down,
                "niche" => niche_fw,
                "atn" => atn_fw,
                "metaweb" => Foodweb(pfim_meta)
            )

            # Dictionary to keep initial S and C values for summary reference
            creation_metrics = Dict()

            for (net_name, fw) in initial_networks

                # summary metrics
                S_init = size(fw.A, 1)
                C_init = sum(fw.A) / (S_init^2)
                creation_metrics[net_name] = (S=S_init, C=C_init)

                # Extract Species Metadata at creation
                # Create a temporary container for this specific record
                temp_metadata = DataFrame()

                # Let the function record metadata into our temporary DataFrame
                record_species_stage!(temp_metadata, i, net_name, "creation", fw.A, df.bodymass, nothing, nothing, B_init)

                # Inject the current t-value into the temporary DataFrame
                temp_metadata[!, :t_val] .= t
                temp_metadata[!, :community] .= community_name

                # Append the completed DataFrame to our global store
                append!(species_metadata_store, temp_metadata, cols=:union)

                # Save the adjacency matrix as well
                push!(adjacency_store, (
                    t=t,
                    net_id=i,
                    community=community_name,
                    net_type=net_name,
                    stage="creation",
                    S=size(fw.A, 1),
                    C=sum(fw.A) / size(fw.A, 1)^2,
                    adjacency=copy(fw.A)
                ))

            end

            # --- 4. Burn-In & Realisation ---
            realised_networks = Dict()

            for (net_name, fw) in initial_networks

                # we treat this model set as not having prior biomasses
                if net_name ∈ ["niche", "down_niche", "metaweb"]
                    realised = realise_network(
                        fw;
                        t=t,
                        threshold=survival_threshold,
                        B0=B_init
                    )
                else
                    realised = realise_network(
                        fw;
                        bodymasses=df.bodymass,
                        t=t,
                        threshold=survival_threshold,
                        B0=B_init
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

            # Record realised networks and associated metadata

            for (net_name, realised) in realised_networks

                A_realised = realised.A
                params = realised.params
                final_biomasses = realised.biomasses
                survivors = realised.survivors

                S_realised = length(survivors)
                C_realised = sum(A_realised) / (S_realised^2)

                # Save realised adjacency matrix

                push!(adjacency_store, (
                    t=t,
                    net_id=i,
                    community=community_name,
                    net_type=net_name,
                    stage="burnin",
                    S=S_realised,
                    C=C_realised,
                    adjacency=copy(A_realised)
                ))

                # Record species metadata after burn-in

                temp_metadata = DataFrame()

                record_species_stage!(
                    temp_metadata,
                    i,
                    net_name,
                    "post_burn_in",
                    A_realised,
                    df.bodymass,
                    params,
                    survivors,
                    final_biomasses
                )

                temp_metadata.t_val .= t
                temp_metadata.community .= community_name

                append!(species_metadata_store, temp_metadata, cols=:union)

                # Network summary

                row = Dict(
                    :t_val => t,
                    :net_id => i,
                    :community => community_name,
                    :net_type => net_name,
                    :C_target => C_targ, :S_creation => creation_metrics[net_name].S,
                    :C_creation => creation_metrics[net_name].C,
                    :C_realised => C_realised
                )

                push!(summary_rows, row)

                # Archive realised network

                network_archive[(community_name, i, net_name)] = Dict("generation" => Dict(
                        "community" => community_name,
                        "burnin_time" => t,
                        "C_target" => C_targ,
                        "bodymass" => copy(df.bodymass),
                        "B0" => copy(B_init)
                    ), "creation" => Dict(
                        "A" => copy(initial_networks[net_name].A),
                        "S" => creation_metrics[net_name].S,
                        "C" => creation_metrics[net_name].C
                    ), "realised" => Dict(
                        "A" => copy(A_realised),
                        "params" => params,
                        "biomasses" => final_biomasses,
                        "survivors" => survivors,
                        "S" => S_realised,
                        "C" => C_realised
                    )
                )
            end
        end
    end
end


# --- 3. Collate and export relevant data ---

archive = Dict(
    # Global metadata
    "metadata" => Dict(
        "seed" => 66,
        "n_networks" => n_networks,
        "survival_threshold" => survival_threshold,
        "burnin_times" => t_values
    ),
    # Realised network collection
    "networks" => network_archive,
    # Species metadata
    "species_metadata" => species_metadata_store,
    # Stored adjacency matrices
    "adjacency_store" => adjacency_store,
    # Network summaries
    "summary" => DataFrame(summary_rows)
)

JLD2.save_object(
    "outputs/realised_network_archive.jld2",
    archive
)

CSV.write("outputs/paleo_species_metadata.csv", archive["species_metadata"])
JLD2.save_object("outputs/adjacency_matrices.jld2", archive["adjacency_store"])