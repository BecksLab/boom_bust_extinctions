#=
################################################################################
RUN EXTINCTION SIMULATIONS
################################################################################

Purpose
-------
Load previously generated realised ecological networks and perform extinction
experiments.

Network generation and burn-in are deliberately separated from this workflow.
This allows extinction scenarios to be modified, expanded, or rerun without
regenerating ecological networks.

Pipeline
--------
1. Load realised network archive
2. Select realised networks
3. Run topological extinction experiments
4. Run dynamic extinction experiments
5. Calculate robustness metrics
6. Export extinction curves and summaries

Inputs
------
realised_network_archive.jld2

Outputs
-------
extinction_summary.csv
extinction_curves.csv


################################################################################
=#

# --- 1. General setup ---------------------------------------------------------

using CSV
using DataFrames
using DifferentialEquations
using EcologicalNetworksDynamics
using Extinctions
using FoodWebTools
using JLD2
using ProgressMeter
using SpeciesInteractionNetworks
using Statistics


include("src/internals.jl")

import Random
Random.seed!(66)

# --- 2. Load realised network archive ----------------------------------------

archive = JLD2.load_object(
    "outputs/realised_network_archive.jld2"
)

network_metadata = archive["networks"]

# --- 3. Storage ---------------------------------------------------------------

summary_rows = Dict[]

topo_curve_store = DataFrame()

dyn_curve_store = DataFrame()

# --- 4. Extinction parameters -------------------------------------------------

# These are now independent of network generation

extinction_times = [5000]
survival_threshold = archive["metadata"]["survival_threshold"]

# --- 5. Run extinction simulations -------------------------------------------

for extinction_t in extinction_times

    println()
    println("==========================================")
    println("Running extinction simulations")
    println("Extinction time = $extinction_t")
    println("==========================================")
    println()


    for ((community, net_id, net_name), network) in network_metadata


        # Extract lineage components
        creation = network["creation"]
        realised = network["realised"]

        # Shared information
        params = realised["params"]
        biomasses = realised["biomasses"]

        # 1. TOPOLOGICAL EXTINCTIONS - CREATION

        A_creation = Matrix(creation["A"])

        N_creation = build_network(A_creation)


        topo_creation_results =
            run_topological_extinctions(
                N_creation,
                network["generation"]["bodymass"]
            )


        topo_creation_robustness =
            compute_robustness(
                topo_creation_results
            )


        topo_creation_curves = Dict(
            k => extinction_breakdown(v)
            for (k, v) in topo_creation_results
        )


        topo_creation_df =
            export_curves(
                topo_creation_curves,
                "topo_creation_$net_name",
                net_id
            )


        topo_creation_df.net_id .= net_id
        topo_creation_df.net_type .= net_name
        topo_creation_df.community .= community
        topo_creation_df.extinction_stage .= "creation"
        topo_creation_df.extinction_method .= "topological"
        topo_creation_df.extinction_time .= extinction_t


        append!(
            topo_curve_store,
            topo_creation_df,
            cols=:union
        )

        # 2. TOPOLOGICAL EXTINCTIONS - REALISED

        A_realised = Matrix(realised["A"])

        N_realised = build_network(A_realised)


        topo_realised_results =
            run_topological_extinctions(
                N_realised,
                params.body_mass
            )


        topo_realised_robustness =
            compute_robustness(
                topo_realised_results
            )


        topo_realised_curves = Dict(
            k => extinction_breakdown(v)
            for (k, v) in topo_realised_results
        )


        topo_realised_df =
            export_curves(
                topo_realised_curves,
                "topo_realised_$net_name",
                net_id
            )


        topo_realised_df.net_id .= net_id
        topo_realised_df.net_type .= net_name
        topo_realised_df.community .= community
        topo_realised_df.extinction_stage .= "realised"
        topo_realised_df.extinction_method .= "topological"
        topo_realised_df.extinction_time .= extinction_t


        append!(
            topo_curve_store,
            topo_realised_df,
            cols=:union
        )

        # 3. DYNAMIC EXTINCTIONS - REALISED

        dyn_results =
            run_dynamic_extinctions(
                params,
                biomasses;
                t=extinction_t
            )

        dyn_robustness =
            compute_robustness(
                dyn_results
            )

        dyn_curves = Dict(
            k => extinction_breakdown(v)
            for (k, v) in dyn_results
        )

        dyn_df =
            export_curves(
                dyn_curves,
                "dyn_realised_$net_name",
                net_id
            )

        dyn_df.net_id .= net_id
        dyn_df.net_type .= net_name
        dyn_df.community .= community
        dyn_df.extinction_stage .= "realised"
        dyn_df.extinction_method .= "dynamic"
        dyn_df.extinction_time .= extinction_t

        append!(
            dyn_curve_store,
            dyn_df,
            cols=:union
        )

        # One summary row per lineage

        summary_row = Dict(:net_id => net_id,
            :net_type => net_name,
            :extinction_time => extinction_t,
            :community => community,
            :S_creation => creation["S"],
            :C_creation => creation["C"],
            :S_realised => realised["S"],
            :C_realised => realised["C"])


        for (k, v) in topo_creation_robustness
            summary_row[
                Symbol("topo_creation_" * k)
            ] = v
        end


        for (k, v) in topo_realised_robustness
            summary_row[
                Symbol("topo_realised_" * k)
            ] = v
        end


        for (k, v) in dyn_robustness
            summary_row[
                Symbol("dyn_realised_" * k)
            ] = v
        end


        push!(
            summary_rows,
            summary_row
        )


    end
end



# --- 6. Export outputs --------------------------------------------------------


summary_df = DataFrame(
    summary_rows
)
curve_df = vcat(
    topo_curve_store,
    dyn_curve_store
)
CSV.write(
    "outputs/extinction_summary.csv",
    summary_df
)
CSV.write(
    "outputs/extinction_curves.csv",
    curve_df
)