#=
Main workflow:
- Generate PFIM reference networks
- Assign body- and biomass
- Burn-in PFIM networks
- Calibrate ADBM + Niche to PFIM attractors
- Burn-in ADBM + Niche (single realisation each)
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
using JLD2
using pfim
using ProgressMeter
using SpeciesInteractionNetworks
using Statistics

include("src/internals.jl")
include("src/adbm.jl")

import Random
Random.seed!(66)

# --- storage ---
rows = Dict[]
topo_curve_store = DataFrame()
dyn_curve_store  = DataFrame()
species_store    = DataFrame()

# --- data ---
traits = CSV.read("data/community.csv", DataFrame)
feeding_rules = CSV.read("data/feeding_rules.csv", DataFrame)

# --- global params ---
n_networks = 2
t = 5000
survival_threshold = 1e-12

# --- calibration ---
adbm_survival_rate = 0.75     # fraction of species surviving burn-in
link_retention = 0.95

# --- body size distribution ---
global_dist = LogNormal(log(30), 1.5)

size_bounds = Dict(
    "tiny"       => (0.1, 10.0),
    "small"      => (10.0, 50.0),
    "medium"     => (50.0, 100.0),
    "large"      => (100.0, 300.0),
    "very_large" => (300.0, 500.0),
    "gigantic"   => (500.0, Inf)
)

# MAIN LOOP

for i in 1:n_networks

    # --- 1. Body sizes ---

    y = collect(String, traits.size)

    bodysize = [
        begin
            lo, hi = size_bounds[s]
            rand(truncated(global_dist, lo, hi))
        end
        for s in y
    ]

    traits[!, :bodymass] = bodysize
    biomass = bodysize .^ (-3/4)

    # --- 2. PFIM networks ---

    mass_rule = (res, con) -> con >= 0.5 * res ? 1 : 0

    pfim_down = PFIM(
        traits,
        feeding_rules;
        return_type = :matrix,
        y = 30.0,
        downsample = true
    )

    pfim_down_size = PFIM(
        traits,
        feeding_rules;
        return_type = :matrix,
        size_col = :bodymass,
        num_size_rule = mass_rule,
        y = 30.0,
        downsample = true
    )

    # --- 3. Realised PFIM networks (burn-in) --- 

    realised_networks = Dict()

    pfim_down_realised = realise_network(
        pfim_down;
        t = t,
        threshold = survival_threshold
    )

    pfim_down_size_realised = realise_network(
        pfim_down_size;
        bodymasses = bodysize,
        t = t,
        threshold = survival_threshold
    )

    if pfim_down_realised !== nothing
        realised_networks["down"] = pfim_down_realised
    end

    if pfim_down_size_realised !== nothing
        realised_networks["down_size"] = pfim_down_size_realised
    end

    # --- 4. Generate matches ADBM + Niche --- 

    matched_networks = copy(realised_networks)

    for (ref_name, ref_net) in realised_networks

        println("\n========================================")
        println("Reference: $ref_name")
        println("========================================")

        target_S = ref_net.S
        target_C = ref_net.C

        # --- latent correction for burn-in loss ---
        S_latent = ceil(Int, target_S / adbm_survival_rate)
        C_latent = target_C / link_retention

        println("Target S: $target_S → Latent S: $S_latent")
        println("Target C: $target_C → Latent C: $C_latent")

        # ADBM

        params = adbm_parameters(traits, bodysize)

        adbm_net = adbmmodel(traits, params, biomass)

        adbm_realised = realise_network(
            Int.(adbm_net.edges.edges);
            t = t,
            threshold = survival_threshold
        )

        if adbm_realised !== nothing
            matched_networks["ADBM_$ref_name"] = adbm_realised
        else
            @warn "ADBM failed for $ref_name"
        end

        # NICHE

        niche_fw = Foodweb(
            :niche;
            S = S_latent,
            C = C_latent
        )

        niche_realised = realise_network(
            Matrix(niche_fw.A);
            t = t,
            threshold = survival_threshold
        )

        if niche_realised !== nothing
            matched_networks["niche_$ref_name"] = niche_realised
        else
            @warn "Niche failed for $ref_name"
        end
    end

    # --- 5. Run simulations --- 

    for (net_name, realised) in matched_networks

        A = realised.A
        params = realised.params
        final_biomasses = realised.biomasses
        survivors = realised.survivors

        N = build_network(A)

        # species data

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
            metabolic_class = MC
        )

        append!(species_store, species_df)

        # topological extinctions

        topo_results = run_topological_extinctions(N, params)
        R_topo = compute_robustness(topo_results)

        topo_curves = Dict(
            k => extinction_breakdown(v)
            for (k, v) in topo_results
        )

        topo_df = export_curves(topo_curves, "topo_$net_name", i)
        append!(topo_curve_store, topo_df)

        # dynamic extinctions

        dyn_results = run_dynamic_extinctions(params, final_biomasses)
        R_dyn = compute_robustness(dyn_results)

        dyn_curves = Dict(
            k => extinction_breakdown(v)
            for (k, v) in dyn_results
        )

        dyn_df = export_curves(dyn_curves, "dyn_$net_name", i)
        append!(dyn_curve_store, dyn_df)

        # summary row

        row = Dict(
            :net_id => i,
            :net_type => net_name,
            :S => length(survivors),
            :C => sum(A) / (length(survivors)^2)
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

# --- 6. Save outputs ---

results_df = DataFrame(rows)
all_curve_df = vcat(topo_curve_store, dyn_curve_store)

CSV.write("outputs/paleo_robustness_summaries.csv", results_df)
CSV.write("outputs/paleo_extinction_curves.csv", all_curve_df)
CSV.write("outputs/paleo_species_metadata.csv", species_store)