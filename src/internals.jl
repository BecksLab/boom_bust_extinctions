"""
    build_network(adj_mat::Matrix{Bool})

Construct a `SpeciesInteractionNetwork` object from a square adjacency matrix.

This function converts a Boolean adjacency matrix into a 
`SpeciesInteractionNetwork{Unipartite{Symbol}, Binary{Bool}}` object, 
where species are assigned generic symbolic identifiers (`:1`, `:2`, ..., `:S`).

# Arguments
- `adj_mat::Matrix{Bool}`: Square adjacency matrix (S × S), where:
    - `adj_mat[i, j] = true` indicates an interaction from species `i` to `j`
    - Rows and columns correspond to the same ordered set of species

# Returns
- `SpeciesInteractionNetwork`: A unipartite binary interaction network

# Assumptions
- The matrix is square (same number of rows and columns)
- Species identities are not preserved (relabelled as generic symbols)

# Errors
- Throws `ArgumentError` if the matrix is not square
"""
function build_network(adj_mat::Matrix{Bool})

    if size(adj_mat, 1) != size(adj_mat, 2)
        throw(ArgumentError("Adjacency matrix must be square"))
    end

    spp_rich = size(adj_mat, 1)
    spp_list = Symbol.(1:spp_rich)

    edges = Binary(.!iszero.(adj_mat))
    nodes = Unipartite(spp_list)

    return SpeciesInteractionNetwork(nodes, edges)
end

"""
    extract_valid_network(A, B, threshold)

Extract the ecologically valid subnetwork based on biomass and connectivity.

Filters species based on:
1. Biomass threshold (alive vs extinct)
2. Network connectivity (removes isolated species)

# Arguments
- `A::Matrix`: Full adjacency matrix
- `B::Vector`: Species biomasses
- `threshold::Float64`: Extinction threshold

# Returns
- `(survivors, A_sub)`:
    - `survivors::Vector{Int}`: Indices of retained species
    - `A_sub::Matrix`: Reduced adjacency matrix
- `(nothing, nothing)` if no valid species remain

# Ecological Interpretation
- Removes extinct species (`B ≤ threshold`)
- Removes species with no trophic links (no prey AND no consumers)

# Notes
- Connectivity filtering ensures only functionally integrated species remain
- This matches assumptions used in topological extinction workflows
"""
function extract_valid_network(A, B, threshold)

    survivors = findall(x -> x > threshold, B)

    if isempty(survivors)
        return nothing, nothing
    end

    A_sub = A[survivors, survivors]

    has_prey = vec(sum(A_sub; dims=2)) .> 0
    has_consumer = vec(sum(A_sub; dims=1)) .> 0

    keep_mask = has_prey .| has_consumer
    survivors = survivors[keep_mask]

    if isempty(survivors)
        return nothing, nothing
    end

    return survivors, A[survivors, survivors]
end

"""
    run_topological_extinctions(N, params)

Run a suite of topological extinction sequences on a network.

Applies multiple extinction strategies based on structural metrics
and returns the resulting extinction sequences.

# Arguments
- `N`: A `SpeciesInteractionNetwork`
- `params`: Model parameters (used for body mass ordering)

# Returns
- `Dict{String, Any}`: Mapping from scenario name to extinction sequence

# Extinction Strategies
- `"degree_high"` / `"degree_low"`: Remove by connectivity
- `"vul_high"` / `"vul_low"`: Remove by vulnerability (in-degree)
- `"gen_high"` / `"gen_low"`: Remove by generality (out-degree)
- `"bm_high"` / `"bm_low"`: Remove by body mass
- `"rand_basal"`: Random removal of basal species only
- `"rand_consumer"`: Random removal of consumers only

# Notes
- These are purely structural (no dynamics)
- Sequences are compatible with `robustness()`
"""
function run_topological_extinctions(N, params)

    results = Dict()

    scenarios = Dict(
        "degree_high"   => extinction(N, "degree", true),
        "degree_low"    => extinction(N, "degree", false),
        "vul_high"      => extinction(N, "vulnerability", true),
        "vul_low"       => extinction(N, "vulnerability", false),
        "gen_high"      => extinction(N, "generality", true),
        "gen_low"       => extinction(N, "generality", false),
        "bm_high"       => extinction(N, Symbol.(sortperm(params.body_mass, rev=true))),
        "bm_low"        => extinction(N, Symbol.(sortperm(params.body_mass, rev=false))),
        "rand_basal"    => extinction(N; protect = :consumer),
        "rand_consumer" => extinction(N; protect = :basal)
    )

    for (name, Ns) in scenarios

        # primary = 1 species removed per step
        primary_counts = collect(0:length(Ns)-1)

        results[name] = (
            networks = Ns,
            primary = primary_counts
        )
    end

    return results
end

"""
    dynamic_extinction_adaptive(params, B_init; kwargs...)

Run adaptive dynamic extinction simulations with sequential perturbations.

At each step:
1. Identify currently viable species
2. Recompute network structure
3. Select a target species based on a chosen criterion
4. Force its extinction
5. Simulate dynamics for `t` time units
6. Record resulting network
7. Repeat until collapse

# Arguments
- `params`: Model parameters (includes adjacency matrix and dynamics)
- `B_init::Vector`: Initial species biomasses

# Keyword Arguments
- `criterion::Symbol`: Metric for species removal
    - `:degree`, `:generality`, `:vulnerability`, `:bodymass`, `:random`
- `descending::Bool`: Whether to remove highest (true) or lowest (false) value
- `t::Int`: Simulation time after each extinction event
- `survival_threshold::Float64`: Biomass threshold for extinction

# Returns
- `Vector`: Sequence of `SpeciesInteractionNetwork` objects

# Ecological Interpretation
- Combines press perturbations (species removal) with relaxation dynamics
- Captures cascading secondary extinctions
- Extinction order adapts to changing network structure

# Notes
- The full dynamical system is preserved (no species removal from ODE system)
- Only the recorded networks are reduced
- Multiple species may go extinct per step due to cascades

# Stopping Conditions
- All species extinct
- No connected subnetwork remains
"""
function dynamic_extinction_adaptive(params, B_init;
    criterion = :degree,
    descending = true,
    t = 5,
    survival_threshold = 1e-3,
    show_progress = true,
    debug = false)

    # --- initial full system ---
    B = copy(B_init)
    A_full = copy(params.A)

    active = collect(1:length(B))
    S0 = length(B)

    Nseq = SpeciesInteractionNetwork[]
    primary_count = Int[0]

    prog = show_progress ? Progress(S0, "Dynamic extinctions") : nothing

    step = 0

    while true

        step += 1

        # --- alive in reduced system ---
        alive_local = findall(B .> survival_threshold)

        if isempty(alive_local)
            break
        end

        alive_global = active[alive_local]

        # --- local adjacency ---
        A_sub = A_full[alive_local, alive_local]

        vulnerability = vec(sum(A_sub, dims=1))
        generality    = vec(sum(A_sub, dims=2))
        degree        = vulnerability .+ generality

        # --- candidate selection ---
        candidates = collect(eachindex(alive_local))

        if isempty(candidates)
            break
        end

        metric =
            criterion == :degree ? degree :
            criterion == :generality ? generality :
            criterion == :vulnerability ? vulnerability :
            criterion == :bodymass ? params.body_mass[alive_global] :
            nothing

        if metric === nothing
            sp_local = rand(candidates)
        else
            vals = metric[candidates]
            target = descending ? maximum(vals) : minimum(vals)
            tied = findall(x -> x == target, vals)
            sp_local = candidates[rand(tied)]
        end

        # --- map to global ---
        sp_local_index = alive_local[sp_local]
        sp_global = active[sp_local_index]

        # --- REMOVE SPECIES (state-space deletion) ---
        deleteat!(B, sp_local_index)
        deleteat!(active, sp_local_index)

        A_full = A_full[setdiff(1:end, sp_local_index), setdiff(1:end, sp_local_index)]

        # --- rebuild ecological model (IMPORTANT FIX) ---
        fw = Foodweb(A_full)   # or your constructor equivalent

        params_local = default_model(
            fw,
            BodyMass(; Z = 10),
            ClassicResponse(; h = 2),
        )

        # --- primary extinction tracking ---
        push!(primary_count, S0 - length(B))

        # --- simulate reduced system ---
        try
            sol = simulate(params_local, B, t;
                show_degenerated = false,
                callback = CallbackSet(
                    extinction_callback(params_local, survival_threshold)
                )
            )

            B = sol.u[end]

            if any(B .< -1e-8)
                @warn "Negative biomass detected (ODE instability)"
            end

            alive_check = count(B .> survival_threshold)
            if alive_check > length(B)
                error("Logical inconsistency in alive count")
            end

        catch e
            @warn "Simulation failed" exception=e
            break
        end

        # --- enforce extinction consistency ---
        B .= max.(B, 0.0)

        # --- rebuild network for output ---
        survivors, A_valid = extract_valid_network(A_full, B, survival_threshold)

        if survivors === nothing
            break
        end

        push!(Nseq, build_network(A_valid))

        if show_progress
            next!(prog)
        end

        if debug
            println("Step $step | remaining: $(length(B))")
        end
    end

    return (
        networks = Nseq,
        primary = primary_count
    )
end

"""
    run_dynamic_extinctions(params, B)

Run all adaptive dynamic extinction scenarios.

# Arguments
- `params`: Model parameters
- `B`: Initial biomasses after burn-in

# Returns
- `Dict{String, Vector}`: Mapping scenario name → network sequence

# Scenarios
- Degree (high/low)
- Vulnerability (high/low)
- Generality (high/low)
- Body mass (high/low)
- Random removal (consumer/basal only)

# Notes
- Each scenario uses adaptive re-ranking at every step
- Results are directly comparable to topological sequences
"""
function run_dynamic_extinctions(params, B; show_progress=true)

    results = Dict()

    scenarios = [
        ("degree_high",   :degree, true),
        ("degree_low",    :degree, false),
        ("vul_high",      :vulnerability, true),
        ("vul_low",       :vulnerability, false),
        ("gen_high",      :generality, true),
        ("gen_low",       :generality, false),
        ("bm_high",       :bodymass, true),
        ("bm_low",        :bodymass, false),
        ("rand_basal",    :random_basal, true),
        ("rand_consumer", :random_consumer, true)
    ]

    prog = show_progress ? Progress(length(scenarios), "Dynamic scenarios") : nothing

    for (name, crit, desc) in scenarios

        results[name] = dynamic_extinction_adaptive(
            params, B;
            criterion = crit,
            descending = desc,
            show_progress = show_progress
        )

        if show_progress
            next!(prog)
        end
    end

    return results
end

"""
    compute_robustness(results_dict)

Compute robustness metrics for a set of extinction sequences.

# Arguments
- `results_dict::Dict`: Mapping scenario name → network sequence

# Returns
- `Dict{String, Float64}`: Mapping scenario name → robustness value (e.g. R50)

# Notes
- Assumes all sequences are compatible with `robustness()`
- Works for both topological and dynamic extinction outputs
"""
function compute_robustness(results_dict)

    R = Dict()

    for (k, v) in results_dict
        # extract only networks
        R[k] = robustness(v.networks)
    end

    return R
end

function extinction_breakdown(result)

    Ns = result.networks
    primary_counts = result.primary

    n = min(length(Ns), length(primary_counts))

    Ns = Ns[1:n]
    primary_counts = primary_counts[1:n]

    # --- richness (state-based) ---
    richness_vals = SpeciesInteractionNetworks.richness.(Ns)

    @assert length(richness_vals) > 1

    if any(diff(richness_vals) .> 0)
        @warn "Non-monotonic richness detected (INVALID CURVE)"
    end

    S0 = richness_vals[1]

    if S0 <= 0
        error("Invalid S0 (zero or negative richness baseline)")
    end

    if any(richness_vals .> S0)
        @warn "Richness exceeds baseline S0 → normalization broken"
    end

    total_loss = (S0 .- richness_vals) ./ S0

    # --- primary loss (event-based) ---
    primary_prop = primary_counts ./ S0

    # enforce physical constraint:
    # primary cannot exceed observed loss
    primary_prop = min.(primary_prop, total_loss)

    # secondary becomes true residual
    secondary_prop = total_loss .- primary_prop

    if any(diff(total_loss) .< 0)
        @warn "Non-monotonic extinction curve detected"
    end

    return (
        primary = primary_prop,
        secondary = secondary_prop,
        total = total_loss
    )
end

function export_curves(curves_dict, type_label, net_id)

    rows = DataFrame()

    for (scenario, c) in curves_dict

        n = length(c.total)

        append!(rows, DataFrame(
            net_id = fill(net_id, n),
            type = fill(type_label, n),
            scenario = fill(scenario, n),
            step = 1:n,
            primary = c.primary,
            secondary = c.secondary,
            total = c.total
        ))
    end

    return rows
end