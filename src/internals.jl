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
function dynamic_extinction_adaptive(params, B0;
    criterion = :degree,
    descending = true,
    t = 5000,
    survival_threshold = 1e-3,
    show_progress = true,
    debug = false
)

    S0 = length(B0)

    A = params.A
    B = copy(B0)

    alive = trues(S0)

    Nseq = SpeciesInteractionNetwork[]

    # --- FIX: track cumulative primary deletions ---
    primary_counts = Int[0]

    prog = show_progress ? Progress(S0, "Dynamic extinctions") : nothing

    step = 0

    while true
        step += 1

        alive_idx = findall(alive)

        if isempty(alive_idx)
            break
        end

        # --- build subnetwork for ranking ---
        A_sub = A[alive_idx, alive_idx]

        vuln = vec(sum(A_sub, dims=1))
        gen  = vec(sum(A_sub, dims=2))
        deg  = vuln .+ gen

        # --- SELECT SPECIES ---
        if criterion == :random_basal
            basal_mask = gen .== 0
            candidates = alive_idx[basal_mask]

            isempty(candidates) && break
            sp = rand(candidates)

        elseif criterion == :random_consumer
            consumer_mask = gen .> 0
            candidates = alive_idx[consumer_mask]

            isempty(candidates) && break
            sp = rand(candidates)

        else
            # --- metric-based selection ---
            metric = if criterion == :degree
                deg
            elseif criterion == :generality
                gen
            elseif criterion == :vulnerability
                vuln
            elseif criterion == :bodymass
                params.body_mass[alive_idx]
            elseif criterion == :random
                nothing
            else
                error("Unknown criterion: $criterion")
            end

            if metric === nothing
                local_choice = rand(1:length(alive_idx))
            else
                vals = metric
                target = descending ? maximum(vals) : minimum(vals)
                tied = findall(x -> x == target, vals)
                local_choice = rand(tied)
            end

            sp = alive_idx[local_choice]
        end

        # --- primary extinction ---
        alive[sp] = false

        # --- FIX: update cumulative primary counts ---
        push!(primary_counts, primary_counts[end] + 1)

        # --- build reduced system for simulation ---
        alive_idx = findall(alive)

        if isempty(alive_idx)
            break
        end

        A_sim = A[alive_idx, alive_idx]
        B_sim = B[alive_idx]
        # OG bodymass
        M_sim = params.M[alive_idx]

        # --- build food web without re-randomising ---
        fw = Foodweb(A_sim)

        # --- reuse original body masses ---
        params_sim = default_model(
            fw,
            BodyMass(M_sim),
            ClassicResponse(; h = 2),
        )

        # --- simulate ---
        try
            sol = simulate(params_sim, B_sim, t;
                show_degenerated = false,
                callback = CallbackSet(
                    extinction_callback(params_sim, survival_threshold)
                )
            )

            B_sim = sol.u[end]

        catch e
            @warn "Simulation failed" exception=e
            break
        end

        # --- write back ---
        B .= 0.0
        B[alive_idx] = B_sim

        # --- apply extinction threshold ---
        extinct_now = findall(x -> x < survival_threshold, B)
        alive[extinct_now] .= false
        B[extinct_now] .= 0.0

        # --- build output network ---
        survivors = findall(alive)

        if isempty(survivors)
            break
        end

        A_valid = A[survivors, survivors]
        push!(Nseq, build_network(A_valid))

        # --- progress ---
        if show_progress
            next!(prog)
        end

        if debug
            println("Step $step | remaining: $(length(survivors))")
        end
    end

    return (
        networks = Nseq,
        primary = primary_counts
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
        R[k] = robustness_integral(v.networks)
    end

    return R
end

function extinction_breakdown(result)

    Ns = result.networks
    primary_counts = result.primary

    n = min(length(Ns), length(primary_counts))

    Ns = Ns[1:n]
    primary_counts = primary_counts[1:n]

    # --- richness ---
    richness_vals = SpeciesInteractionNetworks.richness.(Ns)

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

function robustness_integral(network_sequence::Vector{<:SpeciesInteractionNetwork})
    
    isempty(network_sequence) && error("network_sequence cannot be empty")

    S = SpeciesInteractionNetworks.richness.(network_sequence)
    S0 = S[1]

    target = 0.5 * S0

    for i in 2:length(S)

        if S[i] <= target

            S1, S2 = S[i-1], S[i]

            # number of primary deletions at each step
            x1 = i - 2
            x2 = i - 1

            # linear interpolation between the two points
            if S1 == S2
                return x2 / S0
            end

            frac = (target - S2) / (S1 - S2)

            PDEL = (x2 + frac) / S0

            return PDEL
        end
    end

    # if 50% loss is never reached
    return 0.5
    
end