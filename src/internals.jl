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

    results["degree_high"] = extinction(N, "degree", true)
    results["degree_low"]  = extinction(N, "degree", false)

    results["vul_high"] = extinction(N, "vulnerability", true)
    results["vul_low"]  = extinction(N, "vulnerability", false)

    results["gen_high"] = extinction(N, "generality", true)
    results["gen_low"]  = extinction(N, "generality", false)

    results["bm_high"] = extinction(N, Symbol.(sortperm(params.body_mass, rev=true)))
    results["bm_low"]  = extinction(N, Symbol.(sortperm(params.body_mass, rev=false)))

    results["rand_basal"]    = extinction(N; protect = :consumer)
    results["rand_consumer"] = extinction(N; protect = :basal)

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
    show_progress = true)

    B = copy(B_init)
    A = params.A

    S = length(B)
    removed = falses(S)

    Nseq = SpeciesInteractionNetwork[]

    # --- initial network ---
    survivors_init, A_init = extract_valid_network(A, B, survival_threshold)
    if survivors_init !== nothing
        push!(Nseq, build_network(A_init))
    end

    step = 0
    max_steps = S

    prog = show_progress ? Progress(max_steps, "Dynamic extinctions ($(criterion))") : nothing

    while true

        step += 1

        if step > max_steps
            @warn "Max steps reached — stopping early"
            break
        end

        # --- Identify alive species ---
        alive = findall(i -> (B[i] > survival_threshold) && !removed[i], eachindex(B))

        if isempty(alive)
            break
        end

        # --- Subnetwork ---
        A_sub = A[alive, alive]

        # --- Compute vulnerability ---
        vulnerability = vec(sum(A_sub, dims=1))

        candidates = Int[]
        metric = nothing

        if criterion == :degree
            metric = vec(sum(A_sub, dims=2)) + vec(sum(A_sub, dims=1))
        elseif criterion == :generality
            metric = vec(sum(A_sub, dims=2))
        elseif criterion == :vulnerability
            metric = vec(sum(A_sub, dims=1))
        elseif criterion == :bodymass
            metric = params.body_mass[alive]
        elseif criterion == :random_basal
            candidates = findall(x -> x == 0, vulnerability)
        elseif criterion == :random_consumer
            candidates = findall(x -> x > 0, vulnerability)
        else
            error("Unknown criterion")
        end

        # --- Select species ---
        if criterion in (:random_basal, :random_consumer)

            if isempty(candidates)
                break
            end

            sp = alive[rand(candidates)]

        else
            if metric === nothing || isempty(metric)
                break
            end

            target_val = descending ? maximum(metric) : minimum(metric)
            candidates = findall(x -> x == target_val, metric)

            if isempty(candidates)
                break
            end

            sp = alive[rand(candidates)]
        end

        # --- Force extinction ---
        B[sp] = 0.0
        removed[sp] = true

        # --- Simulation ---
        try
            sol = simulate(params, B, t;
                callback = CallbackSet(
                    extinction_callback(params, survival_threshold)
                )
            )

            B = sol.u[end]

            # prevent regrowth
            B[removed] .= 0.0

        catch e
            @warn "Simulation failed — stopping early" exception=e
            break
        end

        # --- Extract valid network (connectivity filter) ---
        survivors_after, A_valid = extract_valid_network(A, B, survival_threshold)

        if survivors_after === nothing
            break
        end

        push!(Nseq, build_network(A_valid))

        # --- Debug (optional, very useful first run) ---
        # println("Step $step: alive = ", length(survivors_after))

        # --- Progress ---
        if show_progress
            next!(prog; showvalues = [
                (:step, step),
                (:alive, length(survivors_after))
            ])
        end
    end

    return Nseq
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

    for (k, seq) in results_dict
        R[k] = robustness(seq)
    end

    return R
end