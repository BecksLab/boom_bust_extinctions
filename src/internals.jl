"""
    build_network(adj_mat::Matrix{Bool})

Construct a `SpeciesInteractionNetwork` from a Boolean adjacency matrix.

Converts a square adjacency matrix into a unipartite, binary food web where
species are assigned generic symbolic identifiers (`:1`, `:2`, ..., `:S`).

# Arguments
- `adj_mat::Matrix{Bool}`: Square adjacency matrix (S × S), where:
    - `adj_mat[i, j] = true` indicates a trophic interaction (i consumes j)
    - Rows and columns refer to the same ordered species set

# Returns
- `SpeciesInteractionNetwork{Unipartite{Symbol}, Binary{Bool}}`

# Notes
- Species identities are not preserved (relabelled as `:1`, ..., `:S`)
- Interaction direction is preserved
- Assumes adjacency matrix is already ecologically valid

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

Extract an ecologically valid subnetwork based on biomass and connectivity.

Species are retained if they:
1. Have biomass above a survival threshold
2. Participate in at least one trophic interaction

# Arguments
- `A::Matrix`: Full adjacency matrix
- `B::Vector`: Species biomasses
- `threshold::Float64`: Extinction threshold

# Returns
- `(survivors, A_sub)`:
    - `survivors::Vector{Int}`: Indices of retained species
    - `A_sub::Matrix`: Induced subnetwork adjacency matrix
- `(nothing, nothing)` if no valid species remain

# Ecological Interpretation
- Species with `B ≤ threshold` are considered extinct
- Species with no prey AND no consumers are removed as non-functional
- Ensures resulting network is trophically connected

# Notes
- Mirrors assumptions used in topological extinction workflows
- Does NOT distinguish basal vs consumer roles explicitly
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

Run topological (structural) extinction sequences on a food web.

Species are removed sequentially according to structural traits, and
secondary extinctions occur via trophic disconnection (cascade mechanism).

# Arguments
- `N`: `SpeciesInteractionNetwork`
- `params`: Model parameters (used for body mass ordering)

# Returns
- `Dict{String, NamedTuple}`:
    - `networks`: sequence of networks after each extinction step
    - `primary`: cumulative number of primary deletions

# Extinction Scenarios
- `"degree_high"` / `"degree_low"`: connectivity-based removal
- `"vul_high"` / `"vul_low"`: vulnerability (in-degree)
- `"gen_high"` / `"gen_low"`: generality (out-degree)
- `"bm_high"` / `"bm_low"`: body mass ordering
- `"rand_basal"`: random removal of basal species only
- `"rand_consumer"`: random removal of consumers only

# Ecological Assumptions
- Secondary extinctions occur when species lose all prey
- Basal species are defined by generality = 0 (in initial network)
- Random removal respects initial trophic roles (fixed classification)

# Important Implementation Detail
- `remove_disconnected = false`:
    - Disconnected species are NOT removed unless they lose all prey
    - This prevents artificial inflation of secondary extinctions

# Notes
- No population dynamics are simulated
- Fully comparable to classical topological robustness studies
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
    dynamic_extinction_adaptive(params, B0; kwargs...)

Run adaptive dynamical extinction simulations with sequential species removal.

At each step:
1. Select a target species (based on a criterion)
2. Force its extinction (primary removal)
3. Simulate system dynamics for `t` time steps
4. Record secondary extinctions (biomass collapse)
5. Repeat until system collapse

# Arguments
- `params`: Model parameters (includes adjacency matrix and body masses)
- `B0::Vector`: Initial biomasses after burn-in

# Keyword Arguments
- `criterion::Symbol`: Removal rule
    - `:degree`, `:generality`, `:vulnerability`, `:bodymass`
    - `:random_basal`, `:random_consumer`
- `descending::Bool`: Remove highest (true) or lowest (false)
- `t::Int`: Simulation duration after each deletion
- `survival_threshold::Float64`: Biomass extinction threshold
- `show_progress::Bool`
- `debug::Bool`

# Returns
- NamedTuple:
    - `networks`: sequence of surviving subnetworks
    - `primary`: cumulative number of primary deletions

# Ecological Interpretation
- Combines press perturbations with relaxation dynamics
- Secondary extinctions emerge from population dynamics
- Captures indirect effects and trophic cascades

# Key Assumptions
- Body masses are FIXED from initial network (no re-sampling)
- Basal/consumer roles for random deletions are FIXED from initial network
- Species go extinct if `B < survival_threshold`

# Notes
- Multiple species may go extinct per step
- Extinction order is adaptive (recomputed each step)
- This corresponds to the dynamical approach in the Curtsdottir 2011
"""
function dynamic_extinction_adaptive(params, B0;
    criterion = :degree,
    descending = true,
    t = 5000,
    survival_threshold = 1e-12,
    show_progress = true,
    debug = false
)

    S0 = length(B0)

    A = params.A
    B = copy(B0)

    # --- species roles ---
    gen0 = vec(sum(A, dims=2))
    basal0 = gen0 .== 0
    consumer0 = gen0 .> 0

    alive = trues(S0)

    # --- store networks ---
    Nseq = SpeciesInteractionNetwork[]
    push!(Nseq, build_network(Matrix(A)))  # include initial state

    primary_counts = Int[0]

    prog = show_progress ? Progress(S0, "Dynamic extinctions") : nothing
    step = 0

    # --- stopping condition ---
    function continue_condition()
        if criterion == :random_basal
            return any(alive .& basal0)
        elseif criterion == :random_consumer
            return any(alive .& consumer0)
        else
            # all trait-based = consumers only
            return any(alive .& consumer0)
        end
    end

    while continue_condition()
        step += 1

        alive_idx = findall(alive)
        isempty(alive_idx) && break

        # --- build subnetwork for ranking ---
        A_sub = A[alive_idx, alive_idx]

        vuln = vec(sum(A_sub, dims=1))
        gen  = vec(sum(A_sub, dims=2))
        deg  = vuln .+ gen

        # --- SELECT SPECIES ---
        if criterion == :random_basal
            candidates = filter(i -> alive[i] && basal0[i], eachindex(alive))

        elseif criterion == :random_consumer
            candidates = filter(i -> alive[i] && consumer0[i], eachindex(alive))

        else
            # trait-based → consumers only
            candidates = filter(i -> alive[i] && consumer0[i], eachindex(alive))

            # metric evaluated on alive_idx
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

            if !isempty(candidates)
                # map candidates to local indices
                local_inds = findall(i -> alive_idx[i] in candidates, eachindex(alive_idx))

                if metric === nothing
                    local_choice = rand(local_inds)
                else
                    vals = metric[local_inds]
                    target = descending ? maximum(vals) : minimum(vals)
                    tied = local_inds[findall(x -> x == target, vals)]
                    local_choice = rand(tied)
                end

                sp = alive_idx[local_choice]
            else
                break
            end
        end

        isempty(candidates) && break
        sp = rand(candidates)

        # --- primary extinction ---
        alive[sp] = false
        push!(primary_counts, primary_counts[end] + 1)

        # --- simulate ALL remaining species ---
        alive_idx = findall(alive)
        isempty(alive_idx) && break

        A_sim = A[alive_idx, alive_idx]
        B_sim = B[alive_idx]
        M_sim = params.M[alive_idx]

        fw = Foodweb(A_sim)

        params_sim = default_model(
            fw,
            BodyMass(M_sim),
            ClassicResponse(; 
                            h = params.h,
                            #c = params.intraspecific_interference[1],
                            ),
        )

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

        # --- record network ---
        survivors = findall(alive)
        isempty(survivors) && break

        A_valid = A[survivors, survivors]
        push!(Nseq, build_network(A_valid))

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

"""
    extinction_breakdown(result)

Decompose extinction trajectories into primary and secondary components.

# Arguments
- `result`: NamedTuple with:
    - `networks`: sequence of networks
    - `primary`: cumulative primary deletions

# Returns
- NamedTuple:
    - `primary`: proportion of primary removals (PDEL)
    - `secondary`: proportion of secondary extinctions
    - `total`: total extinction proportion

# Definitions
- Primary extinction: directly removed species
- Secondary extinction: species lost due to cascading effects
- Total extinction = primary + secondary

# Notes
- All values are normalised by initial richness (S0)
- Primary is constrained to not exceed total extinction
- Ensures physically consistent decomposition

# Warnings
- Emits warnings if extinction curve is non-monotonic
"""
function extinction_breakdown(result)

    Ns = result.networks
    primary_counts = result.primary

    if isempty(Ns)
        return NaN   
    end

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

        if !hasproperty(c, :total)
        
            rows = DataFrame(
                net_id = Int[],
                type = String[],
                scenario = String[],
                step = Union{Missing, Int}[],
                primary = Union{Missing, Float64}[],
                secondary = Union{Missing, Float64}[],
                total = Union{Missing, Float64}[]
            )
        
        else
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
    end

    return rows
end

"""
    robustness_integral(network_sequence)

Compute robustness (R50) from an extinction sequence.

R50 is defined as the proportion of species that must be removed
(primary deletions) to cause 50% total species loss.

# Arguments
- `network_sequence`: Vector of `SpeciesInteractionNetwork`

# Returns
- `Float64`: R50 value (∈ [0, 0.5])

# Method
- Tracks species richness decline across extinction steps
- Identifies when richness falls below 50% of initial richness (S0)
- Uses linear interpolation between steps to estimate R50

# Interpretation
- R50 = 0.5 → no secondary extinctions (max robustness)
- R50 < 0.5 → cascading extinctions reduce robustness

# Notes
- Assumes 1 primary deletion per step
- Uses step index as proxy for PDEL (as in topological approaches)
"""
function robustness_integral(network_sequence::Vector{<:SpeciesInteractionNetwork})
    
    if isempty(network_sequence)
        return NaN   
    end

    S = SpeciesInteractionNetworks.richness.(network_sequence)
    S0 = S[1]

    target = 0.5 * S0

    for i in 2:length(S)

        if S[i] <= target

            S1, S2 = S[i-1], S[i]

            # PDEL values
            x1 = (i - 2) / S0
            x2 = (i - 1) / S0

            if S1 == S2
                return x2
            end

            # linear interpolation
            frac = (target - S1) / (S2 - S1)

            R50 = x1 + frac * (x2 - x1)

            return R50
        end
    end

    return (length(S) - 1) / S0
    
end

function robustness_auc(network_sequence)
    S = SpeciesInteractionNetworks.richness.(network_sequence)
    S0 = S[1]

    frac = S ./ S0
    x = (0:length(S)-1) ./ S0

    return sum((frac[1:end-1] .+ frac[2:end]) ./ 2 .* diff(x))
end

function compute_trophic_levels(A)
    
    S = size(A, 1) # Species richness.
    out_degree = sum(A; dims = 2)
    D = -(A ./ out_degree) # Diet matrix.
    D[isnan.(D)] .= 0.0
    D[diagind(D)] .= 1.0 .- D[diagind(D)]
    # Solve with the inverse matrix.
    inverse = iszero(det(D)) ? pinv : inv
    tls = inverse(D) * ones(S)

    return tls
end