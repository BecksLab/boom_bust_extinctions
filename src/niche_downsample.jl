using LinearAlgebra

# ==============================================================================
# Helper Functions
# ==============================================================================

# Helper to calculate active species and current connectance (same as original)
function _get_downsample_metrics(mat::AbstractMatrix{Bool}, S::Int)
    active = sum(sum(mat, dims=1) .> 0 .|| sum(mat, dims=2)' .> 0)
    co = sum(mat) / (S^2)
    return active, co
end

"""
    _get_svd_niche_positions(mat::AbstractMatrix{Bool}) -> Vector{Float64}

Extract a 1D niche dimension from the binary adjacency matrix using the first 
singular vectors of Singular Value Decomposition (SVD).
"""
function _get_svd_niche_positions(mat::AbstractMatrix{Bool})
    S = size(mat, 1)
    
    # Convert binary matrix to Float64 for SVD calculation
    mat_float = Float64.(mat)
    
    # Run Singular Value Decomposition: A = U * S * Vt
    # We only need the first singular vector (the primary axis of variance)
    svd_decomp = svd(mat_float)
    u1 = svd_decomp.U[:, 1]
    v1 = svd_decomp.V[:, 1]
    s1 = svd_decomp.S[1]
    
    # Scale coordinates by the strength of the primary singular value
    x_consumer = u1 .* sqrt(s1)
    x_resource = v1 .* sqrt(s1)
    
    # Align to a single unified niche coordinate (average of consumer and resource states)
    n_raw = (x_consumer .+ x_resource) ./ 2.0
    
    # Min-max normalize raw values to the [0.0, 1.0] interval
    min_val = minimum(n_raw)
    max_val = maximum(n_raw)
    
    if max_val > min_val
        n_axis = (n_raw .- min_val) ./ (max_val - min_val)
    else
        # Fallback if SVD collapses or is perfectly uniform
        n_axis = collect(range(0.0, 1.0, length=S))
    end
    
    return n_axis
end

# Single-step probabilistic niche-space prune
"""
    _single_niche_downsample_step(mat::AbstractMatrix{Bool}, sigma_scale::Float64) -> Matrix{Bool}

Perform a single-step probabilistic link pruning on the network based on niche distances.

# Arguments
- `mat::AbstractMatrix{Bool}`: Binary adjacency matrix (Rows = Consumers, Cols = Resources).
- `sigma_scale::Float64`: Strictness parameter. Controls the relative dietary niche breadth. 
  Smaller values compress diet widths, forcing strict specialist-like behavior.
"""
function _single_niche_downsample_step(mat::AbstractMatrix{Bool}, sigma_scale::Float64)
    S = size(mat, 1)
    
    # 1. Extract the latent 1D niche positions using SVD
    n_axis = _get_svd_niche_positions(mat)
    
    # 2. Compute each consumer's diet centroid on the niche axis
    centroids = zeros(Float64, S)
    for i in 1:S
        prey_indices = findall(mat[i, :])
        if !isempty(prey_indices)
            centroids[i] = mean(n_axis[prey_indices])
        else
            centroids[i] = n_axis[i] # Fallback to own position if no prey
        end
    end
    
    # 3. Calculate distance-based probability matrix
    # P(i -> j) = exp( - (d_ij / sigma_i)^2 )
    prob_matrix = zeros(Float64, S, S)
    
    for i in 1:S
        prey_indices = findall(mat[i, :])
        if isempty(prey_indices)
            continue
        end
        
        # Calculate standard deviation of current prey niche positions to act as baseline niche breadth
        prey_positions = n_axis[prey_indices]
        base_sigma = length(prey_positions) > 1 ? std(prey_positions) : 0.1
        if base_sigma == 0.0
            base_sigma = 0.1 # Prevent divide-by-zero for single-prey specialists
        end
        
        # Apply strictness scaling factor
        sigma_i = base_sigma * sigma_scale
        
        for j in prey_indices
            dist = abs(n_axis[j] - centroids[i])
            prob_matrix[i, j] = exp(-(dist / sigma_i)^2)
        end
    end
    
    # Normalize and Clamp
    maxval = maximum(prob_matrix)
    if maxval > 0 && isfinite(maxval)
        prob_matrix ./= maxval
    else
        prob_matrix .= 0.0
    end
    prob_matrix = clamp.(prob_matrix, 0.0, 1.0)

    # Probabilistic draw
    random_draw_matrix = rand(S, S) .<= prob_matrix

    return random_draw_matrix
end

# ==============================================================================
# Public API Function
# ==============================================================================

"""
    downsample_niche_network(int_matrix::AbstractMatrix{Bool}, sigma_scale::Float64; kwargs...) -> Matrix{Bool}

Downsample a food web's interaction matrix based on topological niche-space constraints.

# Arguments
- `int_matrix::AbstractMatrix{Bool}`: Binary adjacency matrix (size `S × S`).
- `sigma_scale::Float64`: Controls the "niche strictness" (typically around 0.5 to 2.0). 

# Keyword Arguments
- `target_co::Union{Nothing, Float64} = nothing`: Target connectance. 
- `min_spp_prop::Float64 = 0.5`: Defensive species conservation limit.
- `max_iter::Int = 50`: Defensive safeguard loop break.
"""
function downsample_niche_network(
    int_matrix::AbstractMatrix{Bool}, 
    sigma_scale::Float64;
    target_co::Union{Nothing, Float64} = nothing,
    min_spp_prop::Float64 = 0.5,
    max_iter::Int = 50
)
    S = size(int_matrix, 1)
    
    # --- Case 1: Single-Step Downsampling ---
    if isnothing(target_co)
        return _single_niche_downsample_step(int_matrix, sigma_scale)
    end

    # --- Case 2: Iterative Connectance-Targeted Downsampling ---
    current_matrix = copy(int_matrix)
    min_species = ceil(Int, S * min_spp_prop)
    _, current_co = _get_downsample_metrics(current_matrix, S)

    if current_co <= target_co
        @warn "Initial connectance ($current_co) is already <= target ($target_co). Returning original network"
        return current_matrix
    end

    iter = 0
    while current_co > target_co && iter < max_iter
        iter += 1
        # In iterative mode, SVD coordinates & centroids are recalculated dynamically 
        # on the newly updated matrix during each step.
        next_matrix = _single_niche_downsample_step(current_matrix, sigma_scale)
        active_spp, next_co = _get_downsample_metrics(next_matrix, S)

        if active_spp < min_species
            @warn "Downsampling halted (iter $iter): Species retention threshold violated ($active_spp/$S active, limit is $min_species)."
            break
        end

        if sum(next_matrix) == 0
            @warn "Downsampling halted (iter $iter): Network collapsed to 0 links."
            break
        end

        current_matrix = next_matrix
        current_co = next_co
    end

    if iter == max_iter && current_co > target_co
        @warn "Reached max iterations ($max_iter) without hitting target connectance. Current Co: $current_co"
    end

    return current_matrix
end