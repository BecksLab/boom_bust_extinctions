using LinearAlgebra

# ==============================================================================
# Helper Functions
# ==============================================================================

function _get_downsample_metrics(mat::AbstractMatrix{Bool}, S::Int)
    active = sum(sum(mat, dims=1) .> 0 .|| sum(mat, dims=2)' .> 0)
    co = sum(mat) / (S^2)
    return active, co
end

"""
    _single_random_downsample_step(mat)

Perform a single probabilistic pruning step in which every existing interaction
has an equal probability of being retained.
"""
function _single_random_downsample_step(mat::AbstractMatrix{Bool})

    S = size(mat,1)

    prob_matrix = zeros(Float64,S,S)

    for i in 1:S
        prey = findall(mat[i,:])
        prob_matrix[i,prey] .= 1.0
    end

    random_draw_matrix = rand(S,S) .<= prob_matrix

    return mat .& random_draw_matrix
end


# ==============================================================================
# Public API
# ==============================================================================

function downsample_random_network(
    int_matrix::AbstractMatrix{Bool};
    target_co::Union{Nothing,Float64}=nothing,
    min_spp_prop::Float64=0.5,
    max_iter::Int=50
)

    S = size(int_matrix,1)

    if isnothing(target_co)
        return _single_random_downsample_step(int_matrix)
    end

    current_matrix = copy(int_matrix)

    min_species = ceil(Int,S*min_spp_prop)

    _, current_co = _get_downsample_metrics(current_matrix,S)

    if current_co <= target_co
        @warn "Initial connectance ($current_co) is already <= target ($target_co). Returning original network."
        return current_matrix
    end

    iter = 0

    while current_co > target_co && iter < max_iter

        iter += 1

        links = findall(current_matrix)

        if isempty(links)
            @warn "Network collapsed to zero links."
            break
        end

        # Remove exactly one random link
        link = rand(links)
        next_matrix = copy(current_matrix)
        next_matrix[link] = false

        active_spp, next_co = _get_downsample_metrics(next_matrix,S)

        if active_spp < min_species
            @warn "Downsampling halted (iter $iter): Species retention threshold violated."
            break
        end

        current_matrix = next_matrix
        current_co = next_co
    end

    return current_matrix

end