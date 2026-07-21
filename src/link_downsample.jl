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
    _single_degree_product_downsample_step(mat)

Perform a single probabilistic pruning step based on the product of consumer
out-degree and resource in-degree.

P(i→j) ∝ kout(i) × kin(j)
"""
function _single_degree_product_downsample_step(mat::AbstractMatrix{Bool})

    S = size(mat,1)

    deg_out = vec(sum(mat,dims=2))
    deg_in  = vec(sum(mat,dims=1))

    prob_matrix = zeros(Float64,S,S)

    for i in 1:S

        prey = findall(mat[i,:])

        for j in prey
            prob_matrix[i,j] = deg_out[i] * deg_in[j]
        end

    end

    maxval = maximum(prob_matrix)

    if maxval > 0 && isfinite(maxval)
        prob_matrix ./= maxval
    else
        prob_matrix .= 0.0
    end

    prob_matrix = clamp.(prob_matrix,0.0,1.0)

    random_draw_matrix = rand(S,S) .<= prob_matrix

    return mat .& random_draw_matrix

end


# ==============================================================================
# Public API
# ==============================================================================

function downsample_degree_product_network(
    int_matrix::AbstractMatrix{Bool};
    target_co::Union{Nothing,Float64}=nothing,
    min_spp_prop::Float64=0.5,
    max_iter::Int=50
)

    S = size(int_matrix,1)

    if isnothing(target_co)
        return _single_degree_product_downsample_step(int_matrix)
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

        next_matrix = _single_degree_product_downsample_step(current_matrix)

        active_spp, next_co = _get_downsample_metrics(next_matrix,S)

        if active_spp < min_species
            @warn "Downsampling halted (iter $iter): Species retention threshold violated ($active_spp/$S active, limit is $min_species)."
            break
        end

        if sum(next_matrix) == 0
            @warn "Downsampling halted (iter $iter): Network collapsed to zero links."
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