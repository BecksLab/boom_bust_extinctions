# General sundry internal functions

using FoodWebTools
using LinearAlgebra
using SpeciesInteractionNetworks
using Statistics


"""
_network_summary(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

    Returns the 'summary statistics' for a network
"""
function _network_summary(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})
    A = _get_matrix(N)
    L = links(N)
    S = SpeciesInteractionNetworks.richness(N)

    _gen = SpeciesInteractionNetworks.generality(N)
    gen = collect(values(_gen))
    _vul = SpeciesInteractionNetworks.vulnerability(N)
    vul = collect(values(_vul))
    ind_maxgen = findmax(gen)[2]
    l_s = L / S
    top = sum(vec(sum(A, dims = 1) .== 0))
    basal = sum(vec(sum(A, dims = 2) .== 0))
    int = (S - (basal + top))
    tl = trophic_level(A)

    chain = chain_metrics(A; max_depth=6)

    # for centrality - freeman centralisation
    c = collect(values(centrality(N)))
    Cmax = maximum(c)


    D = Dict{Symbol,Any}(
        :richness => S,
        :links => L,
        :connectance => connectance(N),
        :distance => distancetobase(N, collect(keys(_gen))[ind_maxgen]),
        :basal => basal / S,
        :top => top / S,
        :intermediate => int / S,
        :herbivory => length(herbivore(N)) / S,
        :omnivory => length(omnivore(N)) / S,
        :cannibal => length(cannibal(N)) / S,
        :predpreyRatio => (basal + int)/(top + int),
        :l_S => l_s,
        :GenSD => std(gen) / l_s,
        :VulSD => std(vul) / l_s,
        :TL => mean(collect(values(tl))),
        :ChLen => chain.ChLen,
        :ChSD  => chain.ChSD,
        :ChNum => chain.ChNum,
        :path => mean(pathlengths(N)),
        :LinkSD => std(values(SpeciesInteractionNetworks.degree(N))) / l_s,
        :S1 => length(findmotif(motifs(Unipartite, 3)[1], remove_cannibals(N)))/(SpeciesInteractionNetworks.richness(N)^2),
        :S2 => length(findmotif(motifs(Unipartite, 3)[2], remove_cannibals(N)))/(SpeciesInteractionNetworks.richness(N)^2),
        :S4 => length(findmotif(motifs(Unipartite, 3)[4], remove_cannibals(N)))/(SpeciesInteractionNetworks.richness(N)^2),
        :S5 => length(findmotif(motifs(Unipartite, 3)[5], remove_cannibals(N)))/(SpeciesInteractionNetworks.richness(N)^2),
        :centrality => sum(Cmax .- c),
        :loops => length(loops(N)) / S,
        :intervals => intervality(A),
        :MaxSim => max_sim(A),
        :Clust => clustering(A),
        :trophicCoherence => trophic_coherence(A),
        :trophicVar => trophic_variance(A),
    )
    return D
end

"""
_get_matrix(N::SpeciesInteractionNetwork)

    Internal function to return a matrix of interactions from a
    SpeciesInteractionNetwork
"""
function _get_matrix(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    species = SpeciesInteractionNetworks.richness(N)
    n = zeros(Bool, (species, species))
    for i in axes(n, 1)
        for j in axes(n, 2)
            if N.edges[i, j] == true
                n[i, j] = true
            end
        end
    end

    return n
end

"""
pathlengths(N::SpeciesInteractionNetwork)

    Returns the shortest pathlengths between all species pairs for a network
"""
function pathlengths(N::SpeciesInteractionNetwork)

    sp = species(N)
    path = Any[]

    for i in eachindex(sp)

        _path = collect(values(shortestpath(N, sp[i])))

        if length(_path) > 0
            append!(path, _path)
        end
    end

    return path
end

"""
herbivore(N::SpeciesInteractionNetwork)

    Returns a vector of species that are herbivores (only consume basal species)
"""
function herbivore(N::SpeciesInteractionNetwork)

    # find basal species
    gen = SpeciesInteractionNetworks.generality(N)
    basal = collect(keys(filter(((k, v),) -> v == 0, gen)))

    sp = species(N)
    herbivores = Any[]

    for i in eachindex(sp)

        prey = collect(successors(N, sp[i]))

        # is the prey a subset of or equal to
        if length(prey) > 0 && prey ⊆ basal
            push!(herbivores, sp[i])
        end
    end

    return herbivores
end

"""
omnivore(N::SpeciesInteractionNetwork)

    Returns a vector of species that are omnivores (feed on  > 1 species of different 
    trophic levels)
"""
function omnivore(N::SpeciesInteractionNetwork)

    omni = Any[]

    A = N.edges.edges
    sp = species(N)

    tl = trophic_level(A; species=sp)

    for i in eachindex(sp)
        prey = collect(successors(N, sp[i]))

        # return trophic level of prey
        _tls = [v for (k, v) in tl if k ∈ prey]

        if length(prey) > 1 && !allequal(_tls)
            push!(omni, sp[i])
        end
    end

    return omni
end

"""
cannibal(N::SpeciesInteractionNetwork)

    Returns a vector of species that are cannibals (feed on themselves)
"""
function cannibal(N::SpeciesInteractionNetwork)

    _cannibal = Any[]
    sp = species(N)

    for i in eachindex(sp)
        prey = collect(successors(N, sp[i]))

        if sp[i] ∈ prey
            push!(_cannibal, sp[i])
        end
    end

    return _cannibal
end

"""
Returns the percentage of species involved in a loop (motifs S3, D2, D4, D5, D6, D7)
"""
function loops(N::SpeciesInteractionNetwork)

    A = _get_matrix(N)

    # empty vector to puch spp index (proxy for id) to
    in_loop_spp = Any[]

    # look at loops that expand all the way to richness of network
    for i in 3:SpeciesInteractionNetworks.richness(N)
        # get the diagonal of the power transformed adj mat (indication of loops)
        _diag = diag(A^i)
        # get indices of spp with val > 0
        spp_ind = findall(!=(0), _diag)

        # add to master list
        if length(spp_ind) > 0
            append!(in_loop_spp, spp_ind)
        end
    end

    # return unique numbers (indices)
    return unique(in_loop_spp)
end

"""
remove_cannibals(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    Identifies and sets cannibalistic link to zero
"""
function remove_cannibals(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    A = _get_matrix(N)

    A[diagind(A)] .= false

    nodes = Unipartite(species(N))
    edges = Binary(A)
    network = SpeciesInteractionNetwork(nodes, edges)

    return network
end

"""
trophic_variance(A::AbstractMatrix{Bool})

Returns variance in trophic levels across species.
"""
function trophic_variance(A::AbstractMatrix{Bool})

    tl = trophic_level(A)
    tls = collect(values(tl))

    return var(tls)

end