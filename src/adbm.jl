# For implementing the ADBM model

"""
  adbm_parameters(df::DataFrame,
    bodymass::Vector{Float64};...)

  Returns the parameters needed for the adbm model. Defaults to the values
  specified in BioEnergeticFoodWebs.jl.

"""
function adbm_parameters(
    df::DataFrame,
    bodymass::Vector{Float64};
    e::Float64 = 1.0,
    a_adbm::Float64 = 0.0189,
    ai::Float64 = -0.491,
    aj::Float64 = -0.465,
    b::Float64 = 0.401,
    h_adbm::Float64 = 1.0,
    hi::Float64 = 1.0,
    hj::Float64 = 1.0,
    n::Float64 = 1.0,
    ni::Float64 = -0.75,
    Hmethod::Symbol = :ratio,
    Nmethod::Symbol = :original,
)

    parameters = Dict{Symbol,Any}(
        :e => e,
        :a_adbm => a_adbm,
        :ai => ai,
        :aj => aj,
        :b => b,
        :h_adbm => h_adbm,
        :hi => hi,
        :hj => hj,
        :n => n,
        :ni => ni,
    )

    #check Hmethod
    if Hmethod ∈ [:ratio, :power]
        parameters[:Hmethod] = Hmethod
    else
        error("Invalid value for Hmethod -- must be :ratio or :power")
    end

    # check Nmethod
    if Nmethod ∈ [:original, :biomass]
        parameters[:Nmethod] = Nmethod
    else
        error("Invalid value for Nmethod -- must be :original or :biomass")
    end

    # add empty cost matrix
    S = size(df, 1)
    parameters[:costMat] = ones(Float64, (S, S))

    # Identify producers - based on tiering class
    is_producer = df.tiering .== "primary"
    parameters[:is_producer] = is_producer

    # add bodymass
    parameters[:bodymass] = bodymass

    return parameters
end


"""
  _get_adbm_terms(S::Int64, parameters::Dict{Symbol,Any}, biomass::Vector{Float64})

  This function takes the parameters for the ADBM model and returns the final
  terms used to determine feeding patterns. It is used internally by adbmmodel().

"""
function _get_adbm_terms(S::Int64, parameters::Dict{Symbol,Any}, biomass::Vector{Float64})
    E = parameters[:e] .* parameters[:bodymass]
    if parameters[:Nmethod] == :original
        N = parameters[:n] .* (parameters[:bodymass] .^ parameters[:ni])
    elseif parameters[:Nmethod] == :biomass
        N = biomass
    end
    A_adbm =
        parameters[:a_adbm] *
        (parameters[:bodymass] .^ parameters[:aj]) *
        (parameters[:bodymass] .^ parameters[:ai])' # a * pred * prey
    for i = 1:S #for each prey
        A_adbm[:, i] = A_adbm[:, i] .* N[i]
    end
    λ = A_adbm
    if parameters[:Hmethod] == :ratio
        H = zeros(Float64, (S, S))
        ratios = (parameters[:bodymass] ./ parameters[:bodymass]')' #PREDS IN ROWS : PREY IN COLS
        for i = 1:S, j = 1:S
            if ratios[j, i] < parameters[:b]
                H[j, i] = parameters[:h_adbm] / (parameters[:b] - ratios[j, i])
            else
                H[j, i] = Inf
            end
        end
    elseif parameters[:Hmethod] == :power
        H =
            parameters[:h_adbm] *
            (parameters[:bodymass] .^ parameters[:hj]) *
            (parameters[:bodymass] .^ parameters[:hi])' # h * pred * prey
    end

    adbmTerms = Dict{Symbol,Any}(:E => E, :λ => λ, :H => H)

    return adbmTerms
end

"""
  _get_feeding_links(S::Int64,E::Vector{Float64}, λ::Array{Float64},
    H::Array{Float64},biomass::Vector{Float64},j)

  This function takes the terms calculated by _get_adbm_terms() and uses them to 
  determine the feeding links of species j. Used internally by adbmmodel().

"""
function _get_feeding_links(
    E::Vector{Float64},
    λ::Array{Float64},
    H::Array{Float64},
    biomass::Vector{Float64},
    j,
)

    profit = E ./ H[j, :]
    # Setting profit of species with zero biomass  to -1.0
    # This prevents them being included in the profitSort
    profit[vec(biomass .== 0.0)] .= -1.0

    profs = sortperm(profit, rev = true) # returns the index

    λSort = λ[j, profs]
    HSort = H[j, profs]
    ESort = E[profs]

    λH = cumsum(λSort .* HSort)
    Eλ = cumsum(ESort .* λSort)

    λH[isnan.(λH)] .= Inf
    Eλ[isnan.(Eλ)] .= Inf

    cumulativeProfit = Eλ ./ (1 .+ λH)

    if all(0 .== cumulativeProfit)
        feeding = []
    else
        feeding = profs[1:maximum(findall(cumulativeProfit .== maximum(cumulativeProfit)))]
    end

    #cumulativeProfit[end] = NaN
    #feeding = profs[(append!([true],cumulativeProfit[1:end-1] .< profitSort[2:end]))]
    return feeding
end


"""
adbmmodel(df::DataFrame, parameters::Dict{Symbol,Any}, biomass::Vector{Float64})

  This function returns the food web based on the ADBM model of Petchey et al. 2008.
  The function takes the paramteres created by adbm_parameters and uses 
  _get_adbm_terms and _get_feeding_links to determine the web structure.

  Note this (and all internal) functions has been ported from the 
  BioEnergeticFoodWebs.jl source code and has been (minimally) modified for the 
  purpose of this project.

  #### References

  Petchey, Owen L., Andrew P. Beckerman, Jens O. Riede, and Philip H. Warren.
  2008. “Size, Foraging, and Food Web Structure.” Proceedings of the National
  Academy of Sciences 105 (11): 4191–96. https://doi.org/10.1073/pnas.0710672105.

"""
function adbmmodel(df::DataFrame, parameters::Dict{Symbol,Any}, biomass::Vector{Float64})

    # data checks
    if "species" ∉ names(df)
        error("Missing $(v) variable as a column in DataFrame, add or rename")
    end

    S = nrow(df)

    adbmMAT = zeros(Bool, (S, S))
    adbmTerms = _get_adbm_terms(S, parameters, biomass)
    E = adbmTerms[:E]
    λ = adbmTerms[:λ]
    H = adbmTerms[:H]
    for j = 1:S
        if !parameters[:is_producer][j]
            if biomass[j] > 0.0
                feeding = _get_feeding_links(E, λ, H, biomass, j)
                adbmMAT[j, feeding] .= 1
            end
        end
    end


    edges = Binary(adbmMAT)
    nodes = Unipartite(Symbol.(df.species)) # return actual species name metadata
    N = SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)
    return N
end