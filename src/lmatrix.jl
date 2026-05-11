using DataFrames
using SpeciesInteractionNetworks

"""
  lmatrix(species, bodymass, is_producer; ...)

    A port of the ATNr L-matrix function (Gauzens et al. 2023). 
    Determines links using an allometric Ricker function. 
    Returns a Unipartite Binary SpeciesInteractionNetwork.
"""
function lmatrix(
    species::Any,
    bodymass::Vector{Float64},
    is_producer::Vector{Bool};
    Ropt::Float64 = 100.0,    # Optimal consumer-prey body mass ratio
    γ::Float64 = 2.0,         # Width of the Ricker consumption window
    threshold::Float64 = 0.01, # Probabilistic cutoff for a "real" link
)

    S = length(species)
    prob_matrix = zeros(Bool, (S, S))

    for i = 1:S
        # The inner loop 'j' is the potential consumer
        for j = 1:S
            
            # 1. Biological Constraints
            # Only consumers can have outgoing links (predation), 
            # and we ignore producers eating other producers.
            if !is_producer[j]
                
                # 2. Calculate Body Mass Ratio (r)
                # How much bigger is the predator (j) than the prey (i)?
                r = bodymass[j] / bodymass[i]

                # 3. The Ricker Function
                # This models the "peak" of feeding efficiency.
                # If r == Ropt, the value is maximized at 1.0.
                # If r is too small (predator too tiny) or too large (prey too tiny), 
                # the value drops toward zero.
                
                # Formula: ((r / Ropt) * exp(1 - (r / Ropt)))^γ
                l_val = ((r / Ropt) * exp(1.0 - (r / Ropt)))^γ

                # 4. Binary Threshold
                # If the probability value exceeds the threshold, the link is formed.
                if l_val > threshold
                    prob_matrix[i, j] = true
                end
            end
        end
    end

    # Wrap the Boolean adjacency matrix into the network object
    return SpeciesInteractionNetwork(Unipartite(Symbol.(species)), Binary(prob_matrix))
end