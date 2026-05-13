# ============================================================
# ECOLOGICAL NETWORK + DYNAMIC NUTRIENT FORCING FRAMEWORK
# ============================================================
# Main workflow:
# 1. Generate food web
# 2. Define nutrient forcing scenario (time-dependent environment)
# 3. Simulate ecological dynamics
# 4. Track biomass + extinctions
# 5. Analyse and visualise outputs
# ============================================================

# ---------------------------
# 1. Load dependencies
# ---------------------------

using CSV
using DataFrames
using DifferentialEquations
using Distributions
using EcologicalNetworksDynamics
using Extinctions
using JLD2
using Plots
using ProgressMeter
using SpeciesInteractionNetworks
using Statistics

# ---------------------------
# 2. Model parameters
# ---------------------------

S = 10                  # number of species
C = 0.1                 # connectance
tmax = 5000             # simulation time
survival_threshold = 1e-12

# ---------------------------
# 3. Build food web
# ---------------------------

fw = Foodweb(:niche; S, C)

# Baseline nutrient + ecological model
# NOTE: supply here is ONLY baseline (overrided it dynamically later)
params = default_model(
    fw,
    BodyMass(; Z = 10),
    NutrientIntake(1; supply = 10.2)
)

# Initial biomass distribution
B0 = rand(Uniform(0.1, 1), S)

# 4. NUTRIENT FORCING FRAMEWORK

# Abstract type for clarity (useful for thesis explanation)
abstract type NutrientScenario end

struct Constant <: NutrientScenario end
struct Linear <: NutrientScenario end
struct Exponential <: NutrientScenario end
struct Sigmoid <: NutrientScenario end
struct Seasonal <: NutrientScenario end

# ------------------------------------------------------------
# Time-dependent nutrient supply function
#
# S0 = baseline nutrient level
# t  = current simulation time
# scenario controls how supply is modified over time
# ------------------------------------------------------------

function supply(t, S0, scenario, tmax)
    if scenario isa Constant
        # No forcing: supply stays at baseline
        return S0
    elseif scenario isa Linear
        # Linear decline from baseline to 20% baseline over the full run
        return max(0.2 * S0, S0 * (1 - 0.8 * t / tmax))
    elseif scenario isa Exponential
        # Exponential decay toward a nonzero minimum
        rate = 0.001
        return S0 * exp(-rate * t) + 0.2 * S0
    elseif scenario isa Sigmoid
        # Smooth logistic transition from high to low supply
        midpoint = tmax / 2
        steepness = 0.02
        min_supply = 0.2 * S0
        return min_supply + (S0 - min_supply) / (1 + exp(steepness * (t - midpoint)))
    elseif scenario isa Seasonal
        # Periodic seasonal variation around baseline supply
        amplitude = 0.3 * S0
        period = 250.0
        return max(0.0, S0 + amplitude * sin(2π * t / period))
    else
        error("Unknown nutrient scenario: $(typeof(scenario))")
    end
end

function run_simulation(m, B0, tmax, S0, scenario)

    # ---------------------------
    # CONSTANT CASE (no forcing)
    # ---------------------------
    if scenario isa Constant

        return simulate(
            m,
            B0,
            tmax;
            N0 = [S0],
            show_degenerated = false
        )
    end

    # ---------------------------
    # DYNAMIC CASE (forcing)
    # ---------------------------
    function update_nutrients!(integrator)
        t = integrator.t
        model = integrator.p.model
        model.nutrients.supply .= supply(t, S0, scenario, tmax)
    end

    cb = PeriodicCallback(update_nutrients!, 1.0)

    return simulate(
        m,
        B0,
        tmax;
        N0 = [S0],
        callback = CallbackSet(
            extinction_callback(m, survival_threshold),
            cb
        ),
        show_degenerated = false
    )
end

# 6. EXPERIMENT SETUP

S = 10
C = 0.1
tmax = 5000
survival_threshold = 1e-12

fw = Foodweb(:niche; S, C)

m = default_model(
    fw,
    BodyMass(; Z = 10),
    NutrientIntake(1; supply = 10.2)
)

B0 = rand(Uniform(0.1, 1), S)


scenarios = [
    Constant(),
    Linear(),
    Exponential(),
    Sigmoid(),
    Seasonal()
]

scenario_names = [
    "Constant",
    "Linear",
    "Exponential",
    "Sigmoid",
    "Seasonal"
]

results = Dict()

for (sc, name) in zip(scenarios, scenario_names)

    println("Running: ", name)

    sol = run_simulation(m, B0, tmax, 10.2, sc)

    results[name] = (sol, sc)
end

biomass = Dict()
supply_history = Dict()

for (name, (sol, sc)) in results

    biomass[name] = hcat(sol.u...)
    supply_history[name] = [supply(t, 10.2, sc, tmax) for t in sol.t]
end

total_biomass = Dict()

for (name, B) in biomass

    total_biomass[name] = vec(sum(B, dims=1))
end

# 7. POST-PROCESSING PLACEHOLDER

burnin = 100  # number of initial timesteps to discard as burn-in

plot()

for name in scenario_names
    plot!(
        total_biomass[name],
        label = name
    )
end

xlabel!("Time")
ylabel!("Total biomass")
title!("Biomass dynamics under nutrient forcing scenarios")
xlims!(burnin, tmax)
ylims!(10, 12.5)

# Plot the actual nutrient supply trajectories for each scenario.
# This makes the forcing explicit and shows the modified supply inputs.
plot_supply = plot()
for name in scenario_names
    plot!(plot_supply, results[name][1].t, supply_history[name], label = name)
end
xlabel!(plot_supply, "Time")
ylabel!(plot_supply, "Nutrient supply")
title!(plot_supply, "Nutrient supply forcing scenarios")

display(plot_supply)

println("Simulation complete.")