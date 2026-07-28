using EcologicalNetworksDynamics
using Random
using CairoMakie

# ============================================================
# SYSTEM SETUP
# ============================================================

S = 10
C = 0.1

Random.seed!(66)

foodweb = Foodweb(:niche; S, C)

# basal species = no incoming links
basal = findall(i -> sum(foodweb.A[i, :]) == 0, 1:S)

println("Basal species: ", basal)

# ============================================================
# FORCING FUNCTIONS
# ============================================================

supply_linear(t) = 8.0 + 0.05 * t

biomass_linear(t) = 0.5 + 0.01 * t

# ============================================================
# COMMON OUTPUT CONTAINER
# ============================================================

function initialize_outputs()

    out_t = Float64[]

    out_B = Matrix{Float64}(undef, 0, S)

    out_N = Float64[]

    return out_t, out_B, out_N

end

# ============================================================
# SCENARIO 1
# CONTINUOUS NUTRIENT FORCING
# ============================================================

function run_continuous_nutrient(; dt = 0.5, Tend = 100.0)

    Random.seed!(66)

    B = fill(0.5, S)

    N = [1.0]

    t = 0.0

    out_t, out_B, out_N = initialize_outputs()

    while t < Tend

        supp = supply_linear(t)

        nutrients = NutrientIntake(
            1;
            turnover = 0.2,
            supply = [supp]
        )

        m = default_model(
            foodweb,
            nutrients,
            Mortality(fill(0.6, S))
        )

        sol = simulate(m, B, dt; N0 = N)

        B = sol[1:S, end]

        N = sol[S+1, end]

        append!(out_t, sol.t .+ t)

        out_B = vcat(out_B, sol[1:S, :]')

        append!(out_N, sol[S+1, :])

        t += dt

    end

    return out_t, out_B, out_N

end

# ============================================================
# SCENARIO 2
# INITIAL NUTRIENT STATE ONLY
# ============================================================

function run_initial_nutrient(; Tend = 100.0)

    Random.seed!(66)

    B0 = fill(0.5, S)

    N0 = [supply_linear(0.0)]

    nutrients = NutrientIntake(
        1;
        turnover = 0.2,
        supply = [supply_linear(0.0)]
    )

    m = default_model(
        foodweb,
        nutrients,
        Mortality(fill(0.6, S))
    )

    sol = simulate(
        m,
        B0,
        Tend;
        N0 = N0
    )

    out_t = sol.t

    out_B = sol[1:S, :]'

    out_N = sol[S+1, :]

    return out_t, out_B, out_N

end

# ============================================================
# SCENARIO 3
# CONTINUOUS BASAL BIOMASS FORCING
# ============================================================

function run_continuous_biomass(; dt = 0.5, Tend = 100.0)

    Random.seed!(66)

    B = fill(0.5, S)

    t = 0.0

    out_t = Float64[]

    out_B = Matrix{Float64}(undef, 0, S)

    m = default_model(
        foodweb,
        Mortality(fill(0.6, S))
    )

    while t < Tend

        sol = simulate(
            m,
            B,
            dt
        )

        B = sol[:, end]

        imposed_biomass = biomass_linear(t + dt)

        B[basal] .= imposed_biomass

        append!(out_t, sol.t .+ t)

        out_B = vcat(out_B, Matrix(sol)')

        t += dt

    end

    return out_t, out_B

end

# ============================================================
# SCENARIO 4
# INITIAL BASAL BIOMASS ONLY
# ============================================================

function run_initial_biomass(; Tend = 100.0)

    Random.seed!(66)

    B0 = fill(0.5, S)

    B0[basal] .= biomass_linear(0.0)

    m = default_model(
        foodweb,
        Mortality(fill(0.6, S))
    )

    sol = simulate(
        m,
        B0,
        Tend
    )

    out_t = sol.t

    out_B = sol'

    return out_t, out_B

end

# ============================================================
# RUN ALL SCENARIOS
# ============================================================

results = Dict()

results["continuous_nutrient"] =
    run_continuous_nutrient()

results["initial_nutrient"] =
    run_initial_nutrient()

results["continuous_biomass"] =
    run_continuous_biomass()

results["initial_biomass"] =
    run_initial_biomass()

# ============================================================
# EXAMPLE PLOTS
# ============================================================

fig = Figure(size = (1200, 800))

ax1 = Axis(
    fig[1,1],
    title = "Continuous nutrient forcing",
    xlabel = "Time",
    ylabel = "Biomass"
)

t, B, N = results["continuous_nutrient"]

for sp in 1:S
    lines!(ax1, t, B[:, sp])
end

ax2 = Axis(
    fig[1,2],
    title = "Initial nutrient condition",
    xlabel = "Time",
    ylabel = "Biomass"
)

t, B, N = results["initial_nutrient"]

for sp in 1:S
    lines!(ax2, t, B[:, sp])
end

ax3 = Axis(
    fig[2,1],
    title = "Continuous biomass forcing",
    xlabel = "Time",
    ylabel = "Biomass"
)

t, B = results["continuous_biomass"]

for sp in 1:S
    lines!(ax3, t, B[:, sp])
end

ax4 = Axis(
    fig[2,2],
    title = "Initial biomass condition",
    xlabel = "Time",
    ylabel = "Biomass"
)

t, B = results["initial_biomass"]

for sp in 1:S
    lines!(ax4, t, B[:, sp])
end

display(fig)

save("../figures/four_scenario_comparison.png", fig)
