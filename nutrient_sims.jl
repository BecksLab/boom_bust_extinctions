using EcologicalNetworksDynamics
using Random
using CairoMakie

# ----------------------------
# MODEL SETUP (1 nutrient)
# ----------------------------
foodweb = Foodweb([0 0; 0 0])
mortality = Mortality([0.6, 1.2])

concentration = [
    1 => (1, 0.8),
    2 => (1, 1.1)
]

nutrients = NutrientIntake(
    1;
    concentration = concentration,
    turnover = 0.9,
    supply = [10.0],
    r = [1, 2],
)

m = default_model(foodweb, nutrients, mortality)

# ----------------------------
# FORCING FUNCTIONS
# ----------------------------
supply_const(t) = 10.0
supply_linear(t) = 8.0 + 0.05 * t
supply_seasonal(t) = 10.0 + 4.0 * sin(2π * t / 20.0)
supply_exp(t) = 8.0 * exp(0.01 * t)

# ----------------------------
# SIMULATION FUNCTION
# ----------------------------
function run_sim(m, supply_func; dt=0.5, Tend=100.0)

    Random.seed!(66)

    B0 = [0.5, 0.5]
    N0 = [1.0]

    t = 0.0
    B = copy(B0)
    N = copy(N0)

    sol_t = Float64[]
    sol_B1 = Float64[]
    sol_B2 = Float64[]
    sol_N = Float64[]

    while t < Tend

        nutrients = NutrientIntake(
            1;
            concentration = concentration,
            turnover = 0.9,
            supply = [supply_func(t)],
            r = [1, 2],
        )

        m_local = default_model(foodweb, nutrients, mortality)

        sol = simulate(m_local, B, dt; N0 = N)

        B = sol[1:2, end]
        N = sol[3, end]

        append!(sol_t, sol.t .+ t)
        append!(sol_B1, sol[1, :])
        append!(sol_B2, sol[2, :])
        append!(sol_N, sol[3, :])

        t += dt
    end

    return sol_t, sol_B1, sol_B2, sol_N
end

# ----------------------------
# RUN ALL SCENARIOS
# ----------------------------
t1, B11, B21, N1 = run_sim(m, supply_const)
t2, B12, B22, N2 = run_sim(m, supply_linear)
t3, B13, B23, N3 = run_sim(m, supply_seasonal)
t4, B14, B24, N4 = run_sim(m, supply_exp)

# ----------------------------
# PLOTTING
# ----------------------------
fig = Figure(; size = (900, 700))

labels = ["Constant", "Linear", "Seasonal", "Exponential"]

data = [
    (t1, B11, B21, N1),
    (t2, B12, B22, N2),
    (t3, B13, B23, N3),
    (t4, B14, B24, N4)
]

for i in 1:4
    t, B1, B2, N = data[i]

    ax = Axis(fig[i, 1],
        xlabel = i == 4 ? "Time" : "",
        ylabel = "Biomass / Nutrient",
        title = labels[i]
    )

    Makie.lines!(ax, t, B1, color = :red, label = "Plant 1")
    Makie.lines!(ax, t, B2, color = :green, label = "Plant 2")
    Makie.lines!(ax, t, N, color = :blue, linestyle = :dot, label = "Nutrient")

    axislegend(ax, position = :rt)
end

fig