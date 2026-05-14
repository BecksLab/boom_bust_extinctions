using EcologicalNetworksDynamics
using Random
using CairoMakie

# ----------------------------
# SYSTEM SETUP
# ----------------------------
S = 10
C = 0.1
foodweb = Foodweb(:niche; S, C)

# basal species (no incoming links)
basal = findall(i -> sum(foodweb.A[:, i]) == 0, 1:S)
println("Basal species: ", basal)

# ----------------------------
# FORCING FUNCTIONS
# ----------------------------
supply_const(t) = 10.0
supply_linear(t) = 8.0 + 0.05 * t
supply_seasonal(t) = 10.0 + 6.0 * sin(2π * t / 20.0)   # slightly stronger forcing
supply_exp(t) = 8.0 * exp(0.01 * t)

scenarios = [
    ("Constant", supply_const),
    ("Linear", supply_linear),
    ("Seasonal", supply_seasonal),
    ("Exponential", supply_exp)
]

# ----------------------------
# SIMULATION FUNCTION (UPDATED)
# ----------------------------
function run_sim(supply_func; dt=0.5, Tend=100.0)

    Random.seed!(66)

    B = fill(0.5, S)
    N = [1.0]

    t = 0.0

    out_t = Float64[]
    out_B = Matrix{Float64}(undef, 0, S)
    out_N = Float64[]

    while t < Tend

        supp = supply_func(t)

        # ----------------------------
        # KEY FIXES HERE
        # ----------------------------
        concentration = [
            i => (1, 0.8 + 0.4 * rand()) for i in 1:S
        ]

        nutrients = NutrientIntake(
            1;
            turnover = 0.2,          # ↓ slower damping → stronger forcing signal
            supply = [supp],
            concentration = concentration,
            r = fill(1.5, S)         # ↑ faster biomass response
        )

        m = default_model(foodweb, nutrients, Mortality(fill(0.6, S)))

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

# ----------------------------
# RUN ALL SCENARIOS
# ----------------------------
results = Dict()

for (name, f) in scenarios
    results[name] = run_sim(f)
end

# ----------------------------
# PLOTTING
# ----------------------------
fig = Figure(resolution = (1100, 900))

for (i, (name, _)) in enumerate(scenarios)

    t, Bmat, N = results[name]

    ax = Axis(fig[i, 1],
        title = name,
        xlabel = i == 4 ? "Time" : "",
        ylabel = "Density"
    )

    # basal species only (clear signal)
    for sp in basal
        lines!(ax, t, Bmat[:, sp], label = "Basal $sp")
    end

    lines!(ax, t, N, color = :blue, linestyle = :dash, label = "Nutrient")

    axislegend(ax, position = :rt)
end

fig