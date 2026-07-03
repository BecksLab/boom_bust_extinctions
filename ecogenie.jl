# ============================================================
# ECOLOGICAL NETWORK + PALEO NUTRIENT FORCING FRAMEWORK
# ============================================================
# Main workflow:
# 1. Load paleo nutrient data from netCDF
# 2. Build food web
# 3. Simulate ecological dynamics with real nutrient forcing
# 4. Track biomass + extinctions
# 5. Analyse and visualise outputs
# ============================================================

# --- 1. Load dependencies ---

using GMT
using EcologicalNetworksDynamics
using Random
using CairoMakie

# --- 2. Create synthetic nutrient scenarios ---

# Load paleo nutrient data (kept for future use)
ecogenie_raw = grdinterpolate("data/ecogenie.nc"; track = (103, 17));

# Create synthetic nutrient scenarios for 5000 timesteps
n_timesteps = 5000
time_vec = 1:n_timesteps

# Scenario 1: Constant nutrient supply
nutrient1_constant = fill(10.0, n_timesteps)

# Scenario 2: Linear increase from 5 to 15
nutrient1_linear = 5.0 .+ (10.0 / n_timesteps) .* time_vec

# Create DataFrame with both scenarios
nutrient_data = DataFrame(
    time = time_vec, 
    nutrient1 = nutrient1_constant, 
    nutrient2 = nutrient1_linear
)

# --- 3. Model parameters ---

S = 10      # number of species
C = 0.1     # connectance
tmax = length(nutrient_data.time)  # simulation time based on data length
survival_threshold = 1e-12

# --- 4. Build food web ---

fw = Foodweb(:niche; S, C)

function run_sim(nutrient1, nutrient2; dt=1, Tend=5000)

    Random.seed!(66)

    # initial biomasses
    B = fill(0.5, S)

    # TWO nutrients now
    N = [1.0, 1.0]

    t = 1

    out_t = Float64[]
    out_B = Matrix{Float64}(undef, 0, S)

    out_N1 = Float64[]
    out_N2 = Float64[]

    # identify basal species
    basal = findall(i -> sum(fw.A[i, :]) == 0, 1:S)

    # create fixed nutrient geometry
    #concentration = [
    #    i => (1, 0.2 + 0.4 * rand()) for i in 1:length(basal)
    #]
#
    #half_saturation = [
    #    i => (1, 0.2 + 0.4 * rand()) for i in 1:length(basal)
    #]

    while t < Tend

        # TWO nutrient supplies

        nutrients = NutrientIntake(
            2;

            # one turnover value per nutrient
            turnover = [0.2, 0.2],

            # supply for BOTH nutrients
            supply = [
                nutrient1[t],
                nutrient2[t]
            ],

            # species × nutrients matrices
            #concentration = concentration,
            #half_saturation = half_saturation,

            # basal species indices
            #r = [i => 1 for i in basal]
        )

        m = default_model(
            fw,
            nutrients,
            Mortality(fill(0.6, S))
        )

        # simulate short interval
        sol = simulate(m, B, dt; N0 = N)

        # update state variables
        B = sol[1:S, end]

        # now TWO nutrient states
        N = sol[S+1:S+2, end]

        # store outputs

        append!(out_t, sol.t .+ t)

        out_B = vcat(out_B, sol[1:S, :]')

        append!(out_N1, sol[S+1, :])
        append!(out_N2, sol[S+2, :])

        t += dt
    end

    return out_t, out_B, out_N1, out_N2
end

# --- 5. Run ---

t, B, N1, N2 = run_sim(
    nutrient_data.nutrient1,
    nutrient_data.nutrient2
)

# --- 6. Plot ---

fig = Figure()

ax = Axis(
    fig[1,1],
    xlabel = "Time",
    ylabel = "Value",
    title = "Dual nutrient forcing"
)

# species biomass
for i in 1:S
    Makie.lines!(ax, t, B[:, i])
end

# nutrients
Makie.lines!(ax, t, N1,
    color = :black,
    linestyle = :dash,
    linewidth = 3,
    label = "Nutrient 1"
)

Makie.lines!(ax, t, N2,
    color = :blue,
    linestyle = :dot,
    linewidth = 3,
    label = "Nutrient 2"
)

axislegend(ax)

save("figures/multi_nutrients.png", fig)