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

# ---------------------------
# 1. Load dependencies
# ---------------------------

using GMT
using EcologicalNetworksDynamics
using Random
using CairoMakie

# ---------------------------
# 2. Create synthetic nutrient scenarios
# ---------------------------

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

# ---------------------------
# 3. Model parameters
# ---------------------------

S = 10      # number of species
C = 0.1     # connectance
tmax = length(nutrient_data.time)  # simulation time based on data length
survival_threshold = 1e-12

# ---------------------------
# 4. Build food web
# ---------------------------

fw = Foodweb(:niche; S, C)

function run_sim(nutrient_val; dt=1, Tend=100.0)

    Random.seed!(66)

    B = fill(0.5, S)
    N = [1.0]

    t = 1

    out_t = Float64[]
    out_B = Matrix{Float64}(undef, 0, S)
    out_N = Float64[]

    while t < Tend

        concentration = [
            i => (1, 0.2 + 0.4 * rand()) for i in 1:length(basal)
        ]

        nutrients = NutrientIntake(
            1;
            turnover = 0.2,
            supply = nutrient_val[t],
            concentration = concentration,
            r = [i => 1 for i in basal]
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
# RUN
# ----------------------------

run_sim(nutrient_data.nutrient1)