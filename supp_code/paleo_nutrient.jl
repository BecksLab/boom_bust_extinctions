#=
Main workflow:
- Generate PFIM reference networks
- Use starting and end EcoGenie points to infer nutrient curve
- Assign bodymass
- Add plankton groups with initial biomass
- Infer biomass for other species by scaling up from plankton biomass
- Burn-in PFIM networks
- Topological extinctions
- Dynamic extinctions
=#

# --- 1. General Set-up ---

using CSV
using DataFrames
using DifferentialEquations
using Distributions
using EcologicalNetworksDynamics
using Extinctions
using JLD2
using pfim
using ProgressMeter
using SpeciesInteractionNetworks
using Statistics
using CairoMakie

include("../src/internals.jl")

import Random
Random.seed!(66)

# --- storage ---
rows = Dict[]
topo_curve_store = DataFrame()
dyn_curve_store = DataFrame()
species_store = DataFrame()

# --- data ---
traits = CSV.read("../data/community_dolomites.csv", DataFrame)
feeding_rules = CSV.read("../data/feeding_rules.csv", DataFrame)

# --- global params ---
t = 5000
survival_threshold = 1e-12

# --- body size distribution ---
global_dist = LogNormal(log(30), 1.5)

size_bounds = Dict(
    "primary" => (0.01, 0.1),
    "tiny" => (0.1, 10.0),
    "small" => (10.0, 50.0),
    "medium" => (50.0, 100.0),
    "large" => (100.0, 300.0),
    "very_large" => (300.0, 500.0),
    "gigantic" => (500.0, Inf)
)

# --- Body sizes ---

y = collect(String, traits.size)

bodysize = [
    begin
        lo, hi = size_bounds[s]
        rand(truncated(global_dist, lo, hi))
    end
    for s in y
]

traits[!, :bodymass] = bodysize
traits.biomass = fill(missing, nrow(traits))

# import plankton df (note this also has ecogenie biomass)
plankton = CSV.read("../data/community_plankton.csv", DataFrame)

df = vcat(traits, plankton)

# BIOMASS

known = .!ismissing.(df.biomass)

b = -3/4

a = exp(mean(log.(df.biomass[known]) .- b .* log.(df.bodymass[known])))

predicted = df.bodymass .^ b

df.biomass[.!known] .= predicted[.!known]

biomass = float.(df.biomass)

select!(df, Not(:biomass, :bodymass))

# --- PFIM networks ---

pfim_down = PFIM(
    df,
    feeding_rules;
    return_type=:matrix,
    y=3.0,
    downsample=true
)

# create food web object
# need to transpose because BEFW has cons as col and resource as row
foodweb = Foodweb(Matrix(transpose(pfim_down)))

S = size(pfim_down, 1)

# basal species (no incoming links)
basal = findall(i -> sum(foodweb.A[i, :]) == 0, 1:S)
println("Basal species: ", basal)

# ----------------------------
# FORCING FUNCTIONS
# ----------------------------

# Inital surface [O2]
poc_pre  = 0.5153117
poc_post = 0.7147954

T = t

supply_linear(t) = poc_pre + ((poc_post-poc_pre)/T) * t
supply_exp(t) = poc_pre + (poc_post - poc_pre) * (t /T)

scenarios = [
    ("Linear", supply_linear),
    ("Exponential", supply_exp)
]

# ----------------------------
# SIMULATION FUNCTION
# ----------------------------
function run_sim(supply_func; dt=0.5, Tend=100.0)

    Random.seed!(66)

    B = biomass
    N = [1.0]

    t = 0.0

    out_t = Float64[]
    out_B = Matrix{Float64}(undef, 0, S)
    out_N = Float64[]

    while t < Tend

        supp = supply_func(t)

        concentration = [
            i => (1, 0.2 + 0.4 * rand()) for i in 1:length(basal)
        ]

        nutrients = NutrientIntake(
            1;
            #turnover=0.2,
            supply=[supp],
            #concentration=concentration,
            r=[i => 1 for i in basal]
        )

        m = default_model(foodweb, nutrients, Mortality(fill(0.6, S)))

        sol = simulate(m, B, dt; N0=N)

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
fig = Figure(; size=(1100, 900))

for (i, (name, _)) in enumerate(scenarios)

    t, Bmat, N = results[name]

    ax = Axis(fig[i, 1],
        title=name,
        xlabel=i == 4 ? "Time" : "",
        ylabel="Density"
    )

    # basal species only (clear signal)
    for sp in 1:S
        Makie.lines!(ax, t, Bmat[:, sp], label="Spp $sp")
    end

    Makie.lines!(ax, t, N, color=:blue, linestyle=:dash, label="Nutrient")

    axislegend(ax, position=:rt)
end

fig

save("../figures/paleo_nutrient_scenarios.png", fig)