# ==============================================================================
# CCR Competition Chemostat Definitions
# 
# Core definitions for pairwise chemostat competition simulations among
# Strategy N (no CCR, strategy_idx=0), Strategy T (transcriptional CCR, strategy_idx=1),
# and Strategy P (post-transcriptional CCR, strategy_idx=2).
# Supports any two-strategy pair: [0,1], [0,2], or [1,2].
#
# This file provides:
#   - Parameter struct and defaults
#   - S1 pulsed inflow function
#   - ODE system (competition_model!) — accelerating mRNA version (active)
#   - ODE system (competition_model!) — original version (commented out)
#   - Initial condition setup (from file or defaults)
#   - Volume fraction grid generation
#   - Linear regression for net growth rate estimation
#   - Ne*s neutral evolution analysis
#   - Competition outcome classification
# ==============================================================================

using DifferentialEquations
using LinearAlgebra
using Statistics
using CSV
using DataFrames
using Random
using Printf
using GLM
using DiffEqCallbacks


# Cache for initial conditions read from file (avoids repeated I/O)
const INITIAL_CONDITIONS_CACHE = Dict{String, Tuple{Dict{Int, Vector{Float64}}, Dict{Int, Vector{Float64}}}}()


# ==============================================================================
# Parameter Struct
# ==============================================================================

mutable struct Parameters
    # --- Substrate uptake ---
    k1::Float64          # S1 maximum specific uptake rate
    k2::Float64          # S2 maximum specific uptake rate
    KS1::Float64         # S1 uptake half-saturation constant
    KS2::Float64         # S2 uptake half-saturation constant

    # --- Metabolism ---
    alpha1::Float64      # X1 metabolic rate constant
    alpha2::Float64      # X2 metabolic rate constant

    # --- Gene expression ---
    kM::Float64          # mRNA synthesis rate constant
    kE::Float64          # Enzyme synthesis rate constant (translational)
    delta1::Float64      # mRNA degradation rate constant
    delta2::Float64      # Enzyme degradation rate constant
    KX1::Float64         # Induction half-saturation constant for M1
    KX2::Float64         # Induction half-saturation constant for M2
    b::Float64           # Basal expression fraction

    # --- CCR parameters ---
    K_CCR1::Float64      # Repression half-saturation constant (Transcriptional, Strategy T)
    K_CCR2::Float64      # Repression half-saturation constant (Post-transcriptional, Strategy P)

    # --- Yield coefficients ---
    Y1::Float64          # Yield coefficient for S1 pathway
    Y2::Float64          # Yield coefficient for S2 pathway

    # --- Operational parameters ---
    D::Float64           # Dilution rate
    Sr1::Float64         # S1 reservoir concentration
    Sr2::Float64         # S2 reservoir concentration (S2_in)

    # --- Wave parameters for S1 feeding ---
    epsilon::Float64     # S1 feeding amplitude
    k_wave::Float64      # Steepness parameter for S1 pulse function
    T1::Float64          # S1 ON duration
    T2::Float64          # S1 OFF duration
end


"""
    Parameters()

Construct default parameter set (Table 2 in Model section).
Dilution rate D can be overridden via the environment variable `CCR_D_VALUE`.
"""
function Parameters()
    D_value = haskey(ENV, "CCR_D_VALUE") ? parse(Float64, ENV["CCR_D_VALUE"]) : 0.04

    return Parameters(
        # Substrate uptake
        1.0,    # k1
        1.0,    # k2
        0.5,    # KS1
        0.5,    # KS2

        # Metabolism
        2.0,    # alpha1
        1.0,    # alpha2

        # Gene expression
        0.2,    # kM
        20.0,   # kE
        0.5,    # delta1
        0.0,    # delta2
        0.01,   # KX1
        0.01,   # KX2
        0.001,  # b

        # CCR parameters
        0.01,   # K_CCR1
        0.01,   # K_CCR2

        # Yield coefficients
        0.7,    # Y1
        0.07,   # Y2

        # Operational parameters
        D_value,# D
        0.0,    # Sr1
        2.0,    # Sr2

        # Wave parameters
        2.0,    # epsilon
        10.0,   # k_wave
        30.0,   # T1
        70.0    # T2
    )
end


# ==============================================================================
# Substrate Inflow Functions
# ==============================================================================

"""
    S1_in(t, params) -> Float64

Preferred carbon source inflow concentration at time `t` (Eq. 4-5 in Model).

Uses a smoothed logistic pulse to avoid numerical discontinuities:
- T1 > 0, T2 = 0  →  constant supply at epsilon
- T1 = 0           →  no supply
- T1 > 0, T2 > 0   →  periodic pulse with period T1 + T2
"""
function S1_in(t, params::Parameters)
    if params.T1 <= 0
        return 0.0
    end

    T0 = params.T1 + params.T2
    if T0 <= 0
        return 0.0
    end

    # Constant supply when no OFF phase
    if params.T1 >= T0
        return params.epsilon
    end

    # Periodic pulse
    t_in_cycle = t % T0

    if 0 < t_in_cycle < params.T1
        t_rise_center = 0.0
        t_fall_center = params.T1
        rise_curve = 1.0 / (1.0 + exp(-params.k_wave * (t_in_cycle - t_rise_center)))
        fall_curve = 1.0 / (1.0 + exp(-params.k_wave * (t_in_cycle - t_fall_center)))
        offset_rise = 1.0 / (1.0 + exp(-params.k_wave * (0.0 - t_rise_center)))
        offset_fall = 1.0 / (1.0 + exp(-params.k_wave * (0.0 - t_fall_center)))
        corrected_shape = (rise_curve - fall_curve) - (offset_rise - offset_fall)
        return params.epsilon * 2.0 * corrected_shape
    else
        return 0.0
    end
end


"""S2 inflow concentration (constant, Eq. 6 in Model)."""
function S2_in(t, params::Parameters)
    return params.Sr2
end


# ==============================================================================
# ODE System: Pairwise Chemostat Competition
#
# State layout:
#   u[1]  = S1  (shared preferred substrate)
#   u[2]  = S2  (shared non-preferred substrate)
#   Per strategy i (offset = 2 + (i-1)*8):
#     u[offset+1..7] = X1, M1, E1, X2, M2, E2, c_minus  (intracellular)
#     u[offset+8]    = c  (biomass concentration)
#
# Strategy indices:
#   0 = Strategy N (No CCR)
#   1 = Strategy T (Transcriptional CCR)
#   2 = Strategy P (Post-transcriptional CCR)
# ==============================================================================

# ------------------------------------------------------------------------------
# Original version (without mRNA acceleration)
# To activate: uncomment this block and comment out the accelerated version below
# ------------------------------------------------------------------------------

function competition_model!(du, u, p, t)
    params, strategy_indices = p

    S1 = u[1]
    S2 = u[2]

    # Substrate dynamics: inflow/dilution (consumption added per strategy below)
    dS1dt = (S1_in(t, params) - S1) * params.D
    dS2dt = (S2_in(t, params) - S2) * params.D

    for (i, strategy_idx) in enumerate(strategy_indices)
        offset = 2 + (i - 1) * 8
        X1      = u[offset + 1]
        M1      = u[offset + 2]
        E1      = u[offset + 3]
        X2      = u[offset + 4]
        M2      = u[offset + 5]
        E2      = u[offset + 6]
        c_minus = u[offset + 7]
        c       = u[offset + 8]

        # Growth rate (Eq. 20)
        rg = (params.k1 * E1 * S1 / (params.KS1 + S1) +
              params.k2 * E2 * S2 / (params.KS2 + S2)) -
             ((1 - params.Y1) * params.alpha1 * X1 +
              (1 - params.Y2) * params.alpha2 * X2)

        # CCR regulation functions (Eq. 12-13, 16-17)
        #   kappa_M: transcriptional CCR modifier for M2 synthesis
        #   kappa_E: post-transcriptional CCR modifier for E2 synthesis
        if strategy_idx == 0       # Strategy N: No CCR
            kappa_M = 1.0
            kappa_E = 1.0
        elseif strategy_idx == 1   # Strategy T: Transcriptional CCR (Eq. 13)
            kappa_M = params.K_CCR1 / (params.K_CCR1 + X1)
            kappa_E = 1.0
        elseif strategy_idx == 2   # Strategy P: Post-transcriptional CCR (Eq. 17)
            kappa_M = 1.0
            kappa_E = params.K_CCR2 / (params.K_CCR2 + X1)
        else
            error("Unknown strategy_idx: $strategy_idx. Supported: 0 (N), 1 (T), 2 (P).")
        end

        # Substrate consumption by this strategy
        dS1dt -= params.k1 * E1 * S1 / (params.KS1 + S1) * c
        dS2dt -= params.k2 * E2 * S2 / (params.KS2 + S2) * c

        # ODEs
        # dX1/dt (Eq. 8)
        du[offset + 1] = params.k1 * E1 * S1 / (params.KS1 + S1) -
                          params.alpha1 * X1 - rg * X1
        # dM1/dt (Eq. 10)
        du[offset + 2] = params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * c_minus -
                          params.delta1 * M1 - rg * M1
        # dE1/dt (Eq. 14)
        du[offset + 3] = params.kE * M1 * c_minus -
                          params.delta2 * E1 - rg * E1
        # dX2/dt (Eq. 9)
        du[offset + 4] = params.k2 * E2 * S2 / (params.KS2 + S2) -
                          params.alpha2 * X2 - rg * X2
        # dM2/dt (Eq. 11) — includes kappa_M for transcriptional CCR
        du[offset + 5] = params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * c_minus * kappa_M -
                          params.delta1 * M2 - rg * M2
        # dE2/dt (Eq. 15) — includes kappa_E for post-transcriptional CCR
        du[offset + 6] = params.kE * M2 * c_minus * kappa_E -
                          params.delta2 * E2 - rg * E2
        # dc_minus/dt (Eq. 18)
        du[offset + 7] = (params.Y1 * params.alpha1 * X1 + params.Y2 * params.alpha2 * X2) +
                          (params.delta1 * M1 + params.delta2 * E1 +
                           params.delta1 * M2 + params.delta2 * E2) -
                          (params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * c_minus +
                           params.kE * M1 * c_minus +
                           params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * c_minus * kappa_M +
                           params.kE * M2 * c_minus * kappa_E) -
                          rg * c_minus
        # dc/dt (Eq. 21)
        du[offset + 8] = (rg - params.D) * c

        # Non-negativity clipping
        for j in (offset+1):(offset+8)
            if u[j] <= 0.0 && du[j] < 0.0
                du[j] = 0.0
            end
        end
    end

    # Non-negativity for substrates
    if u[1] <= 0.0 && dS1dt < 0.0; dS1dt = 0.0; end
    if u[2] <= 0.0 && dS2dt < 0.0; dS2dt = 0.0; end

    du[1] = dS1dt
    du[2] = dS2dt

    return nothing
end



# ==============================================================================
# Linear Regression for Net Growth Rate 
# ==============================================================================

"""
    calculate_net_growth_rates_by_regression(sol, strategies_needing_analysis, 
        strategy_indices, params; ...) -> (net_growth_rates, r_squared_values, regression_info)

Estimate net growth rates via log-linear regression on cycle-endpoint biomass.

Uses biological cycle endpoints (at multiples of T1+T2) when periodic,
or artificial intervals as fallback.
"""
function calculate_net_growth_rates_by_regression(sol, strategies_needing_analysis, strategy_indices, params;
                                                 sampling_interval=nothing,
                                                 n_cycles=100,
                                                 verbose=false)
    # Determine sampling interval
    if sampling_interval === nothing
        cycle_length = params.T1 + params.T2
        if cycle_length > 0
            sampling_interval = cycle_length
            sampling_type = "biological_cycles"
            if verbose
                println("Using biological cycle endpoints (interval = $(cycle_length))")
            end
        else
            sampling_interval = 1.0
            sampling_type = "artificial_intervals"
            if verbose
                println("Using artificial sampling intervals (interval = $(sampling_interval))")
            end
        end
    else
        sampling_type = "user_defined"
        if verbose
            println("Using user-defined sampling interval (interval = $(sampling_interval))")
        end
    end

    # Initialise result containers
    net_growth_rates = Dict{Int, Float64}()
    r_squared_values = Dict{Int, Float64}()
    regression_info  = Dict{Int, Dict{String, Any}}()

    for i in strategies_needing_analysis
        strategy_idx = strategy_indices[i]
        net_growth_rates[strategy_idx] = NaN
        r_squared_values[strategy_idx] = NaN
    end

    # Generate sampling times
    if sampling_type == "biological_cycles"
        cycle_length = sampling_interval
        first_cycle = floor(Int, sol.t[1] / cycle_length)
        last_cycle  = floor(Int, sol.t[end] / cycle_length)
        sampling_times = [(first_cycle + c) * cycle_length
                          for c in 1:(last_cycle - first_cycle)]
        if verbose
            println("Analysis window: [$(sol.t[1]), $(sol.t[end])]")
            println("Cycles $(first_cycle)..$(last_cycle), $(length(sampling_times)) endpoints")
        end
    else
        total_artificial_cycles = floor(Int, (sol.t[end] - sol.t[1]) / sampling_interval)
        sampling_times = [sol.t[1] + c * sampling_interval for c in 0:total_artificial_cycles]
        if verbose
            println("Using $(length(sampling_times)) artificial cycle endpoints")
        end
    end


    # --- Regression per strategy ---
    for analysis_idx in strategies_needing_analysis
        strategy_idx = strategy_indices[analysis_idx]

        if verbose
            println("\n=== ANALYZING STRATEGY $(strategy_idx) ===")
        end

        biomass_offset = 2 + (analysis_idx - 1) * 8 + 8

        regression_info[strategy_idx] = Dict{String, Any}(
            "n_cycles_used"     => length(sampling_times),
            "sampling_times"    => Float64[],
            "sampling_biomass"  => Float64[],
            "log_biomass"       => Float64[],
            "intercept"         => NaN,
            "slope"             => NaN,
            "r_squared"         => NaN,
            "p_value"           => NaN,
            "glm_success"       => false,
            "analysis_method"   => sampling_type,
            "sampling_interval" => sampling_interval
        )

        # Extract biomass at sampling points
        sampling_biomass = Float64[]
        valid_times      = Float64[]

        for (pt, t_sample) in enumerate(sampling_times)
            closest_idx = argmin(abs.(sol.t .- t_sample))
            raw_biomass = sol[biomass_offset, closest_idx]

            if verbose && pt <= 5
                println("Sample $(pt): t=$(round(t_sample, digits=2)), biomass=$(raw_biomass)")
            end

            if raw_biomass > 0
                push!(sampling_biomass, raw_biomass)
                push!(valid_times, t_sample)
            end
        end

        if length(sampling_biomass) < 3
            if verbose
                println("ERROR: Insufficient data points ($(length(sampling_biomass)) < 3)")
            end
            continue
        end

        regression_info[strategy_idx]["sampling_times"]  = valid_times
        regression_info[strategy_idx]["sampling_biomass"] = sampling_biomass

        if verbose
            println("Data: $(length(sampling_biomass)) points, " *
                    "range [$(minimum(sampling_biomass)), $(maximum(sampling_biomass))]")
        end

        # Log-linear regression via GLM
        try
            y_log = log.(sampling_biomass)
            regression_info[strategy_idx]["log_biomass"] = y_log

            try
                relative_times = valid_times .- valid_times[1]
                df = DataFrame(time=relative_times, log_biomass=y_log)
                model = lm(@formula(log_biomass ~ time), df)

                slope     = coef(model)[2]
                intercept = coef(model)[1]
                r_sq      = r2(model)
                p_val     = try; pvalue(model); catch; NaN; end

                net_growth_rates[strategy_idx] = slope
                r_squared_values[strategy_idx] = r_sq

                regression_info[strategy_idx]["intercept"]  = intercept
                regression_info[strategy_idx]["slope"]       = slope
                regression_info[strategy_idx]["r_squared"]   = r_sq
                regression_info[strategy_idx]["p_value"]     = p_val
                regression_info[strategy_idx]["glm_success"] = true

                result_type = isfinite(slope) ? "finite" :
                              isnan(slope)    ? "nan"    :
                              isinf(slope)    ? "infinite" : "other"
                regression_info[strategy_idx]["result_type"] = result_type

                if verbose
                    println("GLM: slope=$(slope), R²=$(r_sq), type=$(result_type)")
                end

            catch glm_error
                if verbose
                    println("GLM failed: $(glm_error)")
                end
                regression_info[strategy_idx]["glm_success"] = false
                regression_info[strategy_idx]["glm_error"]   = string(glm_error)
            end

        catch outer_error
            if verbose
                println("Log-transform error: $(outer_error)")
            end
            regression_info[strategy_idx]["glm_success"] = false
            regression_info[strategy_idx]["setup_error"]  = string(outer_error)
        end
    end

    if verbose
        println("\n=== REGRESSION SUMMARY ===")
        for analysis_idx in strategies_needing_analysis
            strategy_idx = strategy_indices[analysis_idx]
            if haskey(net_growth_rates, strategy_idx) && !isnan(net_growth_rates[strategy_idx])
                s = net_growth_rates[strategy_idx]
                if isfinite(s)
                    println("Strategy $(strategy_idx): slope=$(s) ($(s > 0 ? "GROWING" : "DECLINING"))")
                else
                    println("Strategy $(strategy_idx): slope=$(s) (non-finite)")
                end
            else
                println("Strategy $(strategy_idx): GLM failed")
            end
        end
    end

    return net_growth_rates, r_squared_values, regression_info
end


# ==============================================================================
# Simulation Runner
# ==============================================================================

"""
    run_grid_simulation(params, strategy_indices, u0; Ne, target_tau, analysis_cycles)

Run a single competition simulation and return the solution restricted 
to the final `analysis_cycles` complete cycles.
"""
function run_grid_simulation(params::Parameters, strategy_indices, u0::Vector{Float64};
                             Ne::Float64=1e8,
                             target_tau::Float64=2e8,
                             analysis_cycles::Int=50)
    # Convert generations to physical time (Eq. 27: t = tau * ln(2) / D)
    t_target = (target_tau * log(2)) / params.D
    cycle_length = params.T1 + params.T2

    # Determine analysis time window
    local t_eval
    if cycle_length <= 0
        analysis_start = max(0.0, t_target - analysis_cycles)
        t_eval = collect(analysis_start:1.0:t_target)
    else
        total_cycles = floor(Int, t_target / cycle_length)
        n_analyze    = min(analysis_cycles, total_cycles)
        if n_analyze < 1; return nothing; end
        analysis_start = (total_cycles - n_analyze) * cycle_length
        analysis_end   = total_cycles * cycle_length
        t_eval = collect(analysis_start:1.0:analysis_end)
    end

    # Validate initial condition vector length
    expected_len = 2 + length(strategy_indices) * 8
    if length(u0) != expected_len
        error("u0 length mismatch: expected $(expected_len), got $(length(u0))")
    end

    tspan = (0.0, t_target)
    prob  = ODEProblem(competition_model!, u0, tspan, (params, strategy_indices))

    try
        cb = PositiveDomain(u0; abstol=1e-10, scalefactor=0.5, save=false)
        sol = solve(prob, Tsit5();
                    saveat=t_eval, save_everystep=false, dense=false,
                    abstol=1e-10, reltol=1e-8,
                    maxiters=1_000_000_000_000,
                    callback=cb)
        return sol
    catch e
        println("  GRID SIMULATION FAILED: ", e)
        return nothing
    end
end


# ==============================================================================
# Initial Condition Setup
# ==============================================================================

"""
    setup_initial_conditions(params, strategy_indices, file, use_file, Z)

Standalone version: build u0 with file-based or default initial conditions.
Handles file reading, volume fraction mixing, and fallback to defaults internally.
Used when running the model definitions file independently (without run script).
`Z` is the volume fraction: strategy 1 gets Z, strategy 2 gets 1-Z.
"""

function setup_initial_conditions(params::Parameters, strategy_indices,
                                  initial_conditions_file, use_file_initial_conditions,
                                  volume_fraction)
    num_strategies = length(strategy_indices)
    u0 = [2.0, 2.0]  # Default [S1, S2]

    if use_file_initial_conditions && initial_conditions_file !== nothing && volume_fraction !== nothing
        equilibrium_data = read_initial_conditions_from_file(
            initial_conditions_file, params.T1, params.T2, strategy_indices
        )

        if equilibrium_data !== nothing
            internal_states, external_states = equilibrium_data

            # Mix external substrates by volume fraction (Eq. 26)
            ext1 = external_states[strategy_indices[1]]
            ext2 = external_states[strategy_indices[2]]
            u0[1] = volume_fraction * ext1[1] + (1 - volume_fraction) * ext2[1]
            u0[2] = volume_fraction * ext1[2] + (1 - volume_fraction) * ext2[2]

            for i in 1:num_strategies
                sidx = strategy_indices[i]
                append!(u0, copy(internal_states[sidx]))
                biomass = (i == 1) ?
                    volume_fraction * external_states[strategy_indices[1]][3] :
                    (1 - volume_fraction) * external_states[strategy_indices[2]][3]
                push!(u0, biomass)
            end
            return u0
        else
            println("Warning: Could not read from file, using defaults")
        end
    end

    # Fallback: default initial conditions
    u0 = [2.0, 2.0]
    for i in 1:num_strategies
        append!(u0, [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.88])
        biomass = (volume_fraction !== nothing) ?
            (i == 1 ? volume_fraction : 1 - volume_fraction) : 1.0
        push!(u0, biomass)
    end
    return u0
end


"""
    calculate_initial_conditions_from_file(params, strategy_indices, file, Z)

Build u0 purely from file data with volume fraction mixing.
Returns nothing if file read fails.
Returns nothing if file read fails, letting the caller handle fallback logic.
Called by run_competition_chemostat_grid.jl for pre-computing initial conditions.

"""
function calculate_initial_conditions_from_file(params::Parameters, strategy_indices,
                                               initial_conditions_file, volume_fraction)
    equilibrium_data = read_initial_conditions_from_file(
        initial_conditions_file, params.T1, params.T2, strategy_indices
    )
    if equilibrium_data === nothing
        return nothing
    end

    internal_states, external_states = equilibrium_data

    ext1 = external_states[strategy_indices[1]]
    ext2 = external_states[strategy_indices[2]]
    u0 = [volume_fraction * ext1[1] + (1 - volume_fraction) * ext2[1],
          volume_fraction * ext1[2] + (1 - volume_fraction) * ext2[2]]

    for i in 1:length(strategy_indices)
        sidx = strategy_indices[i]
        append!(u0, copy(internal_states[sidx]))
        biomass = (i == 1) ?
            volume_fraction * external_states[strategy_indices[1]][3] :
            (1 - volume_fraction) * external_states[strategy_indices[2]][3]
        push!(u0, biomass)
    end
    return u0
end


"""Generate default initial conditions with volume fraction split."""
function generate_default_initial_conditions(strategy_indices, volume_fraction)
    u0 = [2.0, 2.0]
    for i in 1:length(strategy_indices)
        append!(u0, [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.88])
        push!(u0, i == 1 ? volume_fraction : 1 - volume_fraction)
    end
    return u0
end


# ==============================================================================
# Volume Fraction Grid
# ==============================================================================

"""
    generate_volume_fraction_grid(grid_size; volume_min=1e-5) -> Vector{Float64}

Generate a symmetric logarithmic grid of volume fractions around 0.5.
Used to test competition outcomes across a range of initial biomass ratios.
"""
function generate_volume_fraction_grid(grid_size::Int, volume_min::Float64=1e-5)
    if grid_size < 1
        error("Grid size must be at least 1")
    end
    if volume_min <= 0 || volume_min >= 0.5
        error("Invalid volume_min: need 0 < volume_min < 0.5, got $(volume_min)")
    end

    grid_size == 1 && return [0.5]
    grid_size == 2 && return [volume_min, 1.0 - volume_min]

    half = grid_size ÷ 2

    if grid_size % 2 == 1
        left  = 10 .^ range(log10(volume_min), log10(0.5), length=half+1)[1:end-1]
        right = 1 .- reverse(left)
        return sort(vcat(left, [0.5], right))
    else
        left  = 10 .^ range(log10(volume_min), log10(0.5), length=half+1)[1:end-1]
        right = 1 .- reverse(left)
        return sort(vcat(left, right))
    end
end


# ==============================================================================
# Read Initial Conditions from Independent Culture Results
# ==============================================================================

"""
    read_initial_conditions_from_file(filepath, T1, T2, strategy_indices)

Read initial conditions from independent culture CSV.
Always uses the fixed reference condition T1=100.0, T2=0.0
(resource-rich environment) regardless of the competition's actual T1/T2.

Returns (internal_states, external_states) or nothing.
"""
function read_initial_conditions_from_file(filepath::String, T1::Float64, T2::Float64,
                                          strategy_indices::Vector{Int})
    fixed_T1 = 100.0
    fixed_T2 = 0.0

    cache_key = "$(filepath)_$(fixed_T1)_$(fixed_T2)"
    if haskey(INITIAL_CONDITIONS_CACHE, cache_key)
        return INITIAL_CONDITIONS_CACHE[cache_key]
    end

    println("Reading initial conditions (T1=$(fixed_T1), T2=$(fixed_T2)) from: $(filepath)")

    try
        df = CSV.read(filepath, DataFrame)
        matching = findall(row -> abs(row.T1 - fixed_T1) < 1e-6 && abs(row.T2 - fixed_T2) < 1e-6,
                           eachrow(df))

        if isempty(matching)
            println("Warning: No matching T1=$(fixed_T1), T2=$(fixed_T2) found in file")
            return nothing
        end

        row = df[matching[1], :]

        internal_result = Dict{Int, Vector{Float64}}()
        external_result = Dict{Int, Vector{Float64}}()

        for strategy_idx in strategy_indices
            internal_cols = [
                "strategy_$(strategy_idx)_X1_cycle_end",
                "strategy_$(strategy_idx)_M1_cycle_end",
                "strategy_$(strategy_idx)_E1_cycle_end",
                "strategy_$(strategy_idx)_X2_cycle_end",
                "strategy_$(strategy_idx)_M2_cycle_end",
                "strategy_$(strategy_idx)_E2_cycle_end",
                "strategy_$(strategy_idx)_c_minus_cycle_end"
            ]
            external_cols = [
                "strategy_$(strategy_idx)_S1_cycle_end",
                "strategy_$(strategy_idx)_S2_cycle_end",
                "strategy_$(strategy_idx)_c_cycle_end"
            ]

            missing_cols = [c for c in [internal_cols; external_cols] if !(c in names(df))]
            if !isempty(missing_cols)
                println("Warning: Missing columns for strategy $(strategy_idx): $(missing_cols)")
                continue
            end

            internal_result[strategy_idx] = [row[c] for c in internal_cols]
            external_result[strategy_idx] = [row[c] for c in external_cols]

            println("Strategy $(strategy_idx) internal: [$(join(string.(internal_result[strategy_idx]), ", "))]")
            println("Strategy $(strategy_idx) external: [$(join(string.(external_result[strategy_idx]), ", "))]")
        end

        if !isempty(internal_result)
            INITIAL_CONDITIONS_CACHE[cache_key] = (internal_result, external_result)
        end

        return isempty(internal_result) ? nothing : (internal_result, external_result)

    catch e
        println("Error reading initial conditions file: $(e)")
        return nothing
    end
end


# ==============================================================================
# Coexistence Classification (Two-Strategy Only)
# ==============================================================================

"""
    classify_coexistence_winner(surviving_strategies, avg_biomass, strategy_indices) -> Int

For two surviving strategies, return a winner code indicating biomass ranking:
  - 12 → strategy 1 > strategy 2
  - 21 → strategy 2 > strategy 1
  - 901 → strategy 0 & strategy 1 coexist, strategy 0 ≥ strategy 1
          (special code because "01" parses as 1, conflicting with exclusion)
  - 10 → strategy 0 & strategy 1 coexist, strategy 1 > strategy 0
"""
function classify_coexistence_winner(surviving_strategies, avg_biomass_last_period, strategy_indices)
    survivors = sort(surviving_strategies)

    if length(survivors) == 1
        return survivors[1]
    elseif length(survivors) == 2
        s1, s2 = survivors[1], survivors[2]
        idx1 = findfirst(x -> x == s1, strategy_indices)
        idx2 = findfirst(x -> x == s2, strategy_indices)
        b1 = avg_biomass_last_period[idx1]
        b2 = avg_biomass_last_period[idx2]

        if b1 > b2
            # Special case: "01" would parse as 1
            if s1 == 0 && s2 == 1
                return 901
            end
            return parse(Int, "$(s1)$(s2)")
        else
            return parse(Int, "$(s2)$(s1)")
        end
    else
        return -40  # Unexpected: more than 2 survivors
    end
end


# ==============================================================================
# Competition Outcome Analysis (Ne*s Framework)
# ==============================================================================

"""
    analyze_grid_simulation_outcome(sol, params, strategy_indices, Ne; ...) -> Dict

Full analysis pipeline for a single competition simulation:
  1. Extinction check (average biomass over last complete cycle)
  2. Log-linear regression for surviving strategies
  3. Ne*s neutral evolution analysis
  4. Final outcome classification
"""
function analyze_grid_simulation_outcome(sol, params::Parameters, strategy_indices, Ne::Float64;
                                        extinction_threshold::Float64=1e-100,
                                        n_cycles_for_regression::Int=50,
                                        verbose=false)
    if sol === nothing || length(sol.t) < 2
        return nothing
    end

    num_strategies = length(strategy_indices)

    # ------------------------------------------------------------------
    # 1. Average biomass over the last complete cycle
    # ------------------------------------------------------------------
    cycle_length = params.T1 + params.T2
    local avg_biomass_last_period

    if cycle_length > 0 && length(sol.t) > 1
        total_time      = sol.t[end] - sol.t[1]
        complete_cycles  = floor(total_time / cycle_length)
        last_start = sol.t[1] + (complete_cycles - 1) * cycle_length
        last_end   = sol.t[1] + complete_cycles * cycle_length
        mask       = (sol.t .>= last_start) .& (sol.t .<= last_end)
        avg_biomass_last_period = [mean(sol[2 + (i-1)*8 + 8, mask]) for i in 1:num_strategies]
    else
        avg_biomass_last_period = [sol[2 + (i-1)*8 + 8, end] for i in 1:num_strategies]
    end

    surviving_indices = findall(avg_biomass_last_period .> extinction_threshold)
    extinct_indices   = findall(avg_biomass_last_period .<= extinction_threshold)

    if verbose
        println("\n=== EXTINCTION ANALYSIS ===")
        for (i, sidx) in enumerate(strategy_indices)
            status = i in surviving_indices ? "SURVIVING" : "EXTINCT"
            println("  Strategy $(sidx): $(avg_biomass_last_period[i]) ($(status))")
        end
    end

    # ------------------------------------------------------------------
    # 2. Initialise per-strategy result tracking
    # ------------------------------------------------------------------
    individual_nes_results = Dict{Int, Dict{String, Any}}()
    for (i, sidx) in enumerate(strategy_indices)
        individual_nes_results[sidx] = Dict(
            "final_biomass"      => avg_biomass_last_period[i],
            "survives"           => i in surviving_indices,
            "survival_reason"    => i in extinct_indices ? "extinction" : "to_be_determined",
            "needs_nes_analysis" => false,
            "net_growth_rate"    => NaN,
            "r_squared"          => NaN,
            "reference_r"        => NaN,
            "s_gen"              => NaN,
            "Ne_s"               => NaN,
            "regression_valid"   => false,
            "result_type"        => "unknown"
        )
    end

    # Early exit: all extinct
    if isempty(surviving_indices)
        if verbose; println("\n=== RESULT: All strategies extinct ==="); end
        return Dict(
            "outcome"                => "Extinction",
            "winner"                 => -1,
            "survivors"              => Int[],
            "s_value"                => NaN,
            "Ne_s_product"           => NaN,
            "individual_nes_results" => individual_nes_results,
            "avg_biomass"            => avg_biomass_last_period,
            "refined_outcome"        => "Extinction"
        )
    end

    # ------------------------------------------------------------------
    # 3. Log-linear regression for surviving strategies
    # ------------------------------------------------------------------
    if verbose; println("\n=== LINEAR REGRESSION ANALYSIS ==="); end

    net_growth_rates, r_squared_values, regression_info = calculate_net_growth_rates_by_regression(
        sol, surviving_indices, strategy_indices, params;
        n_cycles=n_cycles_for_regression, verbose=verbose
    )

    # ------------------------------------------------------------------
    # 4. Classify regression slopes
    # ------------------------------------------------------------------
    strategies_zero_slope    = Int[]
    strategies_nonzero_slope = Int[]
    strategies_invalid_slope = Int[]

    for i in surviving_indices
        sidx = strategy_indices[i]

        if haskey(net_growth_rates, sidx)
            r_val  = net_growth_rates[sidx]
            r_sq   = get(r_squared_values, sidx, NaN)
            rtype  = get(regression_info[sidx], "result_type", "unknown")
            glm_ok = get(regression_info[sidx], "glm_success", false)

            individual_nes_results[sidx]["net_growth_rate"]  = r_val
            individual_nes_results[sidx]["r_squared"]        = r_sq
            individual_nes_results[sidx]["result_type"]      = rtype
            individual_nes_results[sidx]["regression_valid"] = glm_ok

            if verbose
                println("Strategy $(sidx): slope=$(r_val), R²=$(r_sq), type=$(rtype)")
            end

            if isfinite(r_val)
                if r_val == 0.0
                    push!(strategies_zero_slope, sidx)
                    individual_nes_results[sidx]["survival_reason"] = "zero_slope"
                    individual_nes_results[sidx]["s_gen"] = 0.0
                    individual_nes_results[sidx]["Ne_s"]  = 0.0
                else
                    push!(strategies_nonzero_slope, sidx)
                    individual_nes_results[sidx]["needs_nes_analysis"] = true
                end
            else
                push!(strategies_invalid_slope, sidx)
                individual_nes_results[sidx]["survival_reason"] = "invalid_slope_$(rtype)"
                individual_nes_results[sidx]["survives"] = false
            end
        else
            push!(strategies_invalid_slope, sidx)
            individual_nes_results[sidx]["survival_reason"] = "no_regression_results"
            individual_nes_results[sidx]["survives"] = false
        end
    end

    # ------------------------------------------------------------------
    # 5. Ne*s analysis for strategies with non-zero finite slopes
    # ------------------------------------------------------------------
    final_survivors = copy(strategies_zero_slope)

    if !isempty(strategies_nonzero_slope)
        if verbose
            println("\n=== Ne*s ANALYSIS ===")
            println("Zero slope (auto-survive): $(strategies_zero_slope)")
            println("Invalid slope (eliminated): $(strategies_invalid_slope)")
            println("Needing Ne*s: $(strategies_nonzero_slope)")
        end

        # Reference = maximum positive growth rate
        reference_r = 0.0
        for sidx in strategies_nonzero_slope
            r_val = net_growth_rates[sidx]
            if r_val > 0
                reference_r = max(reference_r, r_val)
            end
        end

        if verbose; println("Reference growth rate: $(reference_r)"); end

        for sidx in strategies_nonzero_slope
            r_val = net_growth_rates[sidx]
            individual_nes_results[sidx]["reference_r"] = reference_r

            if r_val > 0 && r_val == reference_r
                # Fastest grower: automatic survival
                individual_nes_results[sidx]["s_gen"]           = 0.0
                individual_nes_results[sidx]["Ne_s"]            = 0.0
                individual_nes_results[sidx]["survival_reason"] = "reference_strategy"
                individual_nes_results[sidx]["survives"]        = true
                push!(final_survivors, sidx)
            else
                # Standard Ne*s test
                s_gen = (r_val - reference_r) * (log(2) / params.D)
                Ne_s  = abs(Ne * s_gen)

                individual_nes_results[sidx]["s_gen"] = s_gen
                individual_nes_results[sidx]["Ne_s"]  = Ne_s

                survives = (Ne_s < 10.0)
                individual_nes_results[sidx]["survives"]        = survives
                individual_nes_results[sidx]["survival_reason"] = survives ? "drift_dominates" : "selection_eliminates"

                if survives
                    push!(final_survivors, sidx)
                end

                if verbose
                    println("Strategy $(sidx): r=$(round(r_val, digits=6)), " *
                            "s=$(round(s_gen, digits=6)), Ne*s=$(round(Ne_s, digits=2)) -> " *
                            "$(survives ? "SURVIVES" : "ELIMINATED")")
                end
            end
        end

        # Store regression details
        for sidx in strategies_nonzero_slope
            if haskey(regression_info, sidx)
                individual_nes_results[sidx]["regression_details"] = regression_info[sidx]
            end
        end
    end

    final_survivors = sort(unique(final_survivors))

    # ------------------------------------------------------------------
    # 6. Final classification
    # ------------------------------------------------------------------
    outcome, winner = if isempty(final_survivors)
        ("Extinction", -1)
    elseif length(final_survivors) == 1
        ("Exclusion", final_survivors[1])
    else
        survivor_str = join(final_survivors, "&")
        ("Coexistence ($(survivor_str))",
         classify_coexistence_winner(final_survivors, avg_biomass_last_period, strategy_indices))
    end

    # Summary Ne*s statistics
    s_value      = NaN
    Ne_s_product = NaN
    for sidx in strategy_indices
        if individual_nes_results[sidx]["needs_nes_analysis"] && !isnan(individual_nes_results[sidx]["Ne_s"])
            current = individual_nes_results[sidx]["Ne_s"]
            if isnan(Ne_s_product) || current > Ne_s_product
                Ne_s_product = current
                s_value      = individual_nes_results[sidx]["s_gen"]
            end
        end
    end

    # Final biomass per strategy
    final_biomass_dict = Dict{String, Float64}()
    for (i, sidx) in enumerate(strategy_indices)
        final_biomass_dict["strategy$(sidx)_final_biomass"] = avg_biomass_last_period[i]
    end

    if verbose
        println("\n" * "="^60)
        println("=== FINAL RESULT ===")
        println("  Outcome: $(outcome)")
        println("  Survivors: $(final_survivors)")
        if !isnan(s_value)
            println("  Max |s|: $(round(abs(s_value), digits=6))")
            println("  Max Ne*s: $(round(Ne_s_product, digits=2))")
        end
        println("="^60)
    end

    return merge(Dict(
        "outcome"                  => outcome,
        "winner"                   => winner,
        "survivors"                => final_survivors,
        "s_value"                  => s_value,
        "Ne_s_product"             => Ne_s_product,
        "avg_biomass"              => avg_biomass_last_period,
        "individual_nes_results"   => individual_nes_results,
        "refined_outcome"          => outcome,
        "strategies_zero_slope"    => strategies_zero_slope,
        "strategies_nonzero_slope" => strategies_nonzero_slope,
        "strategies_invalid_slope" => strategies_invalid_slope
    ), final_biomass_dict)
end