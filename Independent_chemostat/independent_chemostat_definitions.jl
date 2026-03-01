using DifferentialEquations
using LinearAlgebra
using Statistics
using CSV
using DataFrames
using Printf
using DiffEqCallbacks  # provide PositiveDomain callback

"""
Parameters structure for the CCR independent cultivation model.
This model simulates independent (monoculture) growth of different CCR strategies
to compare their steady-state biomass under identical environmental conditions.
"""
mutable struct IndependentParameters
    # Substrate uptake
    k1::Float64          # S1 maximum specific uptake rate
    k2::Float64          # S2 maximum specific uptake rate
    KS1::Float64         # S1 uptake half-saturation constant
    KS2::Float64         # S2 uptake half-saturation constant

    # Metabolism
    alpha1::Float64      # X1 metabolic rate constant
    alpha2::Float64      # X2 metabolic rate constant
 
    # Gene expression
    kM::Float64          # mRNA synthesis rate constant
    kE::Float64          # Enzyme synthesis rate constant (translational)
    delta1::Float64      # mRNA degradation rate constant
    delta2::Float64      # Enzyme degradation rate constant
    KX1::Float64         # Induction half-saturation constant for M1
    KX2::Float64         # Induction half-saturation constant for M2
    b::Float64           # Basal expression fraction 
 
    # CCR parameters
    K_CCR1::Float64      # Repression half-saturation constant (Transcriptional, Strategy T)
    K_CCR2::Float64      # Repression half-saturation constant (Post-transcriptional, Strategy P)
 
    # Yield coefficients
    Y1::Float64          # Yield coefficient for S1 pathway
    Y2::Float64          # Yield coefficient for S2 pathway

    # Operational parameters
    D::Float64           # Dilution rate
    Sr1::Float64         # S1 reservoir concentration
    Sr2::Float64         # S2 reservoir concentration (S2_in)
 
    # Wave parameters for S1 feeding
    epsilon::Float64     # S1 feeding amplitude
    k_wave::Float64      # Steepness parameter for S1 pulse function
    T1::Float64          # S1 ON duration
    T2::Float64          # S1 OFF duration
end
 
"""
Constructor with default parameter values for independent cultivation.
Default values correspond to Table 2 in the Model section.
"""
function IndependentParameters()
    # Get D from environment variable or use default
    D_value = haskey(ENV, "CCR_D_VALUE") ? parse(Float64, ENV["CCR_D_VALUE"]) : 0.04
    
    return IndependentParameters(
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

# ================================================================
# SUBSTRATE INPUT FUNCTIONS
# ================================================================

"""
S1 inflow concentration as a function of time (Eq. 4-5 in Model).
Implements a smooth pulse function with ON duration T1 and OFF duration T2.
"""
function S1_in(t, params::IndependentParameters)
    if params.T1 <= 0
        return 0.0
    end
    
    T0 = params.T1 + params.T2
    if T0 <= 0
        return 0.0
    end
    
    if params.T1 >= T0
        return params.epsilon
    end
    
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
 
"""
S2 inflow concentration (constant, Eq. 6 in Model).
"""
function S2_in(t, params::IndependentParameters)
    return params.Sr2
end

# ================================================================
# ODE SYSTEM
# ================================================================

"""
ODE system for independent cultivation of multiple CCR strategies.
Each strategy is simulated in parallel within the same ODE system but with
separate state variables (no interaction between strategies).

Strategy indices:
- 0: Strategy N (No CCR)
- 1: Strategy T (Transcriptional CCR)
- 2: Strategy P (Post-transcriptional CCR)

State variables per strategy (10 variables, offset = (i-1)*10):
  [S1, S2, X1, M1, E1, X2, M2, E2, c_minus, c]
"""
function independent_cultivation_model!(du, u, p, t)
    params, strategy_indices = p
    
    for (i, strategy_idx) in enumerate(strategy_indices)
        offset = (i-1) * 10

        # 1) Read state variables
        S1 = u[offset+1]; S2 = u[offset+2]; X1 = u[offset+3]
        M1 = u[offset+4]; E1 = u[offset+5]
        X2 = u[offset+6]; M2 = u[offset+7]; E2 = u[offset+8]
        c_minus = u[offset+9]; c = u[offset+10]

        # 2) Growth rate (Eq. 20)
        rg = (params.k1 * E1 * S1 / (params.KS1 + S1) + 
              params.k2 * E2 * S2 / (params.KS2 + S2)) - 
             ((1 - params.Y1) * params.alpha1 * X1 + 
              (1 - params.Y2) * params.alpha2 * X2)

        # 3) CCR regulation functions (Eq. 12-13, 16-17)
        #    kappa_M: transcriptional CCR modifier for M2 synthesis
        #    kappa_E: post-transcriptional CCR modifier for E2 synthesis
        if strategy_idx == 0       # Strategy N: No CCR
            kappa_M = 1.0
            kappa_E = 1.0
        elseif strategy_idx == 1   # Strategy T: Transcriptional CCR (Eq. 13)
            kappa_M = params.K_CCR1 / (params.K_CCR1 + X1)
            kappa_E = 1.0
        elseif strategy_idx == 2   # Strategy P: Post-transcriptional CCR (Eq. 17)
            kappa_M = 1.0
            kappa_E = params.K_CCR2 / (params.K_CCR2 + X1)
        end

        # 4) ODEs
        # dS1/dt (Eq. 2) — independent cultivation: only one strategy per chemostat
        du[offset+1] = (S1_in(t, params) - S1) * params.D - 
                        params.k1 * E1 * S1 / (params.KS1 + S1) * c
        # dS2/dt (Eq. 3)
        du[offset+2] = (S2_in(t, params) - S2) * params.D - 
                        params.k2 * E2 * S2 / (params.KS2 + S2) * c
        # dX1/dt (Eq. 8)
        du[offset+3] = params.k1 * E1 * S1 / (params.KS1 + S1) - 
                        params.alpha1 * X1 - rg * X1
        # dM1/dt (Eq. 10)
        du[offset+4] = params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * c_minus - 
                        params.delta1 * M1 - rg * M1
        # dE1/dt (Eq. 14)
        du[offset+5] = params.kE * M1 * c_minus - 
                        params.delta2 * E1 - rg * E1
        # dX2/dt (Eq. 9)
        du[offset+6] = params.k2 * E2 * S2 / (params.KS2 + S2) - 
                        params.alpha2 * X2 - rg * X2
        # dM2/dt (Eq. 11) — includes kappa_M for transcriptional CCR
        du[offset+7] = params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * c_minus * kappa_M - 
                        params.delta1 * M2 - rg * M2
        # dE2/dt (Eq. 15) — includes kappa_E for post-transcriptional CCR
        du[offset+8] = params.kE * M2 * c_minus * kappa_E - 
                        params.delta2 * E2 - rg * E2
        # dc_minus/dt (Eq. 18)
        du[offset+9] = (params.Y1 * params.alpha1 * X1 + params.Y2 * params.alpha2 * X2) +
                        (params.delta1 * M1 + params.delta2 * E1 +
                         params.delta1 * M2 + params.delta2 * E2) -
                        (params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * c_minus +
                         params.kE * M1 * c_minus +
                         params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * c_minus * kappa_M +
                         params.kE * M2 * c_minus * kappa_E) -
                        rg * c_minus
        # dc/dt (Eq. 21)
        du[offset+10] = (rg - params.D) * c

        # 5) Non-negative constraint check (backup)
        for j in (offset+1):(offset+10)
            if u[j] <= 0.0 && du[j] < 0.0
                du[j] = 0.0
            end
        end
    end
end

# ================================================================
# SIMULATION FUNCTION
# ================================================================

"""
Run independent cultivation simulation for specified CCR strategies.

Arguments:
- params: IndependentParameters struct
- strategy_indices: Vector of strategy indices to simulate (e.g., [0, 1] for N vs T, [1, 2] for T vs P)
- target_tau: Target biological time in generations (default: 1e7, overridden by run file)
- analysis_cycles: Number of final cycles to save for analysis
- initial_biomasses: Initial biomass for each strategy

Returns:
- Dict with "solution" key containing the ODE solution, or nothing if failed
"""
function run_independent_simulation_method_a_cycles(params::IndependentParameters, 
                                                     strategy_indices::Vector{Int};
                                                     target_tau::Float64=1e7,
                                                     analysis_cycles::Int=1,
                                                     initial_biomasses::Vector{Float64}=Float64[])
    # Convert tau to physical time (Eq. 25: t = tau * ln(2) / D)
    t_target = (target_tau * log(2)) / params.D
    cycle_length = params.T1 + params.T2

    # Determine analysis time window (save only the last analysis_cycles cycles)
    local t_eval
    if cycle_length <= 0
        analysis_start_time = max(0.0, t_target - analysis_cycles)
        t_eval = collect(analysis_start_time:0.1:t_target)
    else
        total_complete_cycles = floor(Int, t_target / cycle_length)
        n_cycles_to_analyze = min(analysis_cycles, total_complete_cycles)
        if n_cycles_to_analyze < 1
            return nothing
        end
        analysis_start_time = (total_complete_cycles - n_cycles_to_analyze) * cycle_length
        analysis_end_time = total_complete_cycles * cycle_length
        t_eval = collect(analysis_start_time:0.1:analysis_end_time)
    end

    # Setup initial conditions 
    num_strategies = length(strategy_indices)
    if isempty(initial_biomasses)
        initial_biomasses = fill(0.01, num_strategies)
    end

    # Initial state: S1=2, S2=2, X1=M1=E1=X2=M2=E2=0.02, c_minus=0.88
    initial_vars_template = [2.0, 2.0, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.88]

    u0 = Float64[]
    for i in 1:num_strategies
        strategy_vars = copy(initial_vars_template)
        push!(strategy_vars, initial_biomasses[i])
        append!(u0, strategy_vars)
    end

    # Solve ODE
    tspan = (0.0, t_target)
    prob = ODEProblem(independent_cultivation_model!, u0, tspan, (params, strategy_indices))

    try
        cb = PositiveDomain(u0; abstol=1e-12, scalefactor=0.5, save=false)
        sol = solve(prob, Tsit5();
                    saveat = t_eval,
                    tstops = t_eval,
                    save_everystep = false,
                    dense = false,
                    abstol=1e-12, reltol=1e-10,
                    callback = cb,
                    maxiters = 1e12)

        if sol.retcode != :Success
            return nothing
        end

        return Dict("solution" => sol)
    catch e
        println("  SIMULATION FAILED: ", e)
        return nothing
    end
end

# ================================================================
# ANALYSIS FUNCTIONS
# ================================================================

"""
Analyze simulation outcome by calculating average biomass in the last cycle .

For each strategy, computes the average biomass over the final complete cycle (T1+T2).
Strategies with average biomass below the extinction threshold are classified as extinct.

Returns:
- Dict with outcome, survivors, avg_biomass, and individual results
"""
function analyze_independent_simulation_outcome(sol, params::IndependentParameters, strategy_indices;
                                              extinction_threshold::Float64=1e-100, 
                                              verbose=false)
    if sol === nothing || length(sol.t) < 2
        return nothing
    end
    
    num_strategies = length(strategy_indices)
    cycle_length = params.T1 + params.T2
    local avg_biomass_last_period
    
    # Calculate average biomass in the last complete cycle
    if cycle_length > 0 && length(sol.t) > 1
        total_time = sol.t[end] - sol.t[1]
        complete_cycles = floor(total_time / cycle_length)
        
        if complete_cycles >= 1
            last_cycle_start = sol.t[1] + (complete_cycles - 1) * cycle_length
            last_cycle_end = sol.t[1] + complete_cycles * cycle_length
            
            # Get data points in the last full cycle
            last_cycle_mask = (sol.t .>= last_cycle_start) .& (sol.t .<= last_cycle_end)
            
            if sum(last_cycle_mask) > 1
                avg_biomass_last_period = [mean(sol[(i-1)*10 + 10, last_cycle_mask]) for i in 1:num_strategies]
            else
                # Fallback if data points are scarce
                avg_biomass_last_period = [sol[(i-1)*10 + 10, end] for i in 1:num_strategies]
            end
        else
             avg_biomass_last_period = [sol[(i-1)*10 + 10, end] for i in 1:num_strategies]
        end
    else
        avg_biomass_last_period = [sol[(i-1)*10 + 10, end] for i in 1:num_strategies]
    end
    
    surviving_indices = findall(avg_biomass_last_period .> extinction_threshold)
    
    # Build individual results
    individual_results = Dict{Int, Dict{String, Any}}()
    
    for (i, strategy_idx) in enumerate(strategy_indices)
        is_survivor = i in surviving_indices
        individual_results[strategy_idx] = Dict(
            "final_biomass" => avg_biomass_last_period[i],
            "survives" => is_survivor,
            "survival_reason" => is_survivor ? "above_threshold" : "extinction"
        )
    end
    
    outcome = isempty(surviving_indices) ? "All_Extinct" : "Survivors: $([strategy_indices[i] for i in surviving_indices])"
    
    return Dict(
        "outcome" => outcome, 
        "survivors" => [strategy_indices[i] for i in surviving_indices],
        "avg_biomass" => avg_biomass_last_period,
        "individual_results" => individual_results
    )
end

"""
Compare independent strategies based on their average biomass.

For independent cultivation, we compare which strategy achieves higher biomass
under identical environmental conditions. This is NOT competitive exclusion,
but rather a comparison of growth capacity.

Possible outcomes:
- "Strategy_X_higher_biomass": One strategy has clearly higher biomass
- "Nearly_equal_biomass": Both strategies have nearly identical biomass
- "Both_extinct": Both strategies went extinct
"""
function compare_independent_models(analysis_results, sol, params::IndependentParameters, strategy_indices; 
                                   verbose=false)
    survivors = analysis_results["survivors"]
    individual_results = analysis_results["individual_results"]
    
    comparison_key = "$(strategy_indices[1])_vs_$(strategy_indices[2])"
    winner_key = "$(comparison_key)_winner"
    comparison_results = Dict()
    
    strategyA_idx = strategy_indices[1]
    strategyB_idx = strategy_indices[2]
    
    # Get biomass from analysis results
    biomass_A = get(get(individual_results, strategyA_idx, Dict()), "final_biomass", 0.0)
    biomass_B = get(get(individual_results, strategyB_idx, Dict()), "final_biomass", 0.0)

    if strategyA_idx in survivors && strategyB_idx in survivors
        biomass_diff = abs(biomass_A - biomass_B)
        
        if biomass_diff < 1e-10
            comparison_results[comparison_key] = "Nearly_equal_biomass"
            comparison_results[winner_key] = parse(Int, "$(strategyB_idx)$(strategyA_idx)") 
        elseif biomass_A > biomass_B
            comparison_results[comparison_key] = "Strategy_$(strategyA_idx)_higher_biomass"
            comparison_results[winner_key] = strategyA_idx
        else
            comparison_results[comparison_key] = "Strategy_$(strategyB_idx)_higher_biomass"
            comparison_results[winner_key] = strategyB_idx
        end
    elseif strategyA_idx in survivors
        comparison_results[comparison_key] = "Strategy_$(strategyA_idx)_higher_biomass"
        comparison_results[winner_key] = strategyA_idx
    elseif strategyB_idx in survivors
        comparison_results[comparison_key] = "Strategy_$(strategyB_idx)_higher_biomass"
        comparison_results[winner_key] = strategyB_idx
    else
        comparison_results[comparison_key] = "Both_extinct"
        comparison_results[winner_key] = -1
    end
    
    return comparison_results
end

"""
Extract state variables at the end of the last complete cycle.
Used to provide steady-state initial conditions for co-culture competition experiments.
"""
function extract_last_cycle_end_variables(sol, params::IndependentParameters, strategy_indices; verbose=false)
    cycle_length = params.T1 + params.T2
    last_cycle_variables = Dict{Int, Vector{Float64}}()
    
    if cycle_length > 0 && length(sol.t) > 0
        total_time = sol.t[end] - sol.t[1]
        complete_cycles = floor(total_time / cycle_length)
        last_cycle_end_time = sol.t[1] + complete_cycles * cycle_length
        closest_idx = argmin(abs.(sol.t .- last_cycle_end_time))
        cycle_time = last_cycle_end_time
    else
        closest_idx = length(sol.t)
        cycle_time = sol.t[end]
    end
    
    for (i, strategy_idx) in enumerate(strategy_indices)
        offset = (i-1) * 10
        # Extract all 10 variables: [S1, S2, X1, M1, E1, X2, M2, E2, c_minus, c]
        strategy_variables = [sol[offset+j, closest_idx] for j in 1:10]
        last_cycle_variables[strategy_idx] = strategy_variables
    end
    
    return last_cycle_variables, cycle_time
end