# ==============================================================================
# Competition Chemostat Time-Series Visualization
#
# Standalone module for time-series analysis and visualization
# chemostat competition simulations among Strategy N (no CCR, strategy_idx=0),
# Strategy T (transcriptional CCR, strategy_idx=1), and Strategy P
# (post-transcriptional CCR, strategy_idx=2).
#
# This file provides:
#   - Parameter struct and defaults (matching definitions file)
#   - S1 pulsed inflow function
#   - ODE system (competition_model!)
#   - Initial condition reading from independent culture CSV
#   - Last-cycle QuadGK integration analysis
#   - P1-P5 visualization functions
#
# Workflow:
#   1. Read steady-state ICs from independent culture CSV (T1=100, T2=0)
#   2. Mix with volume fraction Z=0.5 for competition
#   3. Run long-term simulation (target_tau generations)
#   4. Extract last complete cycle for dense output
#   5. QuadGK integration for synthesis rate analysis
#   6. Generate P1-P5 figures
# ==============================================================================

using DifferentialEquations
using DiffEqCallbacks
using QuadGK
using Printf
using Plots
using CSV
using DataFrames


# ==============================================================================
# Placeholder path prefix and IC cache
# ==============================================================================
const IC_FILE_PATH_PREFIX = joinpath(@__DIR__, "initial_data")
const INITIAL_CONDITIONS_CACHE = Dict{String, Tuple{Dict{Int, Vector{Float64}}, Dict{Int, Vector{Float64}}}}()


# ============================================================================
# Parameter Struct
# ============================================================================
mutable struct Parameters
    # Substrate uptake ("Hardware")
    k1::Float64          # S1 maximum specific uptake rate
    k2::Float64          # S2 maximum specific uptake rate
    KS1::Float64         # S1 uptake half-saturation constant
    KS2::Float64         # S2 uptake half-saturation constant
 
    # Metabolism ("Hardware")
    alpha1::Float64      # X1 specific consumption rate
    alpha2::Float64      # X2 specific consumption rate

    # Gene expression ("Software" / Global properties)
    kM::Float64          # Max specific mRNA synthesis rate
    kE::Float64          # Translational rate constant
    delta1::Float64      # mRNA degradation rate
    delta2::Float64      # Protein degradation rate
    KX1::Float64         # Induction half-saturation constant for pathway 1
    KX2::Float64         # Induction half-saturation constant for pathway 2
    b::Float64           # Basal expression fraction
    
    # CCR parameters
    K_CCR1::Float64      # CCR strength for transcriptional level
    K_CCR2::Float64      # CCR strength for post-transcriptional level
 
    # Yield and operational parameters
    Y1::Float64          # Yield coefficient for pathway 1
    Y2::Float64          # Yield coefficient for pathway 2
    D::Float64           # Dilution rate
    Sr1::Float64         # S1 reservoir concentration
    Sr2::Float64         # S2 reservoir concentration
 
    # Wave parameters for S1 feeding
    epsilon::Float64     # S1 feeding amplitude
    k_wave::Float64      # Steepness of switching function (renamed from 'k')
    T1::Float64          # S1 ON duration
    T2::Float64          # S1 OFF duration
end



function Parameters()
    
 
 
    
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
        0.000,  # delta2
        0.01, # KX1
        0.01, # KX2
        0.001,  # b

        # CCR parameters
        0.01,   # K_CCR1
        0.01,   # K_CCR2
        
        # Yield and operational parameters
        0.7,    # Y1
        0.7,   # Y2
        0.06,  # D 
        0.0,    # Sr1
        2.0,    # Sr2
        
        # Wave parameters
        2.0,    # epsilon
        10.0,   # k_wave
        150.0,  # T1
        50.0     # T2  
    )
end





"""
Add T1/T2 phase background coloring to a plot.
- cycle_start: cycle start time (from sol.t[1])
- T1 phase: [cycle_start, cycle_start + T1] → light red
- T2 phase: [cycle_start + T1, cycle_start + T1 + T2] → light blue
"""
function add_phase_background!(p, cycle_start, params;
                               t1_color=:red, t2_color=:lightblue,
                               t1_alpha=0.15, t2_alpha=0.15)
    y_lims = ylims(p)
    y_min, y_max = y_lims

    t1_end = cycle_start + params.T1
    t2_end = t1_end + params.T2

    # T1 phase background (light red)
    plot!(p, [cycle_start, t1_end, t1_end, cycle_start, cycle_start],
          [y_min, y_min, y_max, y_max, y_min],
          fill=true, fillalpha=t1_alpha, fillcolor=t1_color,
          linewidth=0, label="", seriestype=:shape)

    # T2 phase background (light blue)
    plot!(p, [t1_end, t2_end, t2_end, t1_end, t1_end],
          [y_min, y_min, y_max, y_max, y_min],
          fill=true, fillalpha=t2_alpha, fillcolor=t2_color,
          linewidth=0, label="", seriestype=:shape)

    ylims!(p, y_lims)
    return p
end


"""
Add colored background to difference plots based on which strategy is larger.
- When diff > 0 (strategy1 > strategy2): use strategy1's color
- When diff < 0 (strategy2 > strategy1): use strategy2's color

To avoid color mixing with T1/T2 background, call this function twice:
1. First with [:white, :white] and alpha=1.0 to cover the background
2. Then with strategy_colors and desired alpha to add the strategy color

This version merges consecutive segments of the same sign into single polygons
to avoid rendering issues when many small rectangles are drawn in small subplots.
"""
function add_difference_background!(sp, t_eval, diff_data, strategy_colors;
                                    alpha=0.2)
    y_lims = ylims(sp)
    n = length(t_eval)

    # Find continuous regions of positive and negative differences
    # and draw each region as a single filled polygon

    i = 1
    while i <= n
        # Skip points where diff is essentially zero
        if abs(diff_data[i]) < 1e-12
            i += 1
            continue
        end

        # Determine the sign of current region
        current_sign = sign(diff_data[i])

        # Find the end of this continuous region (same sign)
        region_start = i
        while i <= n && sign(diff_data[i]) == current_sign
            i += 1
        end
        region_end = min(i, n)  # Last index included in this region

        # Need at least 2 points to form a region
        if region_end <= region_start
            continue
        end

        # Build polygon for this region: follow the curve, then back along y=0
        t_points = t_eval[region_start:region_end]
        d_points = diff_data[region_start:region_end]

        # Create closed polygon: curve points + return along y=0
        poly_t = vcat(t_points, reverse(t_points), t_points[1])
        poly_y = vcat(d_points, zeros(length(t_points)), d_points[1])

        # Choose color based on sign
        if current_sign > 0
            fill_color = strategy_colors[1]
        else
            fill_color = strategy_colors[2]
        end

        # Draw the filled polygon
        plot!(sp, poly_t, poly_y,
              fill=true, fillalpha=alpha, fillcolor=fill_color,
              linewidth=0, label="", seriestype=:shape)
    end

    ylims!(sp, y_lims)
    return sp
end



# ============================================================================
# Substrate Inflow Functions
# ============================================================================
function S1_in(t, params::Parameters)
    # Handle special cases
    if params.T1 <= 0
        return 0.0
    end
    
    # Step 1: Calculate the total period from T1 (ON) and T2 (OFF) durations
    T0 = params.T1 + params.T2
    if T0 <= 0
        # Handle non-periodic case
        return 0.0
    end
    
    # If T1 >= T0 (i.e., T2 <= 0), return constant epsilon
    if params.T1 >= T0
        return params.epsilon
    end
    
    # Step 2: Calculate the time within the current cycle using the modulo operator
    t_in_cycle = t % T0
    # Step 3: Gated logic - only calculate the pulse if t is strictly within the (0, T1) interval
    if 0 < t_in_cycle < params.T1
        
        # --- Simplified Logic: Rise/Fall centers are fixed at the boundaries ---
        t_rise_center = 0.0
        t_fall_center = params.T1
        # --- Pulse Calculation using the difference of two logistic curves ---
        
        # The rise curve is centered at t=0
        rise_curve = 1.0 / (1.0 + exp(-params.k_wave * (t_in_cycle - t_rise_center)))
        # The fall curve is centered at t=T1
        fall_curve = 1.0 / (1.0 + exp(-params.k_wave * (t_in_cycle - t_fall_center)))
        
        # --- Zero-Offset Correction (for smooth start at zero) ---
        # Calculate the value of the raw pulse at t=0 to create the offset
        offset_rise = 1.0 / (1.0 + exp(-params.k_wave * (0.0 - t_rise_center))) # This is 0.5
        offset_fall = 1.0 / (1.0 + exp(-params.k_wave * (0.0 - t_fall_center)))
        
        # This is the pulse shape, corrected to start at 0. Its peak height is ~0.5.
        corrected_shape = (rise_curve - fall_curve) - (offset_rise - offset_fall)
        
        # --- Amplitude Scaling (to match epsilon) ---
        # Multiply by epsilon for amplitude and by 2.0 to correct the height.
        return params.epsilon * 2.0 * corrected_shape
    else
        # For t=0, t=T1, and the silent T2 phase, the value is exactly 0.
        return 0.0
    end
end


S2_in(t, params::Parameters) = params.Sr2


# ============================================================================
# ODE System: Pairwise Chemostat Competition
# ============================================================================


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

        # Specific growth rate
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





#c- constant
#=
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
        du[offset + 2] = params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * 0.1 -
                          params.delta1 * M1 - rg * M1
        # dE1/dt (Eq. 14)
        du[offset + 3] = params.kE * M1 * 0.1 -
                          params.delta2 * E1 - rg * E1
        # dX2/dt (Eq. 9)
        du[offset + 4] = params.k2 * E2 * S2 / (params.KS2 + S2) -
                          params.alpha2 * X2 - rg * X2
        # dM2/dt (Eq. 11) — includes kappa_M for transcriptional CCR
        du[offset + 5] = params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * 0.1 * kappa_M -
                          params.delta1 * M2 - rg * M2
        # dE2/dt (Eq. 15) — includes kappa_E for post-transcriptional CCR
        du[offset + 6] = params.kE * M2 * 0.1 * kappa_E -
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
=#

# ------------------------------------------------------------------------------
# Accelerating mRNA dynamics version
# K_tau = 1000 accelerates M1 and M2 dynamics only, preserving mass conservation
# in dc_minus/dt by scaling corresponding mRNA synthesis and degradation terms.
# To activate: uncomment this block and comment out the original version above
# ------------------------------------------------------------------------------
#=
function competition_model!(du, u, p, t)
    params, strategy_indices = p

    S1 = u[1]
    S2 = u[2]

    # Substrate dynamics: inflow/dilution (consumption added per strategy below)
    dS1dt = (S1_in(t, params) - S1) * params.D
    dS2dt = (S2_in(t, params) - S2) * params.D

    # Acceleration factor for mRNA dynamics
    K_tau = 100

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

        # ODEs with K_tau acceleration on mRNA dynamics
        # dX1/dt (Eq. 8)
        du[offset + 1] = params.k1 * E1 * S1 / (params.KS1 + S1) -
                          params.alpha1 * X1 - rg * X1
        # dM1/dt (Eq. 10, accelerated)
        du[offset + 2] = K_tau * (params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * c_minus -
                          params.delta1 * M1) - rg * M1
        # dE1/dt (Eq. 14)
        du[offset + 3] = params.kE * M1 * c_minus -
                          params.delta2 * E1 - rg * E1
        # dX2/dt (Eq. 9)
        du[offset + 4] = params.k2 * E2 * S2 / (params.KS2 + S2) -
                          params.alpha2 * X2 - rg * X2
        # dM2/dt (Eq. 11, accelerated, includes kappa_M)
        du[offset + 5] = K_tau * (params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * c_minus * kappa_M -
                          params.delta1 * M2) - rg * M2
        # dE2/dt (Eq. 15, includes kappa_E)
        du[offset + 6] = params.kE * M2 * c_minus * kappa_E -
                          params.delta2 * E2 - rg * E2
        # dc_minus/dt (Eq. 18, with K_tau-consistent mRNA terms)
        du[offset + 7] = (params.Y1 * params.alpha1 * X1 + params.Y2 * params.alpha2 * X2) +
                          (params.delta1 * (M1 + M2) * K_tau + params.delta2 * (E1 + E2)) -
                          (params.kE * M1 * c_minus + params.kE * M2 * c_minus * kappa_E) -
                          (params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * c_minus +
                           params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * c_minus * kappa_M) * K_tau -
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
=#


# ============================================================================
# Read initial conditions from independent culture CSV results
# ============================================================================
function read_initial_conditions_from_file(filepath::String, T1::Float64, T2::Float64, strategy_indices::Vector{Int})
    """
    Read initial conditions for single T1,T2 combination from file
    Returns Dict{Int, Vector{Float64}} mapping strategy_idx to [X1, M1, E1, X2, M2, E2, c_minus]
    PLUS Dict{Int, Vector{Float64}} for external states [S1, S2, c]
    
    MODIFIED: Always uses T1=100.0, T2=0.0 regardless of input parameters
    """
    
    # MODIFICATION: Fix T1,T2 to specific values instead of using input parameters
    fixed_T1 = 100.0
    fixed_T2 = 0.0
    
    # Simple cache key for this specific T1,T2
    cache_key = "$(filepath)_$(fixed_T1)_$(fixed_T2)"
    if haskey(INITIAL_CONDITIONS_CACHE, cache_key)
        return INITIAL_CONDITIONS_CACHE[cache_key]
    end
    
    println("Reading initial conditions for FIXED T1=$(fixed_T1), T2=$(fixed_T2) from: $(filepath)")
    
    try
        # Read CSV file
        df = CSV.read(filepath, DataFrame)
        
        # Find matching row using FIXED values
        matching_rows = findall(row -> abs(row.T1 - fixed_T1) < 1e-6 && abs(row.T2 - fixed_T2) < 1e-6, eachrow(df))
        
        if isempty(matching_rows)
            println("Warning: No matching T1=$(fixed_T1), T2=$(fixed_T2) found in file")
            return nothing
        end
        
        # Use first matching row
        row = df[matching_rows[1], :]
        
        # Extract initial conditions for requested strategies
        internal_result = Dict{Int, Vector{Float64}}()
        external_result = Dict{Int, Vector{Float64}}()
        
        for strategy_idx in strategy_indices
            # Internal variables (original 7 variables)
            # CSV column names must match the output format from independent culture simulations
            internal_cols = [
                "strategy_$(strategy_idx)_X1_cycle_end",
                "strategy_$(strategy_idx)_M1_cycle_end", 
                "strategy_$(strategy_idx)_E1_cycle_end",
                "strategy_$(strategy_idx)_X2_cycle_end",
                "strategy_$(strategy_idx)_M2_cycle_end",
                "strategy_$(strategy_idx)_E2_cycle_end",
                "strategy_$(strategy_idx)_c_minus_cycle_end"
            ]
            
            # External variables (new: S1, S2, biomass)
            external_cols = [
                "strategy_$(strategy_idx)_S1_cycle_end",
                "strategy_$(strategy_idx)_S2_cycle_end",
                "strategy_$(strategy_idx)_c_cycle_end"
            ]
            
            # Check if all columns exist
            missing_cols = [col for col in [internal_cols; external_cols] if !(col in names(df))]
            if !isempty(missing_cols)
                println("Warning: Missing columns for strategy $(strategy_idx): $(missing_cols)")
                continue
            end
            
            # Extract internal values (keep original format)
            internal_values = [row[col] for col in internal_cols]
            internal_result[strategy_idx] = internal_values
            
            # Extract external values (new)
            external_values = [row[col] for col in external_cols]
            external_result[strategy_idx] = external_values
        end
        
        # Cache both results as a tuple
        if !isempty(internal_result)
            INITIAL_CONDITIONS_CACHE[cache_key] = (internal_result, external_result)
        end
        
        return isempty(internal_result) ? nothing : (internal_result, external_result)
        
    catch e
        println("Error reading initial conditions file: $(e)")
        return nothing
    end
end

# ============================================================================
# Build u0 from file data with volume fraction mixing
# ============================================================================
function calculate_initial_conditions_from_file(params::Parameters, strategy_indices, 
                                                initial_conditions_file, volume_fraction)
    """
    Calculate initial conditions from file data with volume fraction mixing
    
    Returns:
    - Vector{Float64}: complete u0 vector, or nothing if failed
    """
    
    # Read equilibrium states from file
    equilibrium_data = read_initial_conditions_from_file(
        initial_conditions_file, params.T1, params.T2, strategy_indices
    )
    
    if equilibrium_data === nothing
        return nothing
    end
    
    internal_states, external_states = equilibrium_data
    
    # Check if we have data for all requested models
    if !all(haskey(internal_states, idx) && haskey(external_states, idx) for idx in strategy_indices)
        println("Warning: File is missing initial condition data for one or more strategies in $(strategy_indices).")
        return nothing
    end

    # Calculate mixed S1, S2 from file data
    external1 = external_states[strategy_indices[1]]
    external2 = external_states[strategy_indices[2]]
    S1_mixed = volume_fraction * external1[1] + (1 - volume_fraction) * external2[1]
    S2_mixed = volume_fraction * external1[2] + (1 - volume_fraction) * external2[2]
    
    u0 = [S1_mixed, S2_mixed]
    
    # Add strategy variables
    num_strategies = length(strategy_indices)
    for i in 1:num_strategies
        strategy_idx = strategy_indices[i]
        
        # Internal variables (7 variables)
        internal_vars = copy(internal_states[strategy_idx])
        append!(u0, internal_vars)
        
        # Biomass (mixed based on volume fraction)
        if i == 1
            biomass = volume_fraction * external_states[strategy_indices[1]][3]
        else
            biomass = (1 - volume_fraction) * external_states[strategy_indices[2]][3]
        end
        push!(u0, biomass)
    end
    
    return u0
end

# ============================================================================
# Last cycle integration analysis using QuadGK
# ============================================================================
function compute_last_cycle_integrals_quadgk(strategy_indices, params, u0, 
                                             sim_time)
    
    println("\n" * "="^80)
    println("LAST CYCLE INTEGRATION ANALYSIS (Two-Stage Optimization)")
    println("="^80)
    
    # 1. Compute last cycle start/end times

    T0 = params.T1 + params.T2
    if T0 <= 0
        println("Warning: T1+T2 is zero or negative. Using sim_time as cycle.")
        T0 = sim_time
    end

    num_complete_cycles = floor(sim_time / T0)
    cycle_start = (num_complete_cycles - 1) * T0
    cycle_end = num_complete_cycles * T0
    cycle_duration = cycle_end - cycle_start
    
    println("\nCycle Information:")
    println("  Simulation time: $sim_time")
    println("  Cycle period: $(T0)")
    println("  Last cycle: [$cycle_start, $cycle_end]")
    println("  Duration: $cycle_duration")
    
    # (IC setup moved to main function)
    
    # ========================================================================
    # STAGE 1: LONG-TERM SIMULATION (Competition dynamics until steady state)
    # ========================================================================
    println("\n⏳ Stage 1: Long-term simulation (0 -> $cycle_start)...")

    prob_longterm = ODEProblem(competition_model!, u0, (0.0, cycle_start), (params, strategy_indices))

    cb_longterm = PositiveDomain(u0; abstol=1e-10, scalefactor=0.5, save=false)

    # Create time points for saving biomass history (one point per time unit)
    t_longterm_save = 0.0:1.0:cycle_start

    sol_longterm = solve(prob_longterm, Tsit5(),
                callback=cb_longterm,
                abstol=1e-10, reltol=1e-8, maxiters=1e12,
                saveat=t_longterm_save)  # Save full trajectory for P1


    if sol_longterm.retcode != :Success
        error("Stage 1 long-term simulation failed: $(sol_longterm.retcode)")
    end

    u0_steady = sol_longterm[end]
    println("✓ Long-term simulation complete. System reached steady state.")

    # ========================================================================
    # STAGE 2: ANALYSIS ( dense outut last cycle)
    # ========================================================================
    println("\n⏳ Stage 2: High-res Simulation ($cycle_start -> $cycle_end)...")
    
    prob_analysis = remake(prob_longterm, u0=u0_steady, tspan=(cycle_start, cycle_end))
    
    cb_analysis = PositiveDomain(u0_steady; abstol=1e-10, scalefactor=0.5, save=false)

    sol = solve(prob_analysis, Tsit5(),
                callback=cb_analysis,
                abstol=1e-10, reltol=1e-8, maxiters=1e12)
                
    if sol.retcode != :Success
        error("Stage 2 analysis failed: $(sol.retcode)")
    end
    
    println("✓ ODE solved successfully (Dense output available for last cycle)")
    
    # ========================================================================
    # QuadGK INTEGRATION
    # ========================================================================
    println("\n⏳ Computing integrals with QuadGK...")
    
    results = Dict()
    
    for (i, strategy_idx) in enumerate(strategy_indices)
        offset = 2 + (i - 1) * 8
        strategy_name = "Strategy_$(strategy_idx)"
        
        println("  Processing $strategy_name...")
        
        # Define integrand functions
        
        # Synthesis rates for integration
        function get_synthesis_rates(t)
            u = sol(t)
            X1, M1, X2, M2 = u[offset+1], u[offset+2], u[offset+4], u[offset+5]
            c_minus = u[offset+7]
            
            x1_regulation = (X1 + params.b * params.KX1) / (params.KX1 + X1)
            x2_regulation = (X2 + params.b * params.KX2) / (params.KX2 + X2)
            
            # MODIFICATION: alpha2 -> kM, alpha3 -> kE
            rM1 = params.kM * x1_regulation * c_minus
            rE1 = params.kE * M1 * c_minus
            
            if strategy_idx == 0
                rM2 = params.kM * x2_regulation * c_minus
                rE2 = params.kE * M2 * c_minus
            elseif strategy_idx == 1
                kappa_M = params.K_CCR1 / (params.K_CCR1 + X1)
                rM2 = params.kM * x2_regulation * c_minus * kappa_M
                rE2 = params.kE * M2 * c_minus
            elseif strategy_idx == 2
                kappa_E = params.K_CCR2 / (params.K_CCR2 + X1)
                rM2 = params.kM * x2_regulation * c_minus
                rE2 = params.kE * M2 * c_minus * kappa_E
            else
                kappa_M = params.K_CCR1 / (params.K_CCR1 + X1)
                kappa_E = params.K_CCR2 / (params.K_CCR2 + X1)
                rM2 = params.kM * x2_regulation * c_minus * kappa_M
                rE2 = params.kE * M2 * c_minus * kappa_E
            end
            
            return rM1, rE1, rM2, rE2
        end
        
        # Pathway secific growth rates for integration (Eq. 33-35)
        function get_pathway_rates(t)
            u = sol(t)
            S1, S2 = u[1], u[2]
            X1, E1, X2, E2 = u[offset+1], u[offset+3], u[offset+4], u[offset+6]
            
            # Compute pathway rates
            rs1 = params.k1 * S1 / (params.KS1 + S1) * E1
            rs2 = params.k2 * S2 / (params.KS2 + S2) * E2
            rX1 = params.alpha1 * X1
            rX2 = params.alpha2 * X2
            
            rg1 = rs1 - (1 - params.Y1) * rX1
            rg2 = rs2 - (1 - params.Y2) * rX2
            
            return rs1, rs2, rg1, rg2
        end
        
        # Investment integrals: ∫(rM + rE)·dt where rM is mRNA synthesis rate, rE is enzyme synthesis rate
        # Represents the total resource cost for building pathway machinery (mRNA + enzyme)
        investment1, _ = quadgk(t -> begin
            rM1, rE1, _, _ = get_synthesis_rates(t)
            rM1 + rE1
        end, cycle_start, cycle_end, rtol=1e-8)

        investment2, _ = quadgk(t -> begin
            _, _, rM2, rE2 = get_synthesis_rates(t)
            rM2 + rE2
        end, cycle_start, cycle_end, rtol=1e-8)
        
        # Gain integrals: ∫rg·dt where rg is the specific growth rate contribution
        # rg1 = rs1 - (1-Y1)·α1·X1, represents net biomass gain from S1 pathway
        # rg2 = rs2 - (1-Y2)·α2·X2, represents net biomass gain from S2 pathway
        gain1, _ = quadgk(t -> get_pathway_rates(t)[3], cycle_start, cycle_end, rtol=1e-8)
        gain2, _ = quadgk(t -> get_pathway_rates(t)[4], cycle_start, cycle_end, rtol=1e-8)
        
        # Biomass average: mean biomass concentration over the cycle
        biomass_integral, _ = quadgk(t -> sol(t)[offset+8], cycle_start, cycle_end, rtol=1e-8)
        biomass_avg = biomass_integral / cycle_duration
        
        # Average substrate concentrations over the cycle
        S1_avg, _ = quadgk(t -> sol(t)[1], cycle_start, cycle_end, rtol=1e-8)
        S2_avg, _ = quadgk(t -> sol(t)[2], cycle_start, cycle_end, rtol=1e-8)
        S1_avg /= cycle_duration
        S2_avg /= cycle_duration

        # Total investment: sum of resource costs for both pathways
        total_investment = investment1 + investment2

        # Total gain: sum of specific growth contributions from both pathways
        total_gain = gain1 + gain2

        # Dilution integral: total biomass lost to dilution over the cycle
        dilution_integral = params.D * cycle_duration

        # Net growth: total gain minus dilution loss
        net_growth = total_gain - dilution_integral
        
        results[strategy_name] = Dict(
            # Cycle timing information
            "cycle_start" => cycle_start,
            "cycle_end" => cycle_end,
            "cycle_duration" => cycle_duration,

            # Investment integrals: resource cost for pathway machinery (mRNA + enzyme synthesis)
            "S1_investment_integral" => investment1,
            "S2_investment_integral" => investment2,
            "total_investment_integral" => total_investment,

            # Gain integrals: net biomass gain from each pathway (∫rg·dt)
            "S1_gain_integral" => gain1,
            "S2_gain_integral" => gain2,
            "total_gain_integral" => total_gain,

            # Investment allocation fractions
            "S1_investment_fraction" => total_investment > 1e-10 ? investment1/total_investment : 0.0,
            "S2_investment_fraction" => total_investment > 1e-10 ? investment2/total_investment : 0.0,

            # Gain contribution fractions
            "S1_gain_fraction" => total_gain > 1e-10 ? gain1/total_gain : 0.0,
            "S2_gain_fraction" => total_gain > 1e-10 ? gain2/total_gain : 0.0,

            # Biomass and growth metrics
            "biomass_average" => biomass_avg,
            "dilution_integral" => dilution_integral,
            "net_growth" => net_growth,
            "survival_status" => net_growth >= 0 ? "Surviving" : "Declining",

            # Average substrate concentrations over the cycle
            "S1_average_concentration" => S1_avg,
            "S2_average_concentration" => S2_avg
        )
    end
    
    println("✓ Integration completed")
    println("="^80)
    
    return results, sol, sol_longterm

end


# ============================================================================
# Print results with updated variable names
# ============================================================================
function print_results(results)
    println("\n" * "="^80)
    println("RESULTS SUMMARY")
    println("="^80)

    strategy_names = sort(collect(keys(results)))

    # Individual strategy results
    for strategy_name in strategy_names
        r = results[strategy_name]
        println("\n[$strategy_name]")
        println("-"^40)

        println("INVESTMENT:")
        println("  ∫(rM1 + rE1)dt = $(@sprintf("%.6f", r["S1_investment_integral"]))")
        println("  ∫(rM2 + rE2)dt = $(@sprintf("%.6f", r["S2_investment_integral"]))")
        println("  Total: $(@sprintf("%.6f", r["total_investment_integral"]))")
        println("  Allocation: S1=$(@sprintf("%.1f", r["S1_investment_fraction"]*100))%, " *
                "S2=$(@sprintf("%.1f", r["S2_investment_fraction"]*100))%")

        println("\nGAIN:")
        println("  ∫rg1·dt = $(@sprintf("%.6f", r["S1_gain_integral"]))")
        println("  ∫rg2·dt = $(@sprintf("%.6f", r["S2_gain_integral"]))")
        println("  Total: $(@sprintf("%.6f", r["total_gain_integral"]))")
        println("  Contribution: S1=$(@sprintf("%.1f", r["S1_gain_fraction"]*100))%, " *
                "S2=$(@sprintf("%.1f", r["S2_gain_fraction"]*100))%")

        println("\nGROWTH:")
        println("  Net Growth: $(@sprintf("%.6f", r["net_growth"]))")
        println("  Status: $(r["survival_status"])")
        println("  ∫c·dt/T = $(@sprintf("%.6f", r["biomass_average"]))")

        println("\nSUBSTRATE:")
        println("  ∫S1·dt/T = $(@sprintf("%.6e", r["S1_average_concentration"]))")
        println("  ∫S2·dt/T = $(@sprintf("%.6e", r["S2_average_concentration"]))")
    end

    # Competition analysis
    if length(strategy_names) >= 2
        println("\n" * "="^80)
        println("COMPETITION ANALYSIS")
        println("="^80)

        net_growths = [results[n]["net_growth"] for n in strategy_names]
        winner_idx = argmax(net_growths)
        winner = strategy_names[winner_idx]

        println("\n🏆 WINNER: $winner")
        println("   Net Growth: $(@sprintf("%.6f", maximum(net_growths)))")
        println("   Advantage: $(@sprintf("%.6f", maximum(net_growths) - minimum(net_growths)))")

        println("\nINVESTMENT STRATEGY:")
        for name in strategy_names
            r = results[name]
            println("  $name: $(@sprintf("%.1f", r["S1_investment_fraction"]*100))% → S1, " *
                    "$(@sprintf("%.1f", r["S2_investment_fraction"]*100))% → S2")
        end

        println("\nGAIN CONTRIBUTION:")
        for name in strategy_names
            r = results[name]
            println("  $name: S1=$(@sprintf("%.1f", r["S1_gain_fraction"]*100))%, " *
                    "S2=$(@sprintf("%.1f", r["S2_gain_fraction"]*100))%, " *
                    "Net Growth=$(@sprintf("%.6f", r["net_growth"]))"*
                    "Avg Biomass=$(@sprintf("%.6f", r["biomass_average"]))")
        end
    end

    println("\n" * "="^80)
end


# ============================================================================
# Basic last cycle plotting function
# ============================================================================
function plot_last_cycle(sol, params, strategy_indices, results)
    """Generate plots for the last cycle"""
    
    # Get cycle information
    first_model = "Strategy_$(strategy_indices[1])"
    cycle_start = results[first_model]["cycle_start"]
    cycle_end = results[first_model]["cycle_end"]
    
    # Time array for plotting
    t_plot = range(cycle_start, cycle_end, length=1000)
    
    # Create figure with subplots
    num_strategies = length(strategy_indices)
    
    # Substrate plot
    p1 = plot(title="Substrates", xlabel="Time", ylabel="Concentration",
              legend=:right, size=(800, 400))
    
    S1_vals = [sol(t)[1] for t in t_plot]
    S2_vals = [sol(t)[2] for t in t_plot]
    
    plot!(p1, t_plot, S1_vals, label="S1", linewidth=2)
    plot!(p1, t_plot, S2_vals, label="S2", linewidth=2)
    
    # S1 inlet
    # S1_in uses params.T1 and params.T2 automatically
    S1_in_vals = [S1_in(t, params) for t in t_plot]
    plot!(p1, t_plot, S1_in_vals, label="S1 inlet", linewidth=2, linestyle=:dash, alpha=0.5)
    
    # Biomass plot
    p2 = plot(title="Biomass", xlabel="Time", ylabel="Biomass",
              legend=:right, size=(800, 400))
    
    colors = [:blue, :red, :green, :orange]
    for (i, strategy_idx) in enumerate(strategy_indices)
        offset = 2 + (i - 1) * 8
        c_vals = [sol(t)[offset+8] for t in t_plot]
        plot!(p2, t_plot, c_vals, label="Strategy $strategy_idx", 
              linewidth=2, color=colors[i])
    end
    
    # Pathway variables
    p3 = plot(title="X1 (Pathway 1)", xlabel="Time", ylabel="Concentration",
              legend=:right, size=(800, 400))
    p4 = plot(title="X2 (Pathway 2)", xlabel="Time", ylabel="Concentration",
              legend=:right, size=(800, 400))
    
    for (i, strategy_idx) in enumerate(strategy_indices)
        offset = 2 + (i - 1) * 8
        X1_vals = [sol(t)[offset+1] for t in t_plot]
        X2_vals = [sol(t)[offset+4] for t in t_plot]
        
        plot!(p3, t_plot, X1_vals, label="Strategy $strategy_idx", 
              linewidth=2, color=colors[i])
        plot!(p4, t_plot, X2_vals, label="Strategy $strategy_idx", 
              linewidth=2, color=colors[i])
    end
    
    # Combine all plots
    p_combined = plot(p1, p2, p3, p4, layout=(2, 2), size=(1600, 1200))
    
    return p_combined
end





# ============================================================================
# HELPER FUNCTION: Add subplot label to top-left corner (outside axis area)
# ============================================================================
"""
Add a subplot label (e.g., "(a)", "(b)") to the top-left corner of a subplot,
positioned outside the axis area (above and to the left of the plot).

Parameters:
- p: The subplot to add the label to
- label: The label string (e.g., "(a)")
- fontsize: Font size for the label (default: 28)
- x_offset: Horizontal offset as fraction of x-range, negative = left of axis (default: -0.08)
- y_offset: Vertical offset as fraction of y-range, >1 = above axis (default: 1.08)
"""
function add_subplot_label!(p, label; fontsize=28, x_offset=-0.08, y_offset=1.1)
    xl = xlims(p)
    yl = ylims(p)
    x_pos = xl[1] + x_offset * (xl[2] - xl[1])
    y_pos = yl[1] + y_offset * (yl[2] - yl[1])
    annotate!(p, [(x_pos, y_pos, text(label, :left, :bottom, fontsize, :bold))])
end

# ============================================================================
# BIOMASS PLOTTING FUNCTION (Long-term competition dynamics)
# ============================================================================

# ============================================================================
# P1: Biomass plot (long-term competition dynamics until steady state)
# ============================================================================

function plot_p1_biomass(sol_longterm, params, strategy_indices;
                         linewidth=3,
                         titlefontsize=22,
                         xlabelfontsize=20,
                         ylabelfontsize=20,
                         tickfontsize=18,
                         legendfontsize=18,
                         title_gap=2)  # Number of newlines after subplot title (controls title-plot spacing)

    # Generate title gap string
    gap_str = repeat("\n", title_gap)

    # Model names and colors (simplified for legend)
    all_strategy_names = ["N", "T", "P"]
    all_strategy_colors = [RGB(0.95,0.8,0.25),
                       RGB(0.35,0.8,0.9),
                       RGB(0.95,0.6,0.5),
                       ]

    strategy_names = [all_strategy_names[i + 1] for i in strategy_indices]
    strategy_colors = [all_strategy_colors[i + 1] for i in strategy_indices]

    # Generate title with strategy comparison
    strategy_comparison = join(strategy_names, " vs ")
    main_title = "Biomass Dynamics During Competition ($strategy_comparison)"

    t_eval = sol_longterm.t

    println("📊 Creating P1: Biomass plots...")

    # Size: 204mm × 70mm @ 300dpi ≈ 2409 × 827 pixels
    p = plot(layout=(1,2), size=(1772, 662),
             left_margin=15Plots.mm, right_margin=10Plots.mm,
             bottom_margin=15Plots.mm,
                top_margin=15Plots.mm,
             #plot_title=main_title, plot_titlefontsize=titlefontsize+2
             )

    for (i, strategy_idx) in enumerate(strategy_indices)
        offset = 2 + (i - 1) * 8
        biomass_idx = offset + 8
        biomass = [sol_longterm[biomass_idx, j] for j in 1:length(t_eval)]

        # Subplot 1: Linear scale
        plot!(p[1], t_eval, biomass,
              color=strategy_colors[i], linewidth=linewidth, label=strategy_names[i],
              title="c (biomass) - Linear Scale" * gap_str, titlefontsize=titlefontsize,
              xlabel="Time", xlabelfontsize=xlabelfontsize,
              ylabel="Biomass Concentration", ylabelfontsize=ylabelfontsize,
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              grid=false, legendfontsize=legendfontsize, legend=:right)

        # Subplot 2: Log scale
        safe_biomass = copy(biomass)
        safe_biomass[safe_biomass .< 1e-10] .= NaN
        plot!(p[2], t_eval, safe_biomass,
              color=strategy_colors[i], linewidth=linewidth, label=strategy_names[i],
              title="c (biomass) - Log Scale" * gap_str, titlefontsize=titlefontsize,
              xlabel="Time", xlabelfontsize=xlabelfontsize,
              ylabel="Biomass Concentration (log₁₀)", ylabelfontsize=ylabelfontsize,
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              legendfontsize=legendfontsize, yscale=:log10, grid=false, legend=:right)
    end

    # Add subplot labels to top-left corner (outside axis area)
    add_subplot_label!(p[1], "(a)", fontsize=titlefontsize, x_offset=-0.07, y_offset=1.01)
    add_subplot_label!(p[2], "(b)", fontsize=titlefontsize)

    println("✅ P1 created!")
    return p
end


# ============================================================================
# P2: Variables in a Single Period at Steady State (4×3 layout)
# ============================================================================
function plot_p2_variables(sol, params, strategy_indices;
                           num_points=1000,
                           linewidth=3,
                           titlefontsize=22,
                           xlabelfontsize=20,
                           ylabelfontsize=20,
                           tickfontsize=18,
                           legendfontsize=18,
                           title_gap=1,  # Number of newlines after subplot title
                           xtick_count=4)  # Number of x-axis ticks (nothing = auto)

    # Generate title gap string
    gap_str = repeat("\n", title_gap)

    # Model names and colors (simplified for legend)
    all_strategy_names = ["N", "T", "P"]
    all_strategy_colors = [RGB(0.95,0.8,0.25),
                       RGB(0.35,0.8,0.9),
                       RGB(0.95,0.6,0.5),
                       ]

    strategy_names = [all_strategy_names[i + 1] for i in strategy_indices]
    strategy_colors = [all_strategy_colors[i + 1] for i in strategy_indices]

    # Generate title with strategy comparison
    strategy_comparison = join(strategy_names, " vs ")
    main_title = "Variables in a Single Period at Steady State ($strategy_comparison)"

    # Time points from dense output
    cycle_start = sol.t[1]
    cycle_end = sol.t[end]
    t_eval = range(cycle_start, cycle_end, length=num_points)
    t_rel = t_eval .- cycle_start  # Relative time starting from 0

    println("📊 Creating P2: Variable plots (4×3 layout)...")

    # ========== Step 1: Extract all data ==========

    # S1_in, S2_in data (substrate input)
    S1_in_data = [S1_in(t, params) for t in t_eval]
    S2_in_data = [S2_in(t, params) for t in t_eval]

    # S1, S2 data (substrate in chemostat)
    S1_data = [sol(t)[1] for t in t_eval]
    S2_data = [sol(t)[2] for t in t_eval]

    # Extract all strategy variables
    all_var_data = Dict()
    for (i, strategy_idx) in enumerate(strategy_indices)
        offset = 2 + (i - 1) * 8
        for v in 1:8
            var_idx = offset + v
            data = [sol(t)[var_idx] for t in t_eval]
            all_var_data[(i, v)] = data
        end
    end

    # ========== Step 2: Calculate Y-axis limits ==========

    # S_in Y-limit
    S_in_min = min(minimum(S1_in_data), minimum(S2_in_data))
    S_in_max = max(maximum(S1_in_data), maximum(S2_in_data))
    S_in_range = S_in_max - S_in_min
    S_in_ylim = (max(0, S_in_min - 0.1*S_in_range), S_in_max + 0.1*S_in_range)

    # S Y-limit
    S_min = min(minimum(S1_data), minimum(S2_data))
    S_max = max(maximum(S1_data), maximum(S2_data))
    S_range = S_max - S_min
    S_ylim = (max(0, S_min - 0.1*S_range), S_max + 0.1*S_range)

    # c (biomass) Y-limit (v=8)
    c_min, c_max = Inf, -Inf
    for (i, _) in enumerate(strategy_indices)
        data = all_var_data[(i, 8)]
        c_min = min(c_min, minimum(data))
        c_max = max(c_max, maximum(data))
    end
    c_range = c_max - c_min
    c_ylim = (max(0, c_min - 0.1*c_range), c_max + 0.1*c_range)

    # Shared Y-limit for X1, M1, E1, X2, M2, E2, c- (v=1,2,3,4,5,6,7)
    shared_var_indices = [1, 2, 3, 4, 5, 6, 7]
    shared_min, shared_max = Inf, -Inf
    for (i, _) in enumerate(strategy_indices)
        for v in shared_var_indices
            data = all_var_data[(i, v)]
            shared_min = min(shared_min, minimum(data))
            shared_max = max(shared_max, maximum(data))
        end
    end
    shared_range = shared_max - shared_min
    shared_ylim = (max(0, shared_min - 0.1*shared_range), shared_max + 0.1*shared_range)

    # ========== Step 3: Create plot (4×3 layout) ==========

    p = plot(layout=(4, 3), size=(1772, 1391),
             left_margin=12Plots.mm, right_margin=8Plots.mm,
             bottom_margin=10Plots.mm,
             )

    # === Row 1: S_in, S, c ===

    # Position 1: S_in (substrate input concentration)
    plot!(p[1], t_rel, S1_in_data, color=:red, linewidth=linewidth, label="S₁,in",
          title="Carbon source (input) " * gap_str, titlefontsize=titlefontsize,
          xlabel="", ylabel="Concentration", ylabelfontsize=ylabelfontsize,
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=S_in_ylim, grid=false, legend=:topright, legendfontsize=legendfontsize)
    plot!(p[1], t_rel, S2_in_data, color=:blue, linewidth=linewidth, label="S₂,in")
    add_phase_background!(p[1], 0.0, params)

    # Position 2: S (substrate in chemostat) - use S_in_ylim to match (a)
    plot!(p[2], t_rel, S1_data, color=:red, linewidth=linewidth, label="S₁",
          title="Carbon source in chemostat" * gap_str, titlefontsize=titlefontsize,
          xlabel="", ylabel="",
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=S_in_ylim, grid=false, legend=:topright, legendfontsize=legendfontsize)
    plot!(p[2], t_rel, S2_data, color=:blue, linewidth=linewidth, label="S₂")
    add_phase_background!(p[2], 0.0, params)

    # Position 3: c (biomass)
    for (i, strategy_idx) in enumerate(strategy_indices)
        plot!(p[3], t_rel, all_var_data[(i, 8)],
              color=strategy_colors[i], linewidth=linewidth, label=strategy_names[i],
              title="c (biomass)" * gap_str, titlefontsize=titlefontsize,
              xlabel="", ylabel="",
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              ylims=c_ylim, grid=false,
              legend=(i == length(strategy_indices) ? :topright : false),
              legendfontsize=legendfontsize)
    end
    add_phase_background!(p[3], 0.0, params)

    # === Row 2: X1, M1, E1 ===

    # Position 4: X1 (internal S1)
    for (i, strategy_idx) in enumerate(strategy_indices)
        plot!(p[4], t_rel, all_var_data[(i, 1)],
              color=strategy_colors[i], linewidth=linewidth, label=strategy_names[i],
              title="X₁ (internal S₁)" * gap_str, titlefontsize=titlefontsize,
              xlabel="", ylabel="Concentration", ylabelfontsize=ylabelfontsize,
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              ylims=shared_ylim, grid=false, legend=false)
    end
    add_phase_background!(p[4], 0.0, params)

    # Position 5: M1 (S1 mRNA)
    for (i, strategy_idx) in enumerate(strategy_indices)
        plot!(p[5], t_rel, all_var_data[(i, 2)],
              color=strategy_colors[i], linewidth=linewidth, label=strategy_names[i],
              title="M₁ (S₁ mRNA)" * gap_str, titlefontsize=titlefontsize,
              xlabel="", ylabel="",
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              ylims=shared_ylim, grid=false, legend=false)
    end
    add_phase_background!(p[5], 0.0, params)

    # Position 6: E1 (S1 enzyme)
    for (i, strategy_idx) in enumerate(strategy_indices)
        plot!(p[6], t_rel, all_var_data[(i, 3)],
              color=strategy_colors[i], linewidth=linewidth, label=strategy_names[i],
              title="E₁ (S₁ enzyme)" * gap_str, titlefontsize=titlefontsize,
              xlabel="", ylabel="",
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              ylims=shared_ylim, grid=false, legend=false)
    end
    add_phase_background!(p[6], 0.0, params)

    # === Row 3: X2, M2, E2 ===

    # Position 7: X2 (internal S2)
    for (i, strategy_idx) in enumerate(strategy_indices)
        plot!(p[7], t_rel, all_var_data[(i, 4)],
              color=strategy_colors[i], linewidth=linewidth, label=strategy_names[i],
              title="X₂ (internal S₂)" * gap_str, titlefontsize=titlefontsize,
              xlabel="", ylabel="Concentration", ylabelfontsize=ylabelfontsize,
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              ylims=shared_ylim, grid=false, legend=false)
    end
    add_phase_background!(p[7], 0.0, params)

    # Position 8: M2 (S2 mRNA)
    for (i, strategy_idx) in enumerate(strategy_indices)
        plot!(p[8], t_rel, all_var_data[(i, 5)],
              color=strategy_colors[i], linewidth=linewidth, label=strategy_names[i],
              title="M₂ (S₂ mRNA)" * gap_str, titlefontsize=titlefontsize,
              xlabel="", ylabel="",
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              ylims=shared_ylim, grid=false, legend=false)
    end
    add_phase_background!(p[8], 0.0, params)

    # Position 9: E2 (S2 enzyme)
    for (i, strategy_idx) in enumerate(strategy_indices)
        plot!(p[9], t_rel, all_var_data[(i, 6)],
              color=strategy_colors[i], linewidth=linewidth, label=strategy_names[i],
              title="E₂ (S₂ enzyme)" * gap_str, titlefontsize=titlefontsize,
              xlabel="", ylabel="",
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              ylims=shared_ylim, grid=false, legend=false)
    end
    add_phase_background!(p[9], 0.0, params)

    # === Row 4: c-, empty, empty ===

    # Position 10: c- (intracellular resource)
    for (i, strategy_idx) in enumerate(strategy_indices)
        plot!(p[10], t_rel, all_var_data[(i, 7)],
              color=strategy_colors[i], linewidth=linewidth, label=strategy_names[i],
              title="c⁻ (intracell. resource)" * gap_str, titlefontsize=titlefontsize,
              xlabel="Time", xlabelfontsize=xlabelfontsize,
              ylabel="Concentration", ylabelfontsize=ylabelfontsize,
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              ylims=shared_ylim, grid=false, legend=false)
    end
    add_phase_background!(p[10], 0.0, params)

    # Position 11, 12: empty
    plot!(p[11], framestyle=:none, grid=false, showaxis=false, ticks=false)
    plot!(p[12], framestyle=:none, grid=false, showaxis=false, ticks=false)

    # Add subplot labels to top-left corner (outside axis area)
    subplot_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)"]
    for i in 1:10
        add_subplot_label!(p[i], subplot_labels[i], fontsize=titlefontsize, x_offset=-0.25, y_offset=1.1)
    end

    println("✅ P2 created!")
    println("   S_in Y-lim: $S_in_ylim")
    println("   S Y-lim: $S_ylim")
    println("   c Y-lim: $c_ylim")
    println("   Shared (X1,M1,E1,X2,M2,E2,c-) Y-lim: $shared_ylim")
    return p
end















# ============================================================================
# P3: Variable Differences in a Single Period at Steady State (4×3 layout)
# ============================================================================
function plot_p3_variable_differences(sol, params, strategy_indices;
                                      num_points=1000,
                                      linewidth=3,
                                      titlefontsize=22,
                                      xlabelfontsize=20,
                                      ylabelfontsize=20,
                                      tickfontsize=18,
                                      legendfontsize=18,
                                      title_gap=1,  # Number of newlines after subplot title
                                      xtick_count=4)  # Number of x-axis ticks (nothing = auto)
    if length(strategy_indices) != 2
        println("⚠️ P3 requires exactly 2 strategies. Skipping.")
        return nothing
    end

    # Generate title gap string
    gap_str = repeat("\n", title_gap)

    # Model names and colors
    all_strategy_names = ["N", "T", "P"]
    all_strategy_colors = [RGB(0.95,0.8,0.25),
                       RGB(0.35,0.8,0.9),
                       RGB(0.95,0.6,0.5),
                       ]

    strategy_names = [all_strategy_names[i + 1] for i in strategy_indices]
    strategy_colors = [all_strategy_colors[i + 1] for i in strategy_indices]
    diff_label = "$(strategy_names[1]) - $(strategy_names[2])"

    # Generate title with strategy comparison
    strategy_comparison = join(strategy_names, " vs ")
    main_title = "Variable Differences in a Single Period at Steady State ($strategy_comparison)"

    cycle_start = sol.t[1]
    cycle_end = sol.t[end]
    t_eval = collect(range(cycle_start, cycle_end, length=num_points))
    t_rel = t_eval .- cycle_start  # Relative time starting from 0

    println("📊 Creating P3: Variable difference plots (4×3 layout)...")

    # ========== Step 1: Extract data ==========

    # S1_in, S2_in data (substrate input) - same as P2
    S1_in_data = [S1_in(t, params) for t in t_eval]
    S2_in_data = [S2_in(t, params) for t in t_eval]

    # S1, S2 data (substrate in chemostat) - same as P2
    S1_data = [sol(t)[1] for t in t_eval]
    S2_data = [sol(t)[2] for t in t_eval]

    # S_in Y-limit
    S_in_min = min(minimum(S1_in_data), minimum(S2_in_data))
    S_in_max = max(maximum(S1_in_data), maximum(S2_in_data))
    S_in_range = S_in_max - S_in_min
    S_in_ylim = (max(0, S_in_min - 0.1*S_in_range), S_in_max + 0.1*S_in_range)

    # S Y-limit
    S_min = min(minimum(S1_data), minimum(S2_data))
    S_max = max(maximum(S1_data), maximum(S2_data))
    S_range = S_max - S_min
    S_ylim = (max(0, S_min - 0.1*S_range), S_max + 0.1*S_range)

    # Calculate variable differences
    offset1 = 2 + (1 - 1) * 8  # Model 1 offset = 2
    offset2 = 2 + (2 - 1) * 8  # Model 2 offset = 10

    diff_data = Dict()
    for v in 1:8
        var_idx1 = offset1 + v
        var_idx2 = offset2 + v
        data1 = [sol(t)[var_idx1] for t in t_eval]
        data2 = [sol(t)[var_idx2] for t in t_eval]
        diff_data[v] = data1 .- data2
    end

    # Calculate Δ(M2+E2)
    M2_data1 = [sol(t)[offset1 + 5] for t in t_eval]
    E2_data1 = [sol(t)[offset1 + 6] for t in t_eval]
    M2_data2 = [sol(t)[offset2 + 5] for t in t_eval]
    E2_data2 = [sol(t)[offset2 + 6] for t in t_eval]
    saved_diff = (M2_data1 .+ E2_data1) .- (M2_data2 .+ E2_data2)

    # ========== Step 2: Calculate Y-axis limits ==========

    # Δc Y-limit (v=8, independent)
    c_diff_min = minimum(diff_data[8])
    c_diff_max = maximum(diff_data[8])
    c_diff_abs_max = max(abs(c_diff_min), abs(c_diff_max))
    c_diff_ylim = (-c_diff_abs_max * 1.1, c_diff_abs_max * 1.1)

    # Shared diff Y-limit for ΔX1, ΔM1, ΔE1, ΔX2, ΔM2, ΔE2, Δc-, Δ(M2+E2)
    shared_diff_indices = [1, 2, 3, 4, 5, 6, 7]
    shared_diff_min, shared_diff_max = Inf, -Inf
    for v in shared_diff_indices
        shared_diff_min = min(shared_diff_min, minimum(diff_data[v]))
        shared_diff_max = max(shared_diff_max, maximum(diff_data[v]))
    end
    shared_diff_min = min(shared_diff_min, minimum(saved_diff))
    shared_diff_max = max(shared_diff_max, maximum(saved_diff))

    shared_diff_abs_max = max(abs(shared_diff_min), abs(shared_diff_max))
    shared_diff_ylim = (-shared_diff_abs_max * 1.1, shared_diff_abs_max * 1.1)

    # ========== Step 3: Create plot (4×3 layout) ==========
    p = plot(layout=(4, 3), size=(1772, 1391),
             left_margin=12Plots.mm, right_margin=8Plots.mm,
             bottom_margin=10Plots.mm,
             )

    # === Row 1: S_in, S, Δc ===

    # Manual ylabel x-offset for row 1 (adjust this value to align with "Difference" labels below)
    ylabel_x_offset = -0.44
    #ylabel_x_offset = -0.375

    # Position 1: S_in (substrate input concentration) - same as P2
    plot!(p[1], t_rel, S1_in_data, color=:red, linewidth=linewidth, label="S₁,in",
         title="Carbon source (input)" * gap_str, titlefontsize=titlefontsize,
          xlabel="", ylabel="", ylabelfontsize=ylabelfontsize,
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=S_in_ylim, grid=false, legend=:topright, legendfontsize=legendfontsize)
    plot!(p[1], t_rel, S2_in_data, color=:blue, linewidth=linewidth, label="S₂,in")
    add_phase_background!(p[1], 0.0, params)
    # Manual ylabel for alignment
    xl1 = xlims(p[1]); yl1 = ylims(p[1])
    annotate!(p[1], xl1[1] + ylabel_x_offset * (xl1[2] - xl1[1]), (yl1[1] + yl1[2]) / 2,
              text("Concentration", :center, ylabelfontsize, rotation=90))

    # Position 2: S (substrate in chemostat) - use S_in_ylim to match (a)
    plot!(p[2], t_rel, S1_data, color=:red, linewidth=linewidth, label="S₁",
          title="Carbon source in chemostat" * gap_str, titlefontsize=titlefontsize,
          xlabel="", ylabel="",
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=S_in_ylim, grid=false, legend=:topright, legendfontsize=legendfontsize)
    plot!(p[2], t_rel, S2_data, color=:blue, linewidth=linewidth, label="S₂")
    add_phase_background!(p[2], 0.0, params)

    # Position 3: Δc (biomass)
    # Layer order: 0.set ylims → 1.T1/T2 background → 2.white base → 3.strategy color → 4.curve
    ylims!(p[3], c_diff_ylim)
    add_phase_background!(p[3], 0.0, params)
    add_difference_background!(p[3], t_rel, diff_data[8], [:white, :white], alpha=1.0)
    add_difference_background!(p[3], t_rel, diff_data[8], strategy_colors, alpha=1.0)
    plot!(p[3], t_rel, diff_data[8], color=:black, linewidth=linewidth, label="",
          title="Δc (biomass)" * gap_str, titlefontsize=titlefontsize,
          xlabel="",  ylabel="Difference", ylabelfontsize=ylabelfontsize,
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=c_diff_ylim, grid=false, legend=:topright, legendfontsize=legendfontsize)
    plot!(p[3], [NaN], [NaN], color=strategy_colors[1], linewidth=linewidth, label=strategy_names[1])
    plot!(p[3], [NaN], [NaN], color=strategy_colors[2], linewidth=linewidth, label=strategy_names[2])

    # === Row 2: ΔX1, ΔM1, ΔE1 ===

    # Position 4: ΔX1 (internal S1)
    ylims!(p[4], shared_diff_ylim)
    add_phase_background!(p[4], 0.0, params)
    add_difference_background!(p[4], t_rel, diff_data[1], [:white, :white], alpha=1.0)
    add_difference_background!(p[4], t_rel, diff_data[1], strategy_colors, alpha=1.0)
    plot!(p[4], t_rel, diff_data[1], color=:black, linewidth=linewidth, label=diff_label,
          title="ΔX₁ (internal S₁)" * gap_str, titlefontsize=titlefontsize,
          xlabel="", ylabel="Difference", ylabelfontsize=ylabelfontsize,
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=shared_diff_ylim, grid=false, legend=false)

    # Position 5: ΔM1 (S1 mRNA)
    ylims!(p[5], shared_diff_ylim)
    add_phase_background!(p[5], 0.0, params)
    add_difference_background!(p[5], t_rel, diff_data[2], [:white, :white], alpha=1.0)
    add_difference_background!(p[5], t_rel, diff_data[2], strategy_colors, alpha=1.0)
    plot!(p[5], t_rel, diff_data[2], color=:black, linewidth=linewidth, label=diff_label,
          title="ΔM₁ (S₁ mRNA)" * gap_str, titlefontsize=titlefontsize,
          xlabel="", ylabel="",
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=shared_diff_ylim, grid=false, legend=false)

    # Position 6: ΔE1 (S1 enzyme)
    ylims!(p[6], shared_diff_ylim)
    add_phase_background!(p[6], 0.0, params)
    add_difference_background!(p[6], t_rel, diff_data[3], [:white, :white], alpha=1.0)
    add_difference_background!(p[6], t_rel, diff_data[3], strategy_colors, alpha=1.0)
    plot!(p[6], t_rel, diff_data[3], color=:black, linewidth=linewidth, label=diff_label,
          title="ΔE₁ (S₁ enzyme)" * gap_str, titlefontsize=titlefontsize,
          xlabel="", ylabel="",
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=shared_diff_ylim, grid=false, legend=false)

    # === Row 3: ΔX2, ΔM2, ΔE2 ===

    # Position 7: ΔX2 (internal S2)
    ylims!(p[7], shared_diff_ylim)
    add_phase_background!(p[7], 0.0, params)
    add_difference_background!(p[7], t_rel, diff_data[4], [:white, :white], alpha=1.0)
    add_difference_background!(p[7], t_rel, diff_data[4], strategy_colors, alpha=1.0)
    plot!(p[7], t_rel, diff_data[4], color=:black, linewidth=linewidth, label=diff_label,
          title="ΔX₂ (internal S₂)" * gap_str, titlefontsize=titlefontsize,
          xlabel="", ylabel="Difference", ylabelfontsize=ylabelfontsize,
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=shared_diff_ylim, grid=false, legend=false)

    # Position 8: ΔM2 (S2 mRNA)
    ylims!(p[8], shared_diff_ylim)
    add_phase_background!(p[8], 0.0, params)
    add_difference_background!(p[8], t_rel, diff_data[5], [:white, :white], alpha=1.0)
    add_difference_background!(p[8], t_rel, diff_data[5], strategy_colors, alpha=1.0)
    plot!(p[8], t_rel, diff_data[5], color=:black, linewidth=linewidth, label=diff_label,
          title="ΔM₂ (S₂ mRNA)" * gap_str, titlefontsize=titlefontsize,
          xlabel="", ylabel="",
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=shared_diff_ylim, grid=false, legend=false)

    # Position 9: ΔE2 (S2 enzyme)
    ylims!(p[9], shared_diff_ylim)
    add_phase_background!(p[9], 0.0, params)
    add_difference_background!(p[9], t_rel, diff_data[6], [:white, :white], alpha=1.0)
    add_difference_background!(p[9], t_rel, diff_data[6], strategy_colors, alpha=1.0)
    plot!(p[9], t_rel, diff_data[6], color=:black, linewidth=linewidth, label=diff_label,
          title="ΔE₂ (S₂ enzyme)" * gap_str, titlefontsize=titlefontsize,
          xlabel="", ylabel="",
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=shared_diff_ylim, grid=false, legend=false)

    # === Row 4: Δc-, Δ(M2+E2), empty ===

    # Position 10: Δc- (intracellular resource)
    ylims!(p[10], shared_diff_ylim)
    add_phase_background!(p[10], 0.0, params)
    add_difference_background!(p[10], t_rel, diff_data[7], [:white, :white], alpha=1.0)
    add_difference_background!(p[10], t_rel, diff_data[7], strategy_colors, alpha=1.0)
    plot!(p[10], t_rel, diff_data[7], color=:black, linewidth=linewidth, label=diff_label,
          title="Δc⁻ (intracell. resource)" * gap_str, titlefontsize=titlefontsize,
          xlabel="Time", xlabelfontsize=xlabelfontsize,
          ylabel="Difference", ylabelfontsize=ylabelfontsize,
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=shared_diff_ylim, grid=false, legend=false)

    # Position 11: Δ(M2+E2) (S2 expression)
    ylims!(p[11], shared_diff_ylim)
    add_phase_background!(p[11], 0.0, params)
    add_difference_background!(p[11], t_rel, saved_diff, [:white, :white], alpha=1.0)
    add_difference_background!(p[11], t_rel, saved_diff, strategy_colors, alpha=1.0)
    plot!(p[11], t_rel, saved_diff, color=:black, linewidth=linewidth, label=diff_label,
          title="Δ(M₂+E₂) (S₂ expression)" * gap_str, titlefontsize=titlefontsize,
          xlabel="Time", xlabelfontsize=xlabelfontsize, ylabel="",
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=shared_diff_ylim, grid=false, legend=false)

    # Position 12: empty
    plot!(p[12], framestyle=:none, grid=false, showaxis=false, ticks=false)

    # Add subplot labels to top-left corner (outside axis area)
    subplot_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)"]
    for i in 1:11
        add_subplot_label!(p[i], subplot_labels[i], fontsize=titlefontsize, x_offset=-0.35, y_offset=1.1)
    end

    println("✅ P3 created!")
    println("   S_in Y-lim: $S_in_ylim")
    println("   S Y-lim: $S_ylim")
    println("   Δc Y-lim: $c_diff_ylim")
    println("   Shared diff Y-lim: $shared_diff_ylim")
    return p
end



# ============================================================================
# P4: Specific Growth Rates in a Single Period at Steady State
# ============================================================================
function plot_p4_growth_rates(sol, params, strategy_indices;
                              num_points=1000,
                              linewidth=3,
                              titlefontsize=22,
                              xlabelfontsize=20,
                              ylabelfontsize=20,
                              tickfontsize=18,
                              legendfontsize=18,
                              title_gap=1,  # Number of newlines after subplot title
                              xtick_count=4)  # Number of x-axis ticks (nothing = auto)
    # Generate title gap string
    gap_str = repeat("\n", title_gap)

    # Model names and colors
    all_strategy_names = ["N", "T", "P"]
    all_strategy_colors = [RGB(0.95,0.8,0.25),
                       RGB(0.35,0.8,0.9),
                       RGB(0.95,0.6,0.5),
                       ]

    strategy_names = [all_strategy_names[i + 1] for i in strategy_indices]
    strategy_colors = [all_strategy_colors[i + 1] for i in strategy_indices]

    # Generate title with strategy comparison
    strategy_comparison = join(strategy_names, " vs ")
    main_title = "Specific Growth Rates in a Single Period at Steady State ($strategy_comparison)"

    cycle_start = sol.t[1]
    cycle_end = sol.t[end]
    t_eval = collect(range(cycle_start, cycle_end, length=num_points))
    t_rel = t_eval .- cycle_start  # Relative time starting from 0

    println("📊 Creating P4: Specific growth rate plots...")

    # ========== Step 1: Calculate specific growth rates ==========
    rg1_data = Dict()
    rg2_data = Dict()
    rg_data = Dict()

    rate_min, rate_max = Inf, -Inf

    for (i, strategy_idx) in enumerate(strategy_indices)
        offset = 2 + (i - 1) * 8

        rg1_vals = Float64[]
        rg2_vals = Float64[]

        for t in t_eval
            u = sol(t)
            S1, S2 = u[1], u[2]
            X1 = u[offset + 1]
            X2 = u[offset + 4]
            E1 = u[offset + 3]
            E2 = u[offset + 6]

            # S1 pathway specific growth rate (consistent with ODE model and QuadGK)
            rs1 = params.k1 * E1 * S1 / (params.KS1 + S1)
            rX1 = params.alpha1 * X1
            rg1 = rs1 - (1 - params.Y1) * rX1

            # S2 pathway specific growth rate (consistent with ODE model and QuadGK)
            rs2 = params.k2 * E2 * S2 / (params.KS2 + S2)
            rX2 = params.alpha2 * X2
            rg2 = rs2 - (1 - params.Y2) * rX2

            push!(rg1_vals, rg1)
            push!(rg2_vals, rg2)
        end

        rg_vals = rg1_vals .+ rg2_vals

        rg1_data[i] = rg1_vals
        rg2_data[i] = rg2_vals
        rg_data[i] = rg_vals

        rate_min = min(rate_min, minimum(rg1_vals), minimum(rg2_vals), minimum(rg_vals))
        rate_max = max(rate_max, maximum(rg1_vals), maximum(rg2_vals), maximum(rg_vals))
    end

    # Shared Y-axis range
    rate_range = rate_max - rate_min
    rate_ylim = (max(0, rate_min - 0.1*rate_range), rate_max + 0.1*rate_range)

    # ========== Step 2: Create plot (3×1 layout) ==========
    p = plot(layout=(3, 1), size=(1772, 1309),
             left_margin=15Plots.mm, right_margin=10Plots.mm,
             bottom_margin=12Plots.mm,
             )

    # Subplot 1: Total specific growth rate rg
    for (i, strategy_idx) in enumerate(strategy_indices)
        plot!(p[1], t_rel, rg_data[i],
              color=strategy_colors[i], linewidth=linewidth, label=strategy_names[i],
              title="rg (Total Specific Growth Rate = rg₁ + rg₂)" * gap_str, titlefontsize=titlefontsize,
              xlabel="", ylabel="Specific Growth Rate", ylabelfontsize=ylabelfontsize,
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              ylims=rate_ylim, grid=false, legend=false)
    end
    hline!(p[1], [params.D], color=:gray, linestyle=:dash, linewidth=2, label="")
    add_phase_background!(p[1], 0.0, params)

    # Subplot 2: S1 pathway specific growth rate rg1
    for (i, strategy_idx) in enumerate(strategy_indices)
        plot!(p[2], t_rel, rg1_data[i],
              color=strategy_colors[i], linewidth=linewidth, label=strategy_names[i],
              title="rg₁ (S₁ Pathway Specific Growth Rate)" * gap_str, titlefontsize=titlefontsize,
              xlabel="", ylabel="Specific Growth Rate", ylabelfontsize=ylabelfontsize,
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              ylims=rate_ylim, grid=false, legend=false)
    end
    add_phase_background!(p[2], 0.0, params)

    # Subplot 3: S2 pathway specific growth rate rg2
    for (i, strategy_idx) in enumerate(strategy_indices)
        plot!(p[3], t_rel, rg2_data[i],
              color=strategy_colors[i], linewidth=linewidth, label=strategy_names[i],
              title="rg₂ (S₂ Pathway Specific Growth Rate)" * gap_str, titlefontsize=titlefontsize,
              xlabel="Time", xlabelfontsize=xlabelfontsize,
              ylabel="Specific Growth Rate", ylabelfontsize=ylabelfontsize,
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              ylims=rate_ylim, grid=false, legend=false)
    end
    add_phase_background!(p[3], 0.0, params)

    # Add subplot labels to top-left corner (outside axis area)
    add_subplot_label!(p[1], "(a)", fontsize=titlefontsize, y_offset=1.1, x_offset=-0.05)
    add_subplot_label!(p[2], "(b)", fontsize=titlefontsize, y_offset=1.1, x_offset=-0.05)
    add_subplot_label!(p[3], "(c)", fontsize=titlefontsize, y_offset=1.1, x_offset=-0.05)

    println("✅ P4 created!")
    println("   Specific growth rate Y-lim: $rate_ylim")
    return p, rg_data, rg1_data, rg2_data, t_eval
end


# ============================================================================
# P5: Specifc growth Rate Differences in a Single Period at Steady State
# ============================================================================
function plot_p5_growth_rate_differences(sol, params, strategy_indices;
                                         num_points=1000,
                                         linewidth=3,
                                         titlefontsize=22,
                                         xlabelfontsize=20,
                                         ylabelfontsize=20,
                                         tickfontsize=18,
                                         legendfontsize=18,
                                         title_gap=1,  # Number of newlines after subplot title
                                         xtick_count=4)  # Number of x-axis ticks (nothing = auto)
    if length(strategy_indices) != 2
        println("⚠️ P5 requires exactly 2 strategies. Skipping.")
        return nothing
    end

    # Generate title gap string
    gap_str = repeat("\n", title_gap)

    # Model names and colors
    all_strategy_names = ["N", "T", "P"]
    all_strategy_colors = [RGB(0.95,0.8,0.25),
                       RGB(0.35,0.8,0.9),
                       RGB(0.95,0.6,0.5),
                       ]

    strategy_names = [all_strategy_names[i + 1] for i in strategy_indices]
    strategy_colors = [all_strategy_colors[i + 1] for i in strategy_indices]
    diff_label = "$(strategy_names[1]) - $(strategy_names[2])"

    # Generate title with strategy comparison
    strategy_comparison = join(strategy_names, " vs ")
    main_title = "Specific Growth Rate Differences in a Single Period at Steady State ($strategy_comparison)"

    cycle_start = sol.t[1]
    cycle_end = sol.t[end]
    t_eval = collect(range(cycle_start, cycle_end, length=num_points))
    t_rel = t_eval .- cycle_start  # Relative time starting from 0

    println("📊 Creating P5: Specifc growth rate difference plots...")

    # ========== Step 1: Calculate specific growth rates ==========
    rg1_data = [Float64[] for _ in 1:2]
    rg2_data = [Float64[] for _ in 1:2]

    for (i, strategy_idx) in enumerate(strategy_indices)
        offset = 2 + (i - 1) * 8

        for t in t_eval
            u = sol(t)
            S1, S2 = u[1], u[2]
            X1 = u[offset + 1]
            X2 = u[offset + 4]
            E1 = u[offset + 3]
            E2 = u[offset + 6]

            # S1 pathway specific growth rate (consistent with ODE model and QuadGK)
            rs1 = params.k1 * E1 * S1 / (params.KS1 + S1)
            rX1 = params.alpha1 * X1
            rg1 = rs1 - (1 - params.Y1) * rX1

            # S2 pathway specific growth rate (consistent with ODE model and QuadGK)
            rs2 = params.k2 * E2 * S2 / (params.KS2 + S2)
            rX2 = params.alpha2 * X2
            rg2 = rs2 - (1 - params.Y2) * rX2

            push!(rg1_data[i], rg1)
            push!(rg2_data[i], rg2)
        end
    end

    # Calculate differences
    rg1_diff = rg1_data[1] .- rg1_data[2]
    rg2_diff = rg2_data[1] .- rg2_data[2]
    rg_diff = (rg1_data[1] .+ rg2_data[1]) .- (rg1_data[2] .+ rg2_data[2])

    # Symmetric Y-axis range centered on 0
    diff_min = min(minimum(rg_diff), minimum(rg1_diff), minimum(rg2_diff))
    diff_max = max(maximum(rg_diff), maximum(rg1_diff), maximum(rg2_diff))
    diff_abs_max = max(abs(diff_min), abs(diff_max))
    diff_ylim = (-diff_abs_max * 1.1, diff_abs_max * 1.1)

    # ========== Step 2: Create plot (3×1 layout) ==========
    p = plot(layout=(3, 1), size=(1772, 1309),
             left_margin=15Plots.mm, right_margin=10Plots.mm,
             bottom_margin=12Plots.mm,
             )

    # Subplot 1: Δrg (Total specific growth rate difference)
    # Layer order: 1.T1/T2 background → 2.white base → 3.strategy color → 4.curve
    ylims!(p[1], diff_ylim)
    add_phase_background!(p[1], 0.0, params)
    add_difference_background!(p[1], t_rel, rg_diff, [:white, :white], alpha=1.0)
    add_difference_background!(p[1], t_rel, rg_diff, strategy_colors, alpha=1.0)
    plot!(p[1], t_rel, rg_diff,
          color=:black, linewidth=linewidth, label="",
          title="Δrg (Total Specific Growth Rate)" * gap_str, titlefontsize=titlefontsize,
          xlabel="", ylabel="Δ Specific Growth Rate", ylabelfontsize=ylabelfontsize,
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=diff_ylim, grid=false, legendfontsize=legendfontsize, legend=:topright)
    plot!(p[1], [NaN], [NaN], color=strategy_colors[1], linewidth=linewidth, label=strategy_names[1])
    plot!(p[1], [NaN], [NaN], color=strategy_colors[2], linewidth=linewidth, label=strategy_names[2])


    # Subplot 2: Δrg1 (S1 pathway difference)
    ylims!(p[2], diff_ylim)
    add_phase_background!(p[2], 0.0, params)
    add_difference_background!(p[2], t_rel, rg1_diff, [:white, :white], alpha=1.0)
    add_difference_background!(p[2], t_rel, rg1_diff, strategy_colors, alpha=1.0)
    plot!(p[2], t_rel, rg1_diff,
          color=:black, linewidth=linewidth, label=diff_label,
          title="Δrg₁ (S₁ Pathway)" * gap_str, titlefontsize=titlefontsize,
          xlabel="", ylabel="Δ Specific Growth Rate", ylabelfontsize=ylabelfontsize,
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
          ylims=diff_ylim, grid=false, legend=false)


    # Subplot 3: Δrg2 (S2 pathway difference)
     ylims!(p[3], diff_ylim)
    #ylims!(p[3], (-0.006, 0.006)) # for acclerating mRNA
    add_phase_background!(p[3], 0.0, params)
    add_difference_background!(p[3], t_rel, rg2_diff, [:white, :white], alpha=1.0)
    add_difference_background!(p[3], t_rel, rg2_diff, strategy_colors, alpha=1.0)
    plot!(p[3], t_rel, rg2_diff,
          color=:black, linewidth=linewidth, label=diff_label,
          title="Δrg₂ (S₂ Pathway)" * gap_str, titlefontsize=titlefontsize,
          xlabel="Time", xlabelfontsize=xlabelfontsize,
          ylabel="Δ Specific Growth Rate", ylabelfontsize=ylabelfontsize,
          xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
           ylims=diff_ylim, grid=false, legend=false)
          #ylims=(-0.006, 0.006), grid=false, legend=false)

    # Add subplot labels to top-left corner (outside axis area)
    add_subplot_label!(p[1], "(a)", fontsize=titlefontsize, y_offset=1.1, x_offset=-0.06)
    add_subplot_label!(p[2], "(b)", fontsize=titlefontsize, y_offset=1.1, x_offset=-0.06)
    add_subplot_label!(p[3], "(c)", fontsize=titlefontsize, y_offset=1.1, x_offset=-0.06)

    println("✅ P5 created!")
    println("   Difference Y-lim: $diff_ylim")
    return p
end


# ============================================================================
# Intelligent file naming function
# ============================================================================
function generate_plot_filename(params, strategy_indices, plot_type;
                                format="pdf", output_dir="plot_results")
    """
    Generate intelligently named plot filename

    Parameters:
    - params: Parameters object containing D, T1, T2
    - strategy_indices: Model index array
    - plot_type: Plot type name (e.g., "P1_biomass", "P2_variables")
    - format: File format (default "pdf")
    - output_dir: Output directory (default "plot_results")

    Returns:
    - Full file path
    """
    # Extract strategy names
    all_strategy_names = ["N", "T", "P", "B"]
    strategy_names_str = join([all_strategy_names[i+1] for i in strategy_indices], "-")

    # Extract key parameters
    D = params.D
    T1 = params.T1
    T2 = params.T2

    # Generate filename
    filename = "$(plot_type)_strategies-$(strategy_names_str)_D-$(D)_T1-$(T1)_T2-$(T2).$(format)"

    # Ensure output directory exists
    if !isdir(output_dir)
        mkpath(output_dir)
        println("📁 Created output directory: $(output_dir)")
    end

    full_path = joinpath(output_dir, filename)

    return full_path
end


# ============================================================================
# Save all plots as PDF files
# ============================================================================
function save_all_plots(p1, p2, p3, p4, p5, params, strategy_indices;
                        format="pdf",
                        output_dir="plot_results",
                        dpi=300)
    """
    Save all P1-P5 plots as PDF files

    Parameters:
    - p1, p2, p3, p4, p5: Plot objects
    - params: Parameters object
    - strategy_indices: Model index array
    - format: File format (default "pdf")
    - output_dir: Output directory
    - dpi: Resolution (only for PNG)
    """

    println("\n" * "="^70)
    println("🎨 Saving all competition analysis plots...")
    println("="^70)

    saved_files = String[]

    # P1: Biomass Dynamics
    println("\n📊 [1/5] Saving P1: Biomass Dynamics...")
    try
        filename1 = generate_plot_filename(params, strategy_indices, "P1_biomass_dynamics",
                                          format=format, output_dir=output_dir)
        savefig(p1, filename1)
        push!(saved_files, filename1)
        println("✅ Saved: $(basename(filename1))")
    catch e
        println("❌ Failed to save P1: $e")
    end

    # P2: Variables
    println("\n📊 [2/5] Saving P2: Variables...")
    try
        filename2 = generate_plot_filename(params, strategy_indices, "P2_variables",
                                          format=format, output_dir=output_dir)
        savefig(p2, filename2)
        push!(saved_files, filename2)
        println("✅ Saved: $(basename(filename2))")
    catch e
        println("❌ Failed to save P2: $e")
    end

    # P3: Variable Differences
    if p3 !== nothing
        println("\n📊 [3/5] Saving P3: Variable Differences...")
        try
            filename3 = generate_plot_filename(params, strategy_indices, "P3_variable_differences",
                                              format=format, output_dir=output_dir)
            savefig(p3, filename3)
            push!(saved_files, filename3)
            println("✅ Saved: $(basename(filename3))")
        catch e
            println("❌ Failed to save P3: $e")
        end
    else
        println("\n⏭️ [3/5] Skipping P3 (requires 2 strategies)")
    end

    # P4: Specific Growth Rates
    println("\n📊 [4/5] Saving P4: Specific Growth Rates...")
    try
        filename4 = generate_plot_filename(params, strategy_indices, "P4_growth_rates",
                                          format=format, output_dir=output_dir)
        savefig(p4, filename4)
        push!(saved_files, filename4)
        println("✅ Saved: $(basename(filename4))")
    catch e
        println("❌ Failed to save P4: $e")
    end

    # P5: Specific Growth Rate Differences
    if p5 !== nothing
        println("\n📊 [5/5] Saving P5: Specific Growth Rate Differences...")
        try
            filename5 = generate_plot_filename(params, strategy_indices, "P5_growth_rate_differences",
                                              format=format, output_dir=output_dir)
            savefig(p5, filename5)
            push!(saved_files, filename5)
            println("✅ Saved: $(basename(filename5))")
        catch e
            println("❌ Failed to save P5: $e")
        end
    else
        println("\n⏭️ [5/5] Skipping P5 (requires 2 strategies)")
    end

    # Print summary
    println("\n" * "="^70)
    println("🎉 All plots saved successfully!")
    println("="^70)
    println("📁 Output directory: $(output_dir)")
    println("📄 File format: $(uppercase(format))")
    println("📊 Number of files saved: $(length(saved_files))")
    println("\nSaved files:")
    for (i, file) in enumerate(saved_files)
        println("  $i. $(basename(file))")
    end
    println("="^70)

    return saved_files
end


# ============================================================================
# Main function with generation-based simulation time and IC loading
# ============================================================================
function main()
    println("\n" * "="^80)
    println("CHEMOSTAT COMPETITION - QUADGK INTEGRATION")
    println("="^80)
    
    # Configuration
    strategy_indices = [0,1]  # Strategies to compare
    
    # Set target generations
    target_tau = 1e4 # 10k generations

    # Fallback ICs (used if file loading fails)
    fallback_initial_conditions = [
        [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.88],
        [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.88]
    ]
    fallback_initial_biomass = [0.01, 0.01]
    
    # ==================================================
    # Tau conversion and IC loading
    # ==================================================
    
    # 1. Create params to get D
    params = Parameters()
    
    # 2. Calculate sim_time from tau
    sim_time = (target_tau * log(2)) / params.D
    println("Target Tau: $(target_tau) generations")
    println("Dilution Rate (D): $(params.D)")
    println("Calculated Simulation Time (t_target): $(sim_time)")

    # 3. Attempt to load ICs from file
    local u0 = nothing
    if length(strategy_indices) == 2
        strategy_pair_str = "$(strategy_indices[1])-$(strategy_indices[2])"
        
        
        filename = "independent_grid_combined_summary_D_$(params.D)_$(strategy_pair_str).csv"
        
        filepath = joinpath(IC_FILE_PATH_PREFIX, filename)
        
        println("\nAttempting to load initial conditions from: $filepath")
        
        try
            u0 = calculate_initial_conditions_from_file(params, strategy_indices, filepath, 0.5)
        catch e
            println("Warning: Error reading initial conditions file. $e")
            u0 = nothing
        end
    end

    # 4. Use fallback ICs if file loading failed
    if u0 === nothing
        println("File load failed or skipped. Using provided/default fallback initial conditions.")
        u0 = [2.0, 2.0] # S1, S2
        for i in 1:length(strategy_indices)
            append!(u0, vcat(fallback_initial_conditions[i], fallback_initial_biomass[i]))
        end
    else
        println("Successfully loaded and mixed initial conditions from file.")
    end
    
    println("Final Initial Conditions: $u0")
    
    println("\nConfiguration:")
    println("  Strategies: $strategy_indices")
    println("  Sim time: $sim_time")
    println("  Period (T1+T2): $(params.T1 + params.T2)")
    println("  ON time T1: $(params.T1)")
    println("  Amplitude ε: $(params.epsilon)")
    
    # ==================================================
    # Run ODE simulation
    # ==================================================
    results, sol, sol_longterm = compute_last_cycle_integrals_quadgk(
        strategy_indices,
        params,
        u0,
        sim_time
    )
    
    # Display numerical results
    print_results(results)

    # Print ΔE1 at the first time point of the last cycle (strategy 1 E1 - strategy 2 E1)
    if length(strategy_indices) == 2
        delta_E1 = sol[1][5] - sol[1][13]
        println("\nΔE1 at first time point (t=$(sol.t[1])): ", delta_E1)
    end
    
    # ==================================================
    # Generate P1-P5 plots
    # ==================================================
    println("\n" * "="^80)
    println("GENERATING PLOTS")
    println("="^80)

    # P1: Biomass Dynamics (long-term competition)
    println("\n⏳ Generating P1: Biomass Dynamics During Competition...")
    p1 = plot_p1_biomass(sol_longterm, params, strategy_indices)
    display(p1)

    # P2: Variables in a Single Period at Steady State
    println("\n⏳ Generating P2: Variables in a Single Period...")
    p2 = plot_p2_variables(sol, params, strategy_indices)
    display(p2)

    # P3: Variable Differences (requires 2 strategies)
    println("\n⏳ Generating P3: Variable Differences...")
    p3 = plot_p3_variable_differences(sol, params, strategy_indices)
    if p3 !== nothing
        display(p3)
    end

    # P4: Specific Growth Rates
    println("\n⏳ Generating P4: Specific Growth Rates...")
    p4, rg_data, rg1_data, rg2_data, t_eval = plot_p4_growth_rates(sol, params, strategy_indices)
    display(p4)

    # P5: Specific Growth Rate Differences (requires 2 strategies)
    println("\n⏳ Generating P5: Specific Growth Rate Differences...")
    p5 = plot_p5_growth_rate_differences(sol, params, strategy_indices)
    if p5 !== nothing
        display(p5)
    end

    println("\n" * "="^80)
    println("✅ All plots generated successfully!")
    println("="^80)

    # ==================================================
    # Save all plots to PDF
    # ==================================================
    println("\n⏳ Saving plots to PDF...")
    saved_files = save_all_plots(p1, p2, p3, p4, p5, params, strategy_indices;
                                 format="pdf", output_dir="plot_results")

    return results, sol, sol_longterm, p1, p2, p3, p4, p5
end


# Run main function
results, sol, sol_longterm, p1, p2, p3, p4, p5 = main();