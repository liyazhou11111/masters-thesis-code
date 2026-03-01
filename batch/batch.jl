# ==============================================================================
# Batch Culture Simulation
#
# Independent batch culture simulations for three CCR strategies:
#   Strategy N (strategy_idx=0): No CCR
#   Strategy T (strategy_idx=1): Transcriptional CCR
#   Strategy P (strategy_idx=2): Post-transcriptional CCR
#
# State layout per strategy:
#   u = [S1, X1, M1, E1, S2, X2, M2, E2, c_minus, c]
#
# Note: Batch culture has no dilution (D=0), so dc/dt = rg * c
#       and dS/dt = -uptake * c (substrate depletion only).
# ==============================================================================

using DifferentialEquations
using Plots
using Printf
using Plots.PlotMeasures


# ==============================================================================
# Parameter Struct
# ==============================================================================
mutable struct Parameters
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
        0.0,    # delta2
        0.01,   # KX1
        0.01,   # KX2
        0.001,  # b
        
        # CCR parameters
        0.01,   # K_CCR1
        0.01,   # K_CCR2
        
        # Yield coefficients
        0.7,    # Y1
        0.7     # Y2
    )
end


# ==============================================================================
# Strategy N: No CCR (Independent batch culture)
# ==============================================================================
function strategy_N_batch!(du, u, p, t)
    S1, X1, M1, E1, S2, X2, M2, E2, c_minus, c = u
    params = p

    # Growth rate (Eq. 20)
    rg = (params.k1 * E1 * S1 / (params.KS1 + S1) +
          params.k2 * E2 * S2 / (params.KS2 + S2)) -
         ((1 - params.Y1) * params.alpha1 * X1 +
          (1 - params.Y2) * params.alpha2 * X2)

    # CCR regulation: Strategy N has no CCR
    kappa_M = 1.0
    kappa_E = 1.0

    # Differential equations (batch: no dilution D, substrate consumed by cells)
    # dS1/dt
    du[1] = -params.k1 * E1 * S1 / (params.KS1 + S1) * c
    # dX1/dt (Eq. 8)
    du[2] = params.k1 * E1 * S1 / (params.KS1 + S1) -
            params.alpha1 * X1 - rg * X1
    # dM1/dt (Eq. 10)
    du[3] = params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * c_minus -
            params.delta1 * M1 - rg * M1
    # dE1/dt (Eq. 14)
    du[4] = params.kE * M1 * c_minus -
            params.delta2 * E1 - rg * E1
    # dS2/dt
    du[5] = -params.k2 * E2 * S2 / (params.KS2 + S2) * c
    # dX2/dt (Eq. 9)
    du[6] = params.k2 * E2 * S2 / (params.KS2 + S2) -
            params.alpha2 * X2 - rg * X2
    # dM2/dt (Eq. 11) — kappa_M = 1.0 for Strategy N
    du[7] = params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * c_minus * kappa_M -
            params.delta1 * M2 - rg * M2
    # dE2/dt (Eq. 15) — kappa_E = 1.0 for Strategy N
    du[8] = params.kE * M2 * c_minus * kappa_E -
            params.delta2 * E2 - rg * E2
    # dc_minus/dt (Eq. 18)
    du[9] = (params.Y1 * params.alpha1 * X1 + params.Y2 * params.alpha2 * X2) +
            (params.delta1 * M1 + params.delta2 * E1 +
             params.delta1 * M2 + params.delta2 * E2) -
            (params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * c_minus +
             params.kE * M1 * c_minus +
             params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * c_minus * kappa_M +
             params.kE * M2 * c_minus * kappa_E) -
            rg * c_minus
    # dc/dt (batch: no dilution)
    du[10] = rg * c

    return nothing
end


# ==============================================================================
# Strategy T: Transcriptional CCR (Independent batch culture)
# ==============================================================================
function strategy_T_batch!(du, u, p, t)
    S1, X1, M1, E1, S2, X2, M2, E2, c_minus, c = u
    params = p

    # Growth rate (Eq. 20)
    rg = (params.k1 * E1 * S1 / (params.KS1 + S1) +
          params.k2 * E2 * S2 / (params.KS2 + S2)) -
         ((1 - params.Y1) * params.alpha1 * X1 +
          (1 - params.Y2) * params.alpha2 * X2)

    # CCR regulation: Strategy T — transcriptional repression of M2 (Eq. 13)
    kappa_M = params.K_CCR1 / (params.K_CCR1 + X1)
    kappa_E = 1.0

    # Differential equations (batch: no dilution D, substrate consumed by cells)
    # dS1/dt
    du[1] = -params.k1 * E1 * S1 / (params.KS1 + S1) * c
    # dX1/dt (Eq. 8)
    du[2] = params.k1 * E1 * S1 / (params.KS1 + S1) -
            params.alpha1 * X1 - rg * X1
    # dM1/dt (Eq. 10)
    du[3] = params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * c_minus -
            params.delta1 * M1 - rg * M1
    # dE1/dt (Eq. 14)
    du[4] = params.kE * M1 * c_minus -
            params.delta2 * E1 - rg * E1
    # dS2/dt
    du[5] = -params.k2 * E2 * S2 / (params.KS2 + S2) * c
    # dX2/dt (Eq. 9)
    du[6] = params.k2 * E2 * S2 / (params.KS2 + S2) -
            params.alpha2 * X2 - rg * X2
    # dM2/dt (Eq. 11) — includes kappa_M for transcriptional CCR
    du[7] = params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * c_minus * kappa_M -
            params.delta1 * M2 - rg * M2
    # dE2/dt (Eq. 15) — kappa_E = 1.0 for Strategy T
    du[8] = params.kE * M2 * c_minus * kappa_E -
            params.delta2 * E2 - rg * E2
    # dc_minus/dt (Eq. 18)
    du[9] = (params.Y1 * params.alpha1 * X1 + params.Y2 * params.alpha2 * X2) +
            (params.delta1 * M1 + params.delta2 * E1 +
             params.delta1 * M2 + params.delta2 * E2) -
            (params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * c_minus +
             params.kE * M1 * c_minus +
             params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * c_minus * kappa_M +
             params.kE * M2 * c_minus * kappa_E) -
            rg * c_minus
    # dc/dt (batch: no dilution)
    du[10] = rg * c

    return nothing
end


# ==============================================================================
# Strategy P: Post-transcriptional CCR (Independent batch culture)
# ==============================================================================
function strategy_P_batch!(du, u, p, t)
    S1, X1, M1, E1, S2, X2, M2, E2, c_minus, c = u
    params = p

    # Growth rate (Eq. 20)
    rg = (params.k1 * E1 * S1 / (params.KS1 + S1) +
          params.k2 * E2 * S2 / (params.KS2 + S2)) -
         ((1 - params.Y1) * params.alpha1 * X1 +
          (1 - params.Y2) * params.alpha2 * X2)

    # CCR regulation: Strategy P — post-transcriptional repression of E2 (Eq. 17)
    kappa_M = 1.0
    kappa_E = params.K_CCR2 / (params.K_CCR2 + X1)

    # Differential equations (batch: no dilution D, substrate consumed by cells)
    # dS1/dt
    du[1] = -params.k1 * E1 * S1 / (params.KS1 + S1) * c
    # dX1/dt (Eq. 8)
    du[2] = params.k1 * E1 * S1 / (params.KS1 + S1) -
            params.alpha1 * X1 - rg * X1
    # dM1/dt (Eq. 10)
    du[3] = params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * c_minus -
            params.delta1 * M1 - rg * M1
    # dE1/dt (Eq. 14)
    du[4] = params.kE * M1 * c_minus -
            params.delta2 * E1 - rg * E1
    # dS2/dt
    du[5] = -params.k2 * E2 * S2 / (params.KS2 + S2) * c
    # dX2/dt (Eq. 9)
    du[6] = params.k2 * E2 * S2 / (params.KS2 + S2) -
            params.alpha2 * X2 - rg * X2
    # dM2/dt (Eq. 11) — kappa_M = 1.0 for Strategy P
    du[7] = params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * c_minus * kappa_M -
            params.delta1 * M2 - rg * M2
    # dE2/dt (Eq. 15) — includes kappa_E for post-transcriptional CCR
    du[8] = params.kE * M2 * c_minus * kappa_E -
            params.delta2 * E2 - rg * E2
    # dc_minus/dt (Eq. 18)
    du[9] = (params.Y1 * params.alpha1 * X1 + params.Y2 * params.alpha2 * X2) +
            (params.delta1 * M1 + params.delta2 * E1 +
             params.delta1 * M2 + params.delta2 * E2) -
            (params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * c_minus +
             params.kE * M1 * c_minus +
             params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * c_minus * kappa_M +
             params.kE * M2 * c_minus * kappa_E) -
            rg * c_minus
    # dc/dt (batch: no dilution)
    du[10] = rg * c

    return nothing
end


# ==============================================================================
# Initial Conditions
# ==============================================================================
function get_initial_conditions()
    return [
        2.0,    # S1
        0.02,   # X1
        0.02,   # M1
        0.02,   # E1
        2.0,    # S2
        0.02,   # X2
        0.02,   # M2
        0.02,   # E2
        0.88,   # c_minus
        0.01    # c
    ]
end


# ==============================================================================
# Run Individual Batch Simulation
# ==============================================================================
function run_single_batch_simulation(strategy_func, strategy_name; sim_time=50.0, num_points=500)
    params = Parameters()
    initial_conditions = get_initial_conditions()
    
    t_span = (0.0, sim_time)
    t_eval = range(t_span[1], t_span[2], length=num_points)
    
    println("Running $(strategy_name) batch culture simulation...")
    
    prob = ODEProblem(strategy_func, initial_conditions, t_span, params)
    
    try
        sol = solve(prob, Tsit5(), 
                   saveat = t_eval,
                   save_everystep = false,
                   dense = false,
                   abstol = 1e-10, reltol = 1e-8,
                   maxiters = 1e12)
        
        println("✅ $(strategy_name) simulation completed successfully!")
        return sol, collect(t_eval), params
    catch e
        println("❌ $(strategy_name) simulation failed: ", e)
        return nothing, nothing, nothing
    end
end



# ==============================================================================
# Main Function: Run All Three Independent Strategies
# ==============================================================================
function run_all_batch_models()
    println("Running three independent batch culture strategies...")
    
    # Strategy N: No CCR
    println("\n" * "="^50)
    println("STRATEGY N: No CCR")
    println("="^50)
    sol1, t_eval1, params1 = run_single_batch_simulation(strategy_N_batch!, "No CCR (N) Strategy")
    
    # Strategy T: Transcriptional CCR
    println("\n" * "="^50)
    println("STRATEGY T: Transcriptional CCR")
    println("="^50)
    sol2, t_eval2, params2 = run_single_batch_simulation(strategy_T_batch!, "Transcriptional CCR (T) Strategy")
    
    # Strategy P: Post-transcriptional CCR
    println("\n" * "="^50)
    println("STRATEGY P: Post-transcriptional CCR")
    println("="^50)
    sol3, t_eval3, params3 = run_single_batch_simulation(strategy_P_batch!, "Post-transcriptional CCR (P) Strategy")
    
    println("\n" * "="^50)
    println("All batch culture simulations completed!")
    println("="^50)
    
    results = ((sol1, t_eval1, params1), (sol2, t_eval2, params2), (sol3, t_eval3, params3))

    return results
end


# ==============================================================================
# Growth Rate Comparison Plot
# ==============================================================================
function add_subplot_label!(p, label; fontsize=28, x_offset=-0.5, y_offset=-1)
    xl = xlims(p)
    yl = ylims(p)
    x_pos = xl[1] + x_offset * (xl[2] - xl[1])
    y_pos = yl[1] + y_offset * (yl[2] - yl[1])
    annotate!(p, [(x_pos, y_pos, text(label, :left, :bottom, fontsize, :bold))])
end

function plot_growth_rate_comparison(results;
                                    linewidth=3,
                                    titlefontsize=22,
                                    xlabelfontsize=20,
                                    ylabelfontsize=20,
                                    tickfontsize=18,
                                    legendfontsize=18,
                                    label_x_offset=-0.37,
                                    label_y_offset=1,
                                    substrate_ylabel_x_offset=[-0.28, -0.33, -0.33])
    # Extract solutions
    (sol1, t_eval1, params1), (sol2, t_eval2, params2), (sol3, t_eval3, params3) = results

    if sol1 === nothing || sol2 === nothing || sol3 === nothing
        println("One or more solutions are missing, cannot create comparison plot")
        return nothing
    end

    # Strategy colors for dc/dt
    dcdt_colors = [RGB(0.95, 0.8, 0.25),   # N — gold
                   RGB(0.35, 0.8, 0.9),     # T — cyan
                   RGB(0.95, 0.6, 0.5)]     # P — salmon

    strategy_short = ["N", "T", "P"]
    solutions = [sol1, sol2, sol3]
    t_evals = [t_eval1, t_eval2, t_eval3]

    # Compute shared x-range across all strategies
    t_min = minimum(first(te) for te in t_evals)
    t_max = maximum(last(te) for te in t_evals)

    # 2×3 layout: top row dc/dt, bottom row substrates
    lay = @layout [a b c
                   d e f]
    p = plot(layout=lay, size=(1772, 800),
             left_margin=15mm,
             right_margin=10mm)

    subplot_labels = ["(a)", "(b)", "(c)"]

    for (i, (sol, t_eval)) in enumerate(zip(solutions, t_evals))
        S1_data = sol[1, :]
        S2_data = sol[5, :]
        c_data = sol[10, :]

        dt = t_eval[2] - t_eval[1]
        dc_dt = diff(c_data) ./ dt
        t_mid = t_eval[1:end-1] .+ dt / 2

        # ---- Top row (subplot i): dc/dt ----
        plot!(p[i], t_mid, dc_dt,
              color=dcdt_colors[i], linewidth=linewidth, label=strategy_short[i],
              ylabel="Growth Rate (dc/dt)",
              ylabelfontsize=ylabelfontsize,
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              legendfontsize=legendfontsize,
              xlims=(t_min, t_max),
              xformatter=_->"",
              top_margin=10mm,
              bottom_margin=-3mm,
              grid=false,
              legend=:topright)

        # ---- Bottom row (subplot i+3): S1 and S2 ----
        plot!(p[i+3], t_eval, S1_data,
              color=:red, linewidth=linewidth, label="S₁",
              xlabel="Time", xlabelfontsize=xlabelfontsize,
              ylabel="",
              xtickfontsize=tickfontsize, ytickfontsize=tickfontsize,
              legendfontsize=legendfontsize,
              xlims=(t_min, t_max),
              top_margin=-3mm,
              bottom_margin=12mm,
              grid=false,
              legend=:topright)

        plot!(p[i+3], t_eval, S2_data,
              color=:blue, linewidth=linewidth, label="S₂")

        # Manual ylabel for precise position control
        xl_s = xlims(p[i+3]); yl_s = ylims(p[i+3])
        annotate!(p[i+3], xl_s[1] + substrate_ylabel_x_offset[i] * (xl_s[2] - xl_s[1]),
                  (yl_s[1] + yl_s[2]) / 2,
                  text("Carbon source Conc.", :center, ylabelfontsize, rotation=90))

        # Add subplot label to top row only
        add_subplot_label!(p[i], subplot_labels[i],
                          fontsize=titlefontsize,
                          x_offset=label_x_offset,
                          y_offset=label_y_offset)

        println("\n Strategy $(strategy_short[i]) Growth Rate:")
        println("  Max dc/dt: $(round(maximum(dc_dt), digits=4))")
        println("  Final dc/dt: $(round(dc_dt[end], digits=4))")
    end

    display(p)
    return p
end


# ==============================================================================
# Plot Filename Generation
# ==============================================================================
function generate_batch_plot_filename(plot_type;
                                     format="pdf",
                                     output_dir="plot_results")
    # Generate filename
    filename = "batch_$(plot_type)_strategies-N-T-P.$(format)"

    # Ensure output directory exists
    if !isdir(output_dir)
        mkpath(output_dir)
        println("Created output directory: $(output_dir)")
    end

    return joinpath(output_dir, filename)
end


# ==============================================================================
# Save Batch Plots
# ==============================================================================
function save_batch_plots(results;
                          format="pdf",
                          output_dir="plot_results",
                          dpi=300,
                          kwargs...)
    println("\n" * "="^70)
    println("Saving batch culture plots...")
    println("="^70)

    saved_files = String[]

    # Growth rate comparison (dc/dt + substrates)
    println("\n Saving: Growth rate comparison...")
    try
        p2 = plot_growth_rate_comparison(results; kwargs...)
        plot!(p2, dpi=dpi)
        filename2 = generate_batch_plot_filename("growth_rate",
                    format=format, output_dir=output_dir)
        savefig(p2, filename2)
        push!(saved_files, filename2)
        println("  Saved: $(basename(filename2))")
    catch e
        println("  Failed: $e")
    end

    # Summary
    println("\n" * "="^70)
    println("All batch plots saved!")
    println("="^70)
    println("Output directory: $(output_dir)")
    println("File format: $(uppercase(format))")
    println("Files saved: $(length(saved_files))")
    for (i, file) in enumerate(saved_files)
        println("  $i. $(basename(file))")
    end
    println("="^70)

    return saved_files
end


# ==============================================================================
# Run all strategies
# ==============================================================================
results = run_all_batch_models();
save_batch_plots(results)