#!/usr/bin/env julia

# ================================================================================
# FILE: run_competition_chemostat_grid.jl
#
# Pairwise competition parameter sweep over (T1, T2) environmental conditions.
# For each (T1, T2), runs competition across a grid of initial volume fractions
# to detect multi-stability and classify outcomes.
#
# Supports strategy pairs: 0-1 (N vs T), 1-2 (T vs P), 0-2 (N vs P)
#
# Usage:
#   julia run_competition_chemostat_grid.jl <T1> <T2>
#
# Environment variables:
#   CCR_D_VALUE              Dilution rate (default: 0.04)
#   CCR_STRATEGY_PAIR        Strategy pair: "0-1", "1-2", "0-2" (default: "1-2")
#   CCR_GRID_SIZE            Number of volume fraction grid points (default: 9)
#   CCR_VOLUME_MIN           Minimum volume fraction (default: 1e-5)
#   CCR_TAU                  Target generations (default: 1e7)
#   CCR_ANALYSIS_CYCLES      Cycles for analysis window (default: 50)
#   CCR_EXTINCTION_THRESHOLD Biomass extinction threshold (default: 1e-100)
#   CCR_NE                   Effective population size (default: 1e8)
#   CCR_REGRESSION_CYCLES    Cycles for regression (default: 50)
#   CCR_INITIAL_CONDITIONS_FILE  Path to equilibrium CSV
#   CCR_USE_FILE_CONDITIONS  Use file-based initial conditions (default: false)
# ================================================================================

using CSV
using DataFrames
using Dates
using Statistics
using GLM
using Base.Threads

println("=== CCR Pairwise Competition Parameter Sweep ===")
println("Julia version: $(VERSION)")
println("Start time: $(now())")

include("ccr_competition_chemostat_definitions.jl")

function main()
    if length(ARGS) < 2
        println("Usage: julia run_competition_chemostat_grid.jl T1_VALUE T2_VALUE")
        exit(1)
    end
    
    T1 = parse(Float64, ARGS[1])
    T2 = parse(Float64, ARGS[2])
    
    # Get parameters from environment
    D = haskey(ENV, "CCR_D_VALUE") ? parse(Float64, ENV["CCR_D_VALUE"]) : 0.04
    
    # Volume fraction grid parameters
    grid_size = haskey(ENV, "CCR_GRID_SIZE") ? parse(Int, ENV["CCR_GRID_SIZE"]) : 9
    volume_min = haskey(ENV, "CCR_VOLUME_MIN") ? parse(Float64, ENV["CCR_VOLUME_MIN"]) : 1e-5
    
    target_tau = haskey(ENV, "CCR_TAU") ? parse(Float64, ENV["CCR_TAU"]) : 1e7
    analysis_cycles = haskey(ENV, "CCR_ANALYSIS_CYCLES") ? parse(Int, ENV["CCR_ANALYSIS_CYCLES"]) : 50
    extinction_threshold = haskey(ENV, "CCR_EXTINCTION_THRESHOLD") ? parse(Float64, ENV["CCR_EXTINCTION_THRESHOLD"]) : 1e-100
    Ne = haskey(ENV, "CCR_NE") ? parse(Float64, ENV["CCR_NE"]) : 1e8
    n_cycles_for_regression = haskey(ENV, "CCR_REGRESSION_CYCLES") ? parse(Int, ENV["CCR_REGRESSION_CYCLES"]) : 50
    
    # Strategy pair configuration
    strategy_pair = haskey(ENV, "CCR_STRATEGY_PAIR") ? strip(replace(ENV["CCR_STRATEGY_PAIR"], r"[\'\"]" => "")) : "1-2"
    
    if strategy_pair == "0-1"
        strategy_indices = [0, 1]
    elseif strategy_pair == "1-2"
        strategy_indices = [1, 2]
    elseif strategy_pair == "0-2"
        strategy_indices = [0, 2]
    else
        error("Invalid CCR_STRATEGY_PAIR: $(strategy_pair). Must be '0-1', '1-2', or '0-2'")
    end
    
    # Initial conditions file support
    initial_conditions_file = haskey(ENV, "CCR_INITIAL_CONDITIONS_FILE") ? ENV["CCR_INITIAL_CONDITIONS_FILE"] : nothing
    use_file_initial_conditions = haskey(ENV, "CCR_USE_FILE_CONDITIONS") ? parse(Bool, ENV["CCR_USE_FILE_CONDITIONS"]) : false
    
    # Output directory
    base_output = haskey(ENV, "CCR_OUTPUT_DIR") ? ENV["CCR_OUTPUT_DIR"] : "results"
    output_dir = joinpath(base_output, "competition_pair_$(strategy_pair)_D_$(D)_grid$(grid_size)")

    mkpath(output_dir)
    
    # Create base parameters
    base_params = Parameters()
    base_params.D = D
    
    println("\n=== Processing Single T1,T2 Combination ===")
    println("T1: $(T1), T2: $(T2)")
    println("Strategy pair: $(strategy_pair) (Strategy $(strategy_indices[1]) vs Strategy $(strategy_indices[2]))")
    println("Grid size: $(grid_size) volume fractions")
    println("Volume min: $(volume_min), max: $(1-volume_min)")
    println("Initial conditions: $(use_file_initial_conditions ? "From file" : "Default")")
    
    # Process this single T1,T2 combination
    process_parameter_combination(T1, T2, base_params, strategy_indices, strategy_pair,
                                grid_size, volume_min,
                                target_tau, analysis_cycles, extinction_threshold,
                                Ne, n_cycles_for_regression, output_dir,
                                initial_conditions_file, use_file_initial_conditions)
    
    println("\n=== Completed: $(now()) ===")
end


function process_parameter_combination(T1, T2, base_params, strategy_indices, strategy_pair,
                                     grid_size, volume_min,
                                     target_tau, analysis_cycles, extinction_threshold,
                                     Ne, n_cycles_for_regression, output_dir,
                                     initial_conditions_file, use_file_initial_conditions)
    """Process a single (T1, T2) combination with volume fraction grid."""
    
    # Setup parameters for this combination
    params = deepcopy(base_params)
    params.T1 = T1
    params.T2 = T2
    
    # Generate volume fraction grid for this parameter combination
    volume_fractions = generate_volume_fraction_grid(grid_size, volume_min)
    
    # Pre-compute all initial conditions for parallel processing
    println("Pre-computing initial conditions for $(length(volume_fractions)) volume fractions...")
    
    all_initial_conditions = Vector{Union{Vector{Float64}, Nothing}}(undef, length(volume_fractions))
    
    for (idx, Z) in enumerate(volume_fractions)
        if use_file_initial_conditions && initial_conditions_file !== nothing
            # Calculate from file
            u0 = calculate_initial_conditions_from_file(params, strategy_indices, 
                                                       initial_conditions_file, Z)
            if u0 === nothing
                println("Warning: Failed to generate initial conditions from file for Z=$(Z), using defaults")
                u0 = generate_default_initial_conditions(strategy_indices, Z)
            end
        else
            # Use defaults
            u0 = generate_default_initial_conditions(strategy_indices, Z)
        end
        
        all_initial_conditions[idx] = u0
    end
    
    # Filter out failed initial conditions
    valid_indices = findall(u0 -> u0 !== nothing, all_initial_conditions)
    if length(valid_indices) < length(volume_fractions)
        println("Warning: $(length(volume_fractions) - length(valid_indices)) initial conditions failed, proceeding with $(length(valid_indices)) valid ones")
    end
    
    # Create grid combinations for valid initial conditions
    all_grid_combinations = []
    for idx in valid_indices
        Z = volume_fractions[idx]
        u0 = all_initial_conditions[idx]
        push!(all_grid_combinations, (
            grid_id = idx,
            volume_fraction = Z,
            initial_conditions = u0,
            strategy_pair = strategy_pair
        ))
    end
    
    total_combinations = length(all_grid_combinations)
    
    # Parallel processing of grid points for this T1,T2 combination
    results_array = Vector{Any}(undef, total_combinations)
    start_time_combo = time()
    
    @threads for combo_idx in 1:total_combinations
        grid_point = all_grid_combinations[combo_idx]
        Z = grid_point.volume_fraction
        u0 = grid_point.initial_conditions
        
        result_dict_to_store = Dict()
        simulation_duration = 0.0
        analysis_duration = 0.0
        
        start_time_point = time()
        sol = nothing
        try
            # Simple call with pre-computed initial conditions
            sol = run_grid_simulation(params, strategy_indices, u0; 
                                      Ne=Ne, target_tau=target_tau, analysis_cycles=analysis_cycles)
            simulation_duration = time() - start_time_point
            
            if sol !== nothing
                analysis_start = time()
                result = analyze_grid_simulation_outcome(sol, params, strategy_indices, Ne;
                                                         extinction_threshold=extinction_threshold,
                                                         n_cycles_for_regression=n_cycles_for_regression,
                                                         verbose=false)
                analysis_duration = time() - analysis_start
                
                if result !== nothing
                    result_dict_to_store = result
                else
                    result_dict_to_store = Dict("outcome" => "Analysis Returned Null", "winner" => -998)
                end
            else
                result_dict_to_store = Dict("outcome" => "Simulation Failed", "winner" => -999)
            end
        catch e
            result_dict_to_store = Dict("outcome" => "Simulation Exception", "winner" => -997)
            println("Exception for Z=$(Z): $(e)")
        end

        # Store result with volume fraction information
        stored_result = (
            grid_id=grid_point.grid_id,
            volume_fraction=Z,
            strategy_pair=strategy_pair,
            result_dict=result_dict_to_store,
            simulation_duration=simulation_duration,
            analysis_duration=analysis_duration,
            total_duration=time() - start_time_point
        )
        
        results_array[combo_idx] = stored_result
    end
    
    combo_elapsed_time = time() - start_time_combo
    println("  ✅ T1=$(T1), T2=$(T2) completed in $(round(combo_elapsed_time/60, digits=2)) minutes")
    
    # Save results for this T1,T2 combination
    save_parameter_combination_results(results_array, T1, T2, params, strategy_indices, strategy_pair,
                                     grid_size, volume_min, target_tau,
                                     analysis_cycles, extinction_threshold, Ne, n_cycles_for_regression,
                                     use_file_initial_conditions, output_dir)
end



function save_parameter_combination_results(detailed_results, T1, T2, params, strategy_indices, strategy_pair,
                                          grid_size, volume_min, target_tau,
                                          analysis_cycles, extinction_threshold, Ne, n_cycles_for_regression,
                                          use_file_initial_conditions, output_dir)
    """Save detailed and summary results for a specific (T1, T2) combination."""

    if isempty(detailed_results)
        println("  Warning: No results to save for T1=$(T1), T2=$(T2)")
        return
    end

    total_combinations = length(detailed_results)
    target_physical_time = (target_tau * log(2)) / params.D
    cycle_length = T1 + T2
    volume_max = 1.0 - volume_min

    # ======================================================================
    # Detail CSV
    # ======================================================================
    results_df = DataFrame(
        grid_id = [r.grid_id for r in detailed_results],
        volume_fraction = [r.volume_fraction for r in detailed_results],
        volume_fraction_complement = [1.0 - r.volume_fraction for r in detailed_results],
        strategy_pair = [r.strategy_pair for r in detailed_results],
    )

    # Fixed columns
    results_df[!, :T1] .= T1
    results_df[!, :T2] .= T2
    results_df[!, :D] .= params.D

    # Analysis output
    results_df[!, :outcome] = [get(r.result_dict, "outcome", "Unknown") for r in detailed_results]
    results_df[!, :refined_outcome] = [get(r.result_dict, "refined_outcome",
                                           get(r.result_dict, "outcome", "Unknown")) for r in detailed_results]
    results_df[!, :winner] = [get(r.result_dict, "winner", -999) for r in detailed_results]
    results_df[!, :survivors] = [join(get(r.result_dict, "survivors", []), ",") for r in detailed_results]

    # Timing
    results_df[!, :simulation_duration_seconds] = [r.simulation_duration for r in detailed_results]
    results_df[!, :analysis_duration_seconds] = [r.analysis_duration for r in detailed_results]
    results_df[!, :total_duration_seconds] = [r.total_duration for r in detailed_results]
    results_df[!, :timestamp] .= string(now())

    # Biological time context
    results_df[!, :target_tau] .= target_tau
    results_df[!, :target_physical_time] .= target_physical_time
    results_df[!, :cycle_length] .= cycle_length

    # Volume fraction parameters
    results_df[!, :grid_size] .= grid_size
    results_df[!, :volume_min] .= volume_min
    results_df[!, :volume_max] .= volume_max

    # Per-strategy results (avg biomass + Ne*s related fields)
    for strategy_idx in strategy_indices
        results_df[!, "strategy_$(strategy_idx)_avg_biomass"] = [
            let avg_biomass = get(r.result_dict, "avg_biomass", [])
                if !isempty(avg_biomass)
                    array_idx = findfirst(x -> x == strategy_idx, strategy_indices)
                    array_idx !== nothing && array_idx <= length(avg_biomass) ? avg_biomass[array_idx] : NaN
                else
                    NaN
                end
            end for r in detailed_results
        ]

        nes_fields = ["survives", "survival_reason", "net_growth_rate", "r_squared", "Ne_s"]
        for field in nes_fields
            results_df[!, "strategy_$(strategy_idx)_$(field)"] = [
                let individual_results = get(r.result_dict, "individual_nes_results", Dict())
                    strategy_result = get(individual_results, strategy_idx, Dict())
                    get(strategy_result, field, field == "survives" ? false : field == "survival_reason" ? "not_analyzed" : NaN)
                end for r in detailed_results
            ]
        end
    end

    # Write detail CSV
    strategy_indices_str = join(strategy_indices, "_")


    filename = "competition_strategies_$(strategy_indices_str)_T1_$(T1)_T2_$(T2)_D_$(params.D)_grid$(grid_size)_pair_$(strategy_pair).csv"

    filepath = joinpath(output_dir, filename)
    CSV.write(filepath, results_df)

    # ======================================================================
    # Summary CSV — consensus-based classification
    # ======================================================================
    successful_results = filter(r -> !occursin("Failed", get(r.result_dict, "outcome", "")) &&
                                     !occursin("Exception", get(r.result_dict, "outcome", "")) &&
                                     !occursin("Null", get(r.result_dict, "outcome", "")),
                                detailed_results)
    num_successful = length(successful_results)

    classification, consensus = classify_consensus_outcome(successful_results, strategy_indices)

    summary_df = DataFrame(
        T1                     = [T1],
        T2                     = [T2],
        D                      = [params.D],
        strategy_pair          = [strategy_pair],
        grid_size              = [grid_size],
        volume_min             = [volume_min],
        volume_max             = [volume_max],
        total_grid_points      = [total_combinations],
        successful_simulations = [num_successful],
        success_rate_percent   = [round(num_successful / total_combinations * 100, digits=1)],
        classification         = [classification],
        consensus              = [consensus],
        target_tau             = [target_tau],
        target_physical_time   = [target_physical_time],
        cycle_length           = [cycle_length],
        timestamp              = [string(now())]
    )

    summary_filename = "competition_strategies_$(strategy_indices_str)_T1_$(T1)_T2_$(T2)_D_$(params.D)_grid$(grid_size)_pair_$(strategy_pair)_summary.csv"
    summary_filepath = joinpath(output_dir, summary_filename)
    CSV.write(summary_filepath, summary_df)

    println("    Saved: $(filename) & summary ($(classification), consensus=$(round(consensus, digits=3)))")
end


# ==============================================================================
# Classify Consensus Outcome
# ==============================================================================

"""
    classify_consensus_outcome(successful_results, strategy_indices) -> (classification, consensus)

Classify competition outcome across all initial conditions (Section 2.4.2):

  - Extinction           — all "Extinction"                          → consensus 1.0
  - Competitive exclusion — all "Exclusion", same winner             → consensus 1.0
  - Stable coexistence    — all "Coexistence", same winner code      → consensus 1.0
  - Neutral coexistence   — all "Coexistence", different winner codes → consensus < 1.0
  - Multi-stability       — mixed outcome types                      → consensus < 1.0

Classification output uses strategy indices directly:
  "1 win", "2 win", "0 win"
  "12 coexist", "21 coexist", "901 coexist", "10 coexist"
  "Neutral coexistence", "Extinction", "Multi-stability"

Strategy name mapping (0→N, 1→T, 2→P) is applied at the plotting stage.
"""
function classify_consensus_outcome(successful_results, strategy_indices)
    if isempty(successful_results)
        return ("No data", 0.0)
    end

    n = length(successful_results)

    # Extract outcome types and winner codes
    outcomes = [get(r.result_dict, "outcome", "Unknown") for r in successful_results]
    winners  = [get(r.result_dict, "winner", -999) for r in successful_results]

    # Map each outcome string to a base type
    base_types = String[]
    for o in outcomes
        if startswith(o, "Extinction")
            push!(base_types, "Extinction")
        elseif startswith(o, "Exclusion")
            push!(base_types, "Exclusion")
        elseif startswith(o, "Coexistence")
            push!(base_types, "Coexistence")
        else
            push!(base_types, "Error")
        end
    end

    unique_types = unique(base_types)

    # --- All initial conditions produce the same outcome type ---
    if length(unique_types) == 1
        the_type = unique_types[1]

        if the_type == "Extinction"
            return ("Extinction", 1.0)

        elseif the_type == "Exclusion"
            unique_winners = unique(winners)
            if length(unique_winners) == 1
                # Competitive exclusion — single winner across all ICs
                return ("$(unique_winners[1]) win", 1.0)
            else
                # All exclusion but different winners → multi-stability
                counts = Dict{Int, Int}()
                for w in winners; counts[w] = get(counts, w, 0) + 1; end
                return ("Multi-stability", maximum(values(counts)) / n)
            end

        elseif the_type == "Coexistence"
            unique_winners = unique(winners)
            if length(unique_winners) == 1
                # Stable coexistence — same biomass ranking across all ICs
                return ("$(unique_winners[1]) coexist", 1.0)
            else
                # All coexistence but different biomass rankings → neutral
                counts = Dict{Int, Int}()
                for w in winners; counts[w] = get(counts, w, 0) + 1; end
                return ("Neutral coexistence", maximum(values(counts)) / n)
            end

        else
            return ("Error", 0.0)
        end

    else
        # --- Mixed outcome types → multi-stability ---
        type_counts = Dict{String, Int}()
        for bt in base_types; type_counts[bt] = get(type_counts, bt, 0) + 1; end
        return ("Multi-stability", maximum(values(type_counts)) / n)
    end
end


# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end