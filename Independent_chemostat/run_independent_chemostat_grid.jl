#!/usr/bin/env julia

# ================================================================================
# FILE: run_independent_chemostat_grid.jl
# PURPOSE: Independent chemostat cultivation with grid-based initial biomass scanning
# 
# This script simulates independent (monoculture) growth of two CCR strategies
# across a grid of initial biomass values to assess robustness of outcomes.
# It compares average biomass achieved by each strategy under identical conditions.
#
# Strategy indices:
# - 0: Strategy N (No CCR)
# - 1: Strategy T (Transcriptional CCR)
# - 2: Strategy P (Post-transcriptional CCR)
# ================================================================================

using CSV
using DataFrames
using Dates
using Statistics
using Base.Threads
using LinearAlgebra

println("=== CCR Independent Chemostat Grid-based Robustness Test ===")
println("Julia version: $(VERSION)")
println("Start time: $(now())")

include("independent_chemostat_definitions.jl")

function main()
    if length(ARGS) < 1
        println("Usage: julia run_independent_chemostat_grid.jl T1_INDEX")
        exit(1)
    end
    
    T1_index = parse(Int, ARGS[1])
    
    # Define parameter ranges (must match job submitter)
    T1_values = collect(0:5:100)  # 41 values: 0, 10, 20, ..., 400
    T2_values = collect(0:10:400)  # 41 values: 0, 10, 20, ..., 400

    if T1_index < 1 || T1_index > length(T1_values)
        error("T1_INDEX must be between 1 and $(length(T1_values))")
    end
    
    T1 = T1_values[T1_index]
    
    # Get parameters from environment
    D = haskey(ENV, "CCR_D_VALUE") ? parse(Float64, ENV["CCR_D_VALUE"]) : 0.04
    
    grid_size = haskey(ENV, "CCR_GRID_SIZE") ? parse(Int, ENV["CCR_GRID_SIZE"]) : 5
    biomass_min_exp = haskey(ENV, "CCR_BIOMASS_MIN_EXP") ? parse(Float64, ENV["CCR_BIOMASS_MIN_EXP"]) : -5.0
    biomass_max_exp = haskey(ENV, "CCR_BIOMASS_MAX_EXP") ? parse(Float64, ENV["CCR_BIOMASS_MAX_EXP"]) : 0.0
    
    target_tau = haskey(ENV, "CCR_TAU") ? parse(Float64, ENV["CCR_TAU"]) : 1e6
    analysis_cycles = haskey(ENV, "CCR_ANALYSIS_CYCLES") ? parse(Int, ENV["CCR_ANALYSIS_CYCLES"]) : 1
    extinction_threshold = haskey(ENV, "CCR_EXTINCTION_THRESHOLD") ? parse(Float64, ENV["CCR_EXTINCTION_THRESHOLD"]) : 1e-100

    strategy_pair = haskey(ENV, "CCR_STRATEGY_PAIR") ? ENV["CCR_STRATEGY_PAIR"] : "0-1"
    
    # Validate strategy pair (only 0-1 and 1-2 are supported)
    if strategy_pair == "0-1"
        strategy_indices = [0, 1]  # N vs T
    elseif strategy_pair == "1-2"
        strategy_indices = [1, 2]  # T vs P
    else
        error("Invalid CCR_MODEL_PAIR: $(strategy_pair). Must be '0-1' (N vs T) or '1-2' (T vs P)")
    end
    
    base_output = haskey(ENV, "CCR_OUTPUT_DIR") ? ENV["CCR_OUTPUT_DIR"] : "results"
    output_dir = joinpath(base_output, "D_$(D)_$(strategy_pair)")


    mkpath(output_dir)
    
    println("Parameters:")
    println("  T1_INDEX = $(T1_index), T1 = $(T1)")
    println("  T2 range: $(T2_values[1]) to $(T2_values[end]) ($(length(T2_values)) values)")
    println("  Grid: $(grid_size)×$(grid_size)")
    println("  Strategy pair: $(strategy_pair) (indices: $(strategy_indices))")
    
    # Generate grid values for initial biomass ( 10^-5 to 10^0)
    biomass_values = 10 .^ range(biomass_min_exp, biomass_max_exp, length=grid_size)
    
    println("\nProcessing T1=$(T1) with $(length(T2_values)) T2 values...")
    
    # Process each T2 value
    for (t2_idx, T2) in enumerate(T2_values)
        println("\n[$(t2_idx)/$(length(T2_values))] Processing T1=$(T1), T2=$(T2)")
        
        # Setup parameters
        params = IndependentParameters()
        params.T1 = T1
        params.T2 = T2
        params.D = D
        
        # Generate grid combinations for initial biomass
        all_grid_combinations = []
        
        if strategy_pair == "0-1"
            for (y_idx, strategy0_initial) in enumerate(biomass_values)
                for (x_idx, strategy1_initial) in enumerate(biomass_values)
                    push!(all_grid_combinations, (
                        x_idx = x_idx, y_idx = y_idx,
                        strategy0_initial = strategy0_initial, strategy1_initial = strategy1_initial,
                        strategy_indices = strategy_indices,
                        grid_id = (y_idx - 1) * grid_size + x_idx
                    ))
                end
            end
        elseif strategy_pair == "1-2"
            for (y_idx, strategy1_initial) in enumerate(biomass_values)
                for (x_idx, strategy2_initial) in enumerate(biomass_values)
                    push!(all_grid_combinations, (
                        x_idx = x_idx, y_idx = y_idx,
                        strategy1_initial = strategy1_initial, strategy2_initial = strategy2_initial,
                        strategy_indices = strategy_indices,
                        grid_id = (y_idx - 1) * grid_size + x_idx
                    ))
                end
            end
        end
        
        total_combinations = length(all_grid_combinations)
        results_array = Vector{Any}(undef, total_combinations)
        start_time_t2 = time()
        
        @threads for combo_idx in 1:total_combinations
            grid_point = all_grid_combinations[combo_idx]
            
            result_dict_to_store = Dict()
            start_time_point = time()
            
            try
                local initial_biomasses
                if strategy_pair == "0-1"
                    initial_biomasses = [grid_point.strategy0_initial, grid_point.strategy1_initial]
                elseif strategy_pair == "1-2"
                    initial_biomasses = [grid_point.strategy1_initial, grid_point.strategy2_initial]
                end
                
                # Run Simulation
                sim_result = run_independent_simulation_method_a_cycles(
                    params, strategy_indices;
                    target_tau=target_tau,
                    analysis_cycles=analysis_cycles,
                    initial_biomasses=initial_biomasses
                )
                
                if sim_result !== nothing
                    sol = sim_result["solution"]
                    
                    # Analysis
                    analysis_result = analyze_independent_simulation_outcome(
                        sol, params, strategy_indices;
                        extinction_threshold=extinction_threshold,
                        verbose=false
                    )
                    
                    if analysis_result !== nothing
                        # Comparison
                        comparison_result = compare_independent_models(
                            analysis_result, sol, params, strategy_indices; verbose=false
                        )
                        
                        # Last Cycle Variables
                        last_cycle_vars, cycle_time = extract_last_cycle_end_variables(
                            sol, params, strategy_indices; verbose=false
                        )
                        
                        result_dict_to_store = merge(
                            analysis_result,
                            comparison_result,
                            Dict("last_cycle_variables" => last_cycle_vars,
                                 "last_cycle_time" => cycle_time)
                        )
                    else
                        result_dict_to_store = Dict("outcome" => "Analysis Failed", "winner" => -998)
                    end
                else
                    result_dict_to_store = Dict("outcome" => "Simulation Failed", "winner" => -999)
                end
                
            catch e
                result_dict_to_store = Dict("outcome" => "Exception: $(e)", "winner" => -997)
            end
            
            elapsed_time = time() - start_time_point
            
            # Store initial biomass values along with results
            if strategy_pair == "0-1"
                results_array[combo_idx] = (
                    x_idx = grid_point.x_idx,
                    y_idx = grid_point.y_idx,
                    grid_id = grid_point.grid_id,
                    strategy0_initial = grid_point.strategy0_initial,
                    strategy1_initial = grid_point.strategy1_initial,
                    result_dict = result_dict_to_store,
                    elapsed_time = elapsed_time
                )
            elseif strategy_pair == "1-2"
                results_array[combo_idx] = (
                    x_idx = grid_point.x_idx,
                    y_idx = grid_point.y_idx,
                    grid_id = grid_point.grid_id,
                    strategy1_initial = grid_point.strategy1_initial,
                    strategy2_initial = grid_point.strategy2_initial,
                    result_dict = result_dict_to_store,
                    elapsed_time = elapsed_time
                )
            end
        end
        
        elapsed_t2 = time() - start_time_t2
        println("    Completed $(total_combinations) grid points in $(round(elapsed_t2, digits=1))s")
        
        # Save results
        save_detailed_results(results_array, T1, T2, params, strategy_indices, output_dir)
        save_summary_results(results_array, T1, T2, params, strategy_indices, 
                           grid_size, biomass_min_exp, biomass_max_exp, output_dir)
    end
    
    println("\n=== Completed all T2 values for T1=$(T1) ===")
end


function save_detailed_results(results_array, T1, T2, params, strategy_indices, output_dir)
    comparison_key = "$(strategy_indices[1])_vs_$(strategy_indices[2])"
    winner_key = "$(comparison_key)_winner"
    strategy_pair_str = "$(strategy_indices[1])_$(strategy_indices[2])"
    
    df = DataFrame()
    df[!, :grid_id] = [r.grid_id for r in results_array]
    df[!, :x_idx] = [r.x_idx for r in results_array]
    df[!, :y_idx] = [r.y_idx for r in results_array]
    
    # Initial conditions
    if strategy_indices == [0, 1]
        df[!, :strategy0_initial] = [r.strategy0_initial for r in results_array]
        df[!, :strategy1_initial] = [r.strategy1_initial for r in results_array]
    elseif strategy_indices == [1, 2]
        df[!, :strategy1_initial] = [r.strategy1_initial for r in results_array]
        df[!, :strategy2_initial] = [r.strategy2_initial for r in results_array]
    end
    
    # Params
    df[!, :T1] .= T1
    df[!, :T2] .= T2
    df[!, :D] .= params.D
    
    # Results
    df[!, :outcome] = [get(r.result_dict, "outcome", "Unknown") for r in results_array]
    df[!, :comparison_result] = [get(r.result_dict, comparison_key, "Unknown") for r in results_array]
    df[!, :comparison_winner] = [get(r.result_dict, winner_key, -999) for r in results_array]
    df[!, :elapsed_time] = [r.elapsed_time for r in results_array]
    
    # Per-strategy stats
    for strategy_idx in strategy_indices
        df[!, "strategy_$(strategy_idx)_survives"] = [
            get(get(r.result_dict, "individual_results", Dict()), strategy_idx, Dict()) |> x -> get(x, "survives", false) 
            for r in results_array
        ]
        df[!, "strategy_$(strategy_idx)_final_biomass"] = [
            get(get(r.result_dict, "individual_results", Dict()), strategy_idx, Dict()) |> x -> get(x, "final_biomass", NaN) 
            for r in results_array
        ]
        
        # Last cycle variables
        for (j, var_name) in enumerate(["S1", "S2", "X1", "M1", "E1", "X2", "M2", "E2", "c_minus", "c"])
            df[!, "strategy_$(strategy_idx)_$(var_name)"] = [
                get(get(r.result_dict, "last_cycle_variables", Dict()), strategy_idx, zeros(10))[j]
                for r in results_array
            ]
        end
    end

    filename = "independent_grid_detailed_strategies_$(strategy_pair_str)_T1_$(T1)_T2_$(T2)_D_$(params.D).csv"
    CSV.write(joinpath(output_dir, filename), df)
    println("    Detailed results saved: $(filename)")
end


function save_summary_results(results_array, T1, T2, params, strategy_indices, 
                            grid_size, biomass_min_exp, biomass_max_exp, output_dir)
    
    successful_results = filter(r -> !contains(get(r.result_dict, "outcome", "Failed"), "Failed") && 
                                    !contains(get(r.result_dict, "outcome", "Exception"), "Exception"), results_array)
    num_successful = length(successful_results)
    total_results = length(results_array)
    
    # Generate comparison keys
    comparison_key = "$(strategy_indices[1])_vs_$(strategy_indices[2])"
    winner_key = "$(comparison_key)_winner"
    
    """
    Classify outcome based on all grid point results.
    
    Classification logic:
    1. If all results are the same Strategy_X_higher_biomass → Strategy_X_higher_biomass (consensus = 1.0)
    2. If all results are Both_extinct → Both_extinct (consensus = 1.0)
    3. If all results are Nearly_equal_biomass → check biomass variation across initial conditions:
       - If biomass variation < 1e-10 → Stable_equal (unique steady state)
       - If biomass variation ≥ 1e-10 → Neutral_equal (steady state depends on initial conditions)
    4. If mixed results → Multi-stability (consensus < 1.0)
    """
    function classify_outcome(comp_key, win_key)
        if isempty(successful_results)
            return "No_Data", -999, 0.0
        end
        
        outcomes = [get(r.result_dict, comp_key, "Unknown") for r in successful_results]
        winners = [get(r.result_dict, win_key, -999) for r in successful_results]
        
        # Count each outcome type
        outcome_counts = Dict{String, Int}()
        for o in outcomes
            outcome_counts[o] = get(outcome_counts, o, 0) + 1
        end
        
        if isempty(outcome_counts)
            return "No_Data", -999, 0.0
        end
        
        # Get unique outcomes
        unique_outcomes = collect(keys(outcome_counts))
        
        # Case: All results are the same
        if length(unique_outcomes) == 1
            single_outcome = unique_outcomes[1]
            
            # Find winner for this outcome
            winner_for_outcome = [w for (o, w) in zip(outcomes, winners) if o == single_outcome]
            dominant_winner = !isempty(winner_for_outcome) ? winner_for_outcome[1] : -999
            
            # Special case: All Nearly_equal_biomass → determine Stable_equal vs Neutral_equal
            if single_outcome == "Nearly_equal_biomass"
                # Collect biomass values for the first strategy to check variation
                first_strategy_idx = strategy_indices[1]
                biomass_values = Float64[]
                for r in successful_results
                    individual_results = get(r.result_dict, "individual_results", Dict())
                    strategy_result = get(individual_results, first_strategy_idx, Dict())
                    biomass = get(strategy_result, "final_biomass", NaN)
                    if !isnan(biomass) && isfinite(biomass)
                        push!(biomass_values, biomass)
                    end
                end
                
                if length(biomass_values) >= 2
                    biomass_range = maximum(biomass_values) - minimum(biomass_values)
                    if biomass_range < 1e-10
                        return "Stable_equal", dominant_winner, 1.0
                    else
                        return "Neutral_equal", dominant_winner, 1.0
                    end
                else
                    # Not enough data to determine, default to Stable_equal
                    return "Stable_equal", dominant_winner, 1.0
                end
            else
                # Both_extinct or Strategy_X_higher_biomass
                return single_outcome, dominant_winner, 1.0
            end
        else
            # Case: Mixed results → Multi-stability
            # Calculate consensus as the fraction of the most common outcome
            sorted_outcomes = sort(collect(outcome_counts), by=x->x[2], rev=true)
            dominant_outcome, count = sorted_outcomes[1]
            consensus_value = count / num_successful
            
            # Find winner for dominant outcome
            winner_for_dominant = [w for (o, w) in zip(outcomes, winners) if o == dominant_outcome]
            dominant_winner = !isempty(winner_for_dominant) ? winner_for_dominant[1] : -999
            
            return "Multi-stability", dominant_winner, consensus_value
        end
    end
    
    # Calculate biomass variance for each strategy (to assess sensitivity to initial conditions)
    strategy_biomass_variances = Dict{Int, Float64}()
    for strategy_idx in strategy_indices
        biomass_values = []
        for r in successful_results
            individual_results = get(r.result_dict, "individual_results", Dict())
            strategy_result = get(individual_results, strategy_idx, Dict())
            biomass = get(strategy_result, "final_biomass", NaN)
            if !isnan(biomass) && isfinite(biomass)
                push!(biomass_values, biomass)
            end
        end
        
        if length(biomass_values) > 1
            strategy_biomass_variances[strategy_idx] = var(biomass_values)
        else
            strategy_biomass_variances[strategy_idx] = NaN
        end
    end
    
    # Get representative last cycle variables (from first successful result)
    representative_cycle_vars = Dict{Int, Vector{Float64}}()
    if !isempty(successful_results)
        rep_result = successful_results[1]
        representative_cycle_vars = get(rep_result.result_dict, "last_cycle_variables", Dict())
    end
    
    # Classify outcome
    final_outcome, win_result, consensus_value = classify_outcome(comparison_key, winner_key)
    
    # Create summary DataFrame
    strategy_pair_str = "$(strategy_indices[1])_$(strategy_indices[2])"
    summary_df = DataFrame(
        T1 = [T1], 
        T2 = [T2], 
        D = [params.D], 
        K_CCR1 = [params.K_CCR1], 
        K_CCR2 = [params.K_CCR2],
        strategy_pair = ["$(strategy_indices[1])-$(strategy_indices[2])"],
        
        # Grid info
        grid_size = [grid_size],
        biomass_min_exp = [biomass_min_exp],
        biomass_max_exp = [biomass_max_exp],
        total_grid_points = [total_results],
        successful_points = [num_successful],
        success_rate = [num_successful / total_results],
        
        # Comparison results
        outcome = [final_outcome],
        winner = [win_result],
        consensus_value = [consensus_value],
        
        timestamp = [string(now())]
    )
    
    # Add biomass variances
    for strategy_idx in strategy_indices
        summary_df[!, "strategy_$(strategy_idx)_biomass_variance"] = [get(strategy_biomass_variances, strategy_idx, NaN)]
    end
    
    # Add representative last cycle variables
    for strategy_idx in strategy_indices
        if haskey(representative_cycle_vars, strategy_idx)
            vars = representative_cycle_vars[strategy_idx]
            for (i, var_name) in enumerate(["S1", "S2", "X1", "M1", "E1", "X2", "M2", "E2", "c_minus", "c"])
                summary_df[!, "strategy_$(strategy_idx)_$(var_name)_cycle_end"] = [length(vars) >= i ? vars[i] : NaN]
            end
        else
            for var_name in ["S1", "S2", "X1", "M1", "E1", "X2", "M2", "E2", "c_minus", "c"]
                summary_df[!, "strategy_$(strategy_idx)_$(var_name)_cycle_end"] = [NaN]
            end
        end
    end
    
    # Save summary file
    filename = "independent_grid_summary_strategies_$(strategy_pair_str)_T1_$(T1)_T2_$(T2)_D_$(params.D).csv"
    filepath = joinpath(output_dir, filename)
    CSV.write(filepath, summary_df)
    println("    Summary results saved: $(filename)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end