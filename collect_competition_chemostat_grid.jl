#!/usr/bin/env julia

# ================================================================================
# FILE: collect_competition_chemostat_grid.jl
#
# PURPOSE: Collect and combine competition chemostat grid SUMMARY results
#
# ENVIRONMENT VARIABLES:
#   CCR_D_VALUE        - Dilution rate (default: 0.07)
#   CCR_STRATEGY_PAIR  - Strategy pair (default: "1-2")
#   CCR_GRID_SIZE      - Grid size (default: 10)
#   CCR_INPUT_SUFFIX   - Optional suffix for input directory name
#                        e.g., "alpha1_10" → reads from competition_pair_1-2_D_0.04_grid9_alpha1_10
#   CCR_OUTPUT_SUFFIX  - Optional suffix for output filename
#                        e.g., "alpha1_10" → outputs competition_summary_pair_12_D_0.04_grid9_alpha1_10_combined.csv
#                        (If not set but CCR_INPUT_SUFFIX is set, uses CCR_INPUT_SUFFIX as output suffix too)
# ================================================================================

using CSV
using DataFrames
using Statistics
using Dates

"""
Extract parameters from competition SUMMARY filename.
Expected format: competition_strategies_X_Y_T1_A.A_T2_B.B_D_C.C_gridN_pair_XY_summary.csv
"""
function extract_summary_params_from_filename(filename::String)
    # Default values
    T1, T2, D = NaN, NaN, NaN
    strategies = ""
    grid_size = NaN
    strategy_pair = ""
    
    # Extract strategy indices
    strategy_match = match(r"strategies_([0-9_]+)_T1", filename)
    if strategy_match !== nothing
        strategies = strategy_match.captures[1]
    end
    
    # Extract T1
    t1_match = match(r"T1_([0-9]+\.?[0-9]*)", filename)
    if t1_match !== nothing
        T1 = parse(Float64, t1_match.captures[1])
    end
    
    # Extract T2
    t2_match = match(r"T2_([0-9]+\.?[0-9]*)", filename)
    if t2_match !== nothing
        T2 = parse(Float64, t2_match.captures[1])
    end
    
    # Extract D
    d_match = match(r"D_([0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)", filename)
    if d_match !== nothing
        D = parse(Float64, d_match.captures[1])
    end
    
    # Extract grid size
    grid_match = match(r"grid([0-9]+)_pair", filename)
    if grid_match !== nothing
        grid_size = parse(Int, grid_match.captures[1])
    end
    
    # Extract strategy pair
    pair_match = match(r"pair_([0-9-]+)_summary", filename)
    if pair_match !== nothing
        strategy_pair = pair_match.captures[1]
    end
    
    return T1, T2, D, strategies, grid_size, strategy_pair
end

"""
Collect all competition SUMMARY files from input directory.
"""
function collect_competition_summary_results(input_dir::String, output_file::String)
    println("=== Competition Summary Collector ===")
    println("Input directory: $(input_dir)")
    println("Output file: $(output_file)")
    
    # Find all competition SUMMARY CSV files
    all_files = readdir(input_dir)
    summary_files = filter(x -> contains(x, "competition_strategies_") && 
                               contains(x, "_summary.csv") && 
                               contains(x, "_pair_"), all_files)
    
    if length(summary_files) == 0
        println("Available files in directory (first 20):")
        for f in all_files[1:min(20, length(all_files))]
            println("  $(f)")
        end
        error("No competition summary files found in $(input_dir)")
    end
    
    println("Found $(length(summary_files)) competition summary files")
    
    # Initialize combined DataFrame
    combined_df = DataFrame()
    successful_count = 0
    failed_count = 0
    total_parameter_combinations = 0
    
    # Process each summary file
    for (i, filename) in enumerate(summary_files)
        file_path = joinpath(input_dir, filename)
        
        try
            # Read the summary file
            df = CSV.read(file_path, DataFrame)
            
            if nrow(df) > 0
                # Extract parameters from filename
                T1_file, T2_file, D_file, strategies_str, grid_size, strategy_pair = extract_summary_params_from_filename(filename)
                
                # Add filename-derived parameters
                df[!, :filename_T1] = fill(T1_file, nrow(df))
                df[!, :filename_T2] = fill(T2_file, nrow(df))
                df[!, :filename_D] = fill(D_file, nrow(df))
                df[!, :filename_strategies] = fill(strategies_str, nrow(df))
                df[!, :filename_grid_size] = fill(grid_size, nrow(df))
                df[!, :filename_strategy_pair] = fill(strategy_pair, nrow(df))
                df[!, :source_filename] = fill(filename, nrow(df))
                
                # Add analysis metadata
                df[!, :analysis_method] = fill("volume_grid_deterministic_with_nes", nrow(df))
                df[!, :data_type] = fill("parameter_combination_summary", nrow(df))
                
                # For first successful file, initialize structure
                if successful_count == 0
                    combined_df = copy(df)
                    println("Columns from first summary file ($(length(names(df))) total):")
                    for col in names(df)
                        println("  $(col)")
                    end
                else
                    # Append to combined DataFrame (handling column differences)
                    append!(combined_df, df, cols=:union)
                end
                
                successful_count += 1
                total_parameter_combinations += nrow(df)
                
                # Progress indicator
                if i % 50 == 0 || i == length(summary_files)
                    println("Processed $(i)/$(length(summary_files)) files... ($(total_parameter_combinations) parameter combinations so far)")
                end
                
            else
                println("Warning: Empty summary file $(filename)")
                failed_count += 1
            end
            
        catch e
            println("Warning: Failed to read $(filename): $(e)")
            failed_count += 1
        end
    end
    
    println("\nCollection Results:")
    println("Successfully processed: $(successful_count) summary files")
    println("Failed to process: $(failed_count) summary files")
    println("Total parameter combinations: $(total_parameter_combinations)")
    
    if nrow(combined_df) == 0
        error("No summary data was successfully collected!")
    end
    
    # Sort by parameters for better organization
    sort_cols = []
    if "T1" in names(combined_df); push!(sort_cols, :T1); end
    if "T2" in names(combined_df); push!(sort_cols, :T2); end
    if "D" in names(combined_df); push!(sort_cols, :D); end
    if "strategy_pair" in names(combined_df); push!(sort_cols, :strategy_pair); end
    
    if !isempty(sort_cols)
        sort!(combined_df, sort_cols)
    end
    
    # Create output directory if needed
    output_dir = dirname(output_file)
    if !isempty(output_dir) && !isdir(output_dir)
        mkpath(output_dir)
    end
    
    # Save combined data
    CSV.write(output_file, combined_df)
    
    println("\n=== Collection Summary ===")
    println("Total summary records: $(nrow(combined_df))")
    println("Total columns: $(ncol(combined_df))")
    println("Combined summary data saved to: $(output_file)")
    
    # Show parameter ranges
    if "T1" in names(combined_df)
        valid_T1 = filter(!isnan, combined_df.T1)
        valid_T2 = filter(!isnan, combined_df.T2)
        valid_D = filter(!isnan, combined_df.D)
        
        if !isempty(valid_T1)
            T1_values = sort(unique(valid_T1))
            println("T1 range: $(minimum(T1_values)) to $(maximum(T1_values)) ($(length(T1_values)) values)")
        end
        if !isempty(valid_T2)
            T2_values = sort(unique(valid_T2))
            println("T2 range: $(minimum(T2_values)) to $(maximum(T2_values)) ($(length(T2_values)) values)")
        end
        if !isempty(valid_D)
            println("D values: $(sort(unique(valid_D)))")
        end
    end
    
    println("\nCompetition summary collection completed!")
    return combined_df
end

"""
Get environment parameters.
"""
function get_env_params()
    D = haskey(ENV, "CCR_D_VALUE") ? parse(Float64, ENV["CCR_D_VALUE"]) : 0.07
    strategy_pair = haskey(ENV, "CCR_STRATEGY_PAIR") ? strip(replace(ENV["CCR_STRATEGY_PAIR"], r"[\'\"]" => "")) : "1-2"
    grid_size = haskey(ENV, "CCR_GRID_SIZE") ? parse(Int, ENV["CCR_GRID_SIZE"]) : 10
    input_suffix = haskey(ENV, "CCR_INPUT_SUFFIX") ? strip(ENV["CCR_INPUT_SUFFIX"]) : ""
    output_suffix = haskey(ENV, "CCR_OUTPUT_SUFFIX") ? strip(ENV["CCR_OUTPUT_SUFFIX"]) : input_suffix
    return D, strategy_pair, grid_size, input_suffix, output_suffix
end

"""
Main function.

Usage:
  julia collect_competition_chemostat_grid.jl <input_dir> <output_file>
  julia collect_competition_chemostat_grid.jl   # fallback to environment variables
"""
function main()
    if length(ARGS) >= 2
        # --- Mode 1: Command line arguments (preferred) ---
        input_dir = ARGS[1]
        output_file = ARGS[2]
    elseif length(ARGS) == 1
        error("Usage: julia collect_competition_chemostat_grid.jl <input_dir> <output_file>")
    else
        # --- Mode 2: Environment variables (backward compatible) ---
        D, strategy_pair, grid_size, input_suffix, output_suffix = get_env_params()
        
        println("Using environment variables:")
        println("  D = $(D), Strategy pair = $(strategy_pair), Grid size = $(grid_size)")
        
        base_path = haskey(ENV, "CCR_OUTPUT_DIR") ? ENV["CCR_OUTPUT_DIR"] : "results"

        input_dir_name = "competition_pair_$(strategy_pair)_D_$(D)_grid$(grid_size)"
        if !isempty(input_suffix)
            input_dir_name *= "_$(input_suffix)"
        end
        input_dir = joinpath(base_path, input_dir_name)
        
        output_dir = haskey(ENV, "CCR_COLLECT_OUTPUT_DIR") ? ENV["CCR_COLLECT_OUTPUT_DIR"] : "collected_results"

        mkpath(output_dir)
        pair_str = replace(strategy_pair, "-" => "")
        output_filename = "competition_summary_pair_$(pair_str)_D_$(D)_grid$(grid_size)"
        if !isempty(output_suffix)
            output_filename *= "_$(output_suffix)"
        end
        output_filename *= "_combined.csv"
        output_file = joinpath(output_dir, output_filename)
    end
    
    if !isdir(input_dir)
        error("Input directory not found: $(input_dir)")
    end
    
    try
        combined_df = collect_competition_summary_results(input_dir, output_file)
        println("\nCompetition summary collection completed!")
        return 0
    catch e
        println("Error: $(e)")
        return 1
    end
end

# Run main function when executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end