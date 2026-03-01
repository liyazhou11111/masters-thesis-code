#!/usr/bin/env julia

# ================================================================================
# FILE: collect_independent_chemostat_summaries.jl
# PURPOSE: Aggregates summary CSV results from independent chemostat grid analysis
# ================================================================================

using CSV
using DataFrames
using Dates

function collect_summary_files(data_dir::String)
    """Collection of all summary CSV files for specific parameters"""
    
    # Get target parameters from environment variables (or defaults)
    TARGET_D = haskey(ENV, "CCR_D_VALUE") ? parse(Float64, ENV["CCR_D_VALUE"]) : 0.04
    strategy_pair = haskey(ENV, "CCR_STRATEGY_PAIR") ? ENV["CCR_STRATEGY_PAIR"] : "0-1"
    
    println("=== Collecting Independent Chemostat Grid Summaries ===")
    println("Source Directory: $(data_dir)")
    println("Target Parameters:")
    println("  D = $(TARGET_D)")
    println("  Strategy Pair = $(strategy_pair)")
    
    # Find summary files
    if !isdir(data_dir)
        error("Directory does not exist: $(data_dir)")
    end

    all_files = readdir(data_dir, join=true)
    summary_files = filter(f -> contains(basename(f), "independent_grid_summary_strategies_") && 
                              endswith(f, ".csv"), all_files)
    
    println("Found $(length(summary_files)) potential summary files.")
    
    if isempty(summary_files)
        println("❌ No summary files found. Check if simulations finished.")
        return nothing
    end

    # Collect data
    all_data = DataFrame[]
    valid_count = 0
    
    for file in summary_files
        try
            df = CSV.read(file, DataFrame)
            
            # Basic validation: check if DataFrame is not empty
            if nrow(df) > 0
                # Check if D value matches target (floating point comparison with tolerance)
                file_D = df.D[1]
                
                if abs(file_D - TARGET_D) < 1e-9
                    push!(all_data, df)
                    valid_count += 1
                    # Print progress every 50 files
                    if valid_count % 50 == 0
                        print(".")
                    end
                end
            end
        catch e
            println("\n❌ Failed to read: $(basename(file)) - $e")
        end
    end
    println("\nSuccessfully collected $(valid_count) valid files matching target parameters.")
    
    if isempty(all_data)
        error("No files matched the target parameters (D=$(TARGET_D)).")
    end
    
    # Combine all dataframes
    combined_df = vcat(all_data...)
    
    # Sort by T1 then T2 for cleaner output (useful for heatmaps)
    sort!(combined_df, [:T1, :T2])
    
    # Define Output Directory
    base_output_dir = haskey(ENV, "CCR_COLLECT_OUTPUT_DIR") ? ENV["CCR_COLLECT_OUTPUT_DIR"] : "collected_results"
    mkpath(base_output_dir)
    
    # Create filename
    output_filename = "independent_grid_combined_summary_D_$(TARGET_D)_$(strategy_pair).csv"
    output_file = joinpath(base_output_dir, output_filename)
    
    # Save
    CSV.write(output_file, combined_df)
    
    println("\n" * "="^60)
    println("✅ COMBINATION COMPLETE")
    println("Total rows: $(nrow(combined_df))")
    println("Saved to: $(output_file)")
    println("="^60)
    
    return combined_df
end

function main()
    # Determine input directory
    if length(ARGS) >= 1
        data_dir = ARGS[1]
    else
        # Auto-construct path from environment variables if not provided as argument
        D_str = haskey(ENV, "CCR_D_VALUE") ? ENV["CCR_D_VALUE"] : "0.04"
        D_val = parse(Float64, D_str)
        
        strategy_pair = haskey(ENV, "CCR_STRATEGY_PAIR") ? ENV["CCR_STRATEGY_PAIR"] : "0-1"
        
        # Construct path based on run_independent_chemostat_grid.jl logic
        base_path = haskey(ENV, "CCR_OUTPUT_DIR") ? ENV["CCR_OUTPUT_DIR"] : "results"
        folder_name = "D_$(D_val)_$(strategy_pair)"
        data_dir = joinpath(base_path, folder_name)
        
        println("No directory argument provided.")
        println("Attempting to infer path from environment: $(data_dir)")
    end
    
    collect_summary_files(data_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end