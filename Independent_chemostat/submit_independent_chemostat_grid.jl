#!/usr/bin/env julia

# ================================================================================
# FILE: submit_independent_chemostat_grid.jl
# PURPOSE: SLURM job array submitter for independent chemostat simulations
# 
# This script generates and submits a SLURM job array that runs independent
# chemostat simulations across a parameter sweep of T1 values.
# ================================================================================

using Printf
using Dates

function main()
    println("=== CCR Independent Chemostat Grid Sweep Job Array Submitter ===")
    println("Generated on: $(now())")

    # Parameter ranges
    T1_values = haskey(ENV, "CCR_T1_RANGE") ? eval(Meta.parse(ENV["CCR_T1_RANGE"])) : collect(0:5:100)
    T2_values = haskey(ENV, "CCR_T2_RANGE") ? eval(Meta.parse(ENV["CCR_T2_RANGE"])) : collect(200:10:400)

    D_value = haskey(ENV, "CCR_D_VALUE") ? parse(Float64, ENV["CCR_D_VALUE"]) : 0.04
    
    grid_size = haskey(ENV, "CCR_GRID_SIZE") ? parse(Int, ENV["CCR_GRID_SIZE"]) : 5
    biomass_min_exp = haskey(ENV, "CCR_BIOMASS_MIN_EXP") ? parse(Float64, ENV["CCR_BIOMASS_MIN_EXP"]) : -5.0
    biomass_max_exp = haskey(ENV, "CCR_BIOMASS_MAX_EXP") ? parse(Float64, ENV["CCR_BIOMASS_MAX_EXP"]) : 0.0
    
    target_tau = haskey(ENV, "CCR_TAU") ? parse(Float64, ENV["CCR_TAU"]) : 1e6
    analysis_cycles = haskey(ENV, "CCR_ANALYSIS_CYCLES") ? parse(Int, ENV["CCR_ANALYSIS_CYCLES"]) : 1
    extinction_threshold = haskey(ENV, "CCR_EXTINCTION_THRESHOLD") ? parse(Float64, ENV["CCR_EXTINCTION_THRESHOLD"]) : 1e-100
    
    strategy_pair = haskey(ENV, "CCR_STRATEGY_PAIR") ? ENV["CCR_STRATEGY_PAIR"] : "0-1"
    
    # Validate strategy pair (only 0-1 and 1-2 are supported)
    valid_pairs = ["0-1", "1-2"]
    if !(strategy_pair in valid_pairs)
        error("Invalid CCR_STRATEGY_PAIR: $(strategy_pair). Must be '0-1' (N vs T) or '1-2' (T vs P)")
    end

    total_jobs = length(T1_values)
    
    # Setup directories
    base_project_dir = haskey(ENV, "CCR_PROJECT_DIR") ? ENV["CCR_PROJECT_DIR"] : "independent_grid"
    base_output_dir = haskey(ENV, "CCR_OUTPUT_DIR") ? ENV["CCR_OUTPUT_DIR"] : "results"

    
    sweep_folder_name = "independent_grid_strategies_$(strategy_pair)_D_$(D_value)"
    submission_dir = joinpath(base_project_dir, sweep_folder_name)
    results_dir = joinpath(base_output_dir, "D_$(D_value)_$(strategy_pair)")
    
    if isdir(submission_dir)
        println("WARNING: Directory exists: $(submission_dir)")
    end
    
    mkpath(submission_dir)
    mkpath(results_dir)
    mkpath(joinpath(submission_dir, "logs"))

    # Copy necessary files
    needed_files = [
        "independent_chemostat_definitions.jl",
        "run_independent_chemostat_grid.jl"
    ]
    
    for file in needed_files
        if isfile(file)
            cp(file, joinpath(submission_dir, file), force=true)
        else
            error("Required file not found: $(file)")
        end
    end

    # Create parameter file (T1 index mapping)
    param_file = joinpath(submission_dir, "T1_parameters.txt")
    open(param_file, "w") do f
        for (i, T1) in enumerate(T1_values)
            println(f, "$i $T1")
        end
    end

    # Print configuration summary
    println("\nConfiguration:")
    println("  Strategy pair: $(strategy_pair)")
    println("  T1 values: $(length(T1_values)) ($(T1_values[1]) to $(T1_values[end]), step=$(T1_values[2]-T1_values[1]))")
    println("  T2 values: $(length(T2_values)) ($(T2_values[1]) to $(T2_values[end]), step=$(T2_values[2]-T2_values[1]))")
    println("  D = $(D_value)")
    println("  Grid size: $(grid_size)×$(grid_size)")
    println("  Total jobs: $(total_jobs)")

    # Computational requirements
    time_limit = "15:00:00"
    concurrent_limit = 20
    memory_per_job = "12G"
    cpus_per_task = 6

    # Generate SLURM Script
    slurm_script = """#!/bin/bash
#SBATCH --job-name=CCR_Indep_$(strategy_pair)
#SBATCH --account=uoa04294
#SBATCH --time=$(time_limit)
#SBATCH --mem=$(memory_per_job)
#SBATCH --cpus-per-task=$(cpus_per_task)
#SBATCH --array=1-$(total_jobs)%$(concurrent_limit)
#SBATCH --partition=milan
#SBATCH --output=logs/slurm_%A_%a.out
#SBATCH --error=logs/slurm_%A_%a.err

module load Julia/1.11.3-GCC-12.3.0-VTune

export OPENBLAS_NUM_THREADS=1
    
# Export parameters
export CCR_D_VALUE=$(D_value)
export CCR_OUTPUT_DIR=$(results_dir)
export CCR_STRATEGY_PAIR="$(strategy_pair)"
export CCR_GRID_SIZE=$(grid_size)
export CCR_BIOMASS_MIN_EXP=$(biomass_min_exp)
export CCR_BIOMASS_MAX_EXP=$(biomass_max_exp)
export CCR_TAU=$(target_tau)
export CCR_ANALYSIS_CYCLES=$(analysis_cycles)
export CCR_EXTINCTION_THRESHOLD=$(extinction_threshold)

mkdir -p logs

PARAM_LINE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" T1_parameters.txt)
read -r T1_INDEX T1_VALUE <<< "\$PARAM_LINE"

echo "Running Independent Chemostat Grid"
echo "T1_INDEX: \$T1_INDEX, T1_VALUE: \$T1_VALUE"
echo "Strategy pair: $(strategy_pair)"

julia --threads=\$SLURM_CPUS_PER_TASK run_independent_chemostat_grid.jl \$T1_INDEX
"""
    
    script_file = joinpath(submission_dir, "submit_independent_grid.sh")
    open(script_file, "w") do f
        write(f, slurm_script)
    end
    chmod(script_file, 0o755)

    # Submit job array
    try
        original_dir = pwd()
        cd(submission_dir)
        result = read(`sbatch submit_independent_grid.sh`, String)
        println("\n✅ SUBMITTED SUCCESSFULLY!")
        println(strip(result))
        cd(original_dir)
    catch e
        println("\n❌ ERROR: Failed to submit job array: $(e)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end