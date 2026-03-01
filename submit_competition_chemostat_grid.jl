#!/usr/bin/env julia

# ================================================================================
# FILE: submit_competition_chemostat_grid.jl (Large Scale Parameter Sweep)
#
# PURPOSE: Submit 451 individual jobs, one per T1,T2 combination
# Each job processes one T1,T2 with grid_size² initial conditions
# ================================================================================

using Printf
using Dates

function main()
    println("=== CCR Large Scale Parameter Sweep Job Submitter ===")
    println("Generated on: $(now())")
    
    
    # Define parameter ranges (must match run script)
    
    T1_values = collect(range(0, 100, length=21))
    T2_values = collect(range(0, 400, length=41))
    

    
    # Fixed parameters from environment
    D_value = haskey(ENV, "CCR_D_VALUE") ? parse(Float64, ENV["CCR_D_VALUE"]) : 0.04
      
    grid_size = haskey(ENV, "CCR_GRID_SIZE") ? parse(Int, ENV["CCR_GRID_SIZE"]) : 9
    volume_min = haskey(ENV, "CCR_VOLUME_MIN") ? parse(Float64, ENV["CCR_VOLUME_MIN"]) : 1e-5
    
    Ne = haskey(ENV, "CCR_NE") ? parse(Float64, ENV["CCR_NE"]) : 1e8
    target_tau = haskey(ENV, "CCR_TAU") ? parse(Float64, ENV["CCR_TAU"]) : 1e7
    analysis_cycles = haskey(ENV, "CCR_ANALYSIS_CYCLES") ? parse(Int, ENV["CCR_ANALYSIS_CYCLES"]) : 50
    extinction_threshold = haskey(ENV, "CCR_EXTINCTION_THRESHOLD") ? parse(Float64, ENV["CCR_EXTINCTION_THRESHOLD"]) : 1e-100
    n_cycles_for_regression = haskey(ENV, "CCR_REGRESSION_CYCLES") ? parse(Int, ENV["CCR_REGRESSION_CYCLES"]) : 50
    
    strategy_pair = haskey(ENV, "CCR_STRATEGY_PAIR") ? strip(replace(ENV["CCR_STRATEGY_PAIR"], r"[\'\"]" => "")) : "1-2"
    
    # Initial conditions file support
    initial_conditions_file = haskey(ENV, "CCR_INITIAL_CONDITIONS_FILE") ? ENV["CCR_INITIAL_CONDITIONS_FILE"] : ""
    use_file_initial_conditions = haskey(ENV, "CCR_USE_FILE_CONDITIONS") ? parse(Bool, ENV["CCR_USE_FILE_CONDITIONS"]) : false

    # Create parameter combinations
    param_combinations = [(T1, T2) for T1 in T1_values for T2 in T2_values]
    total_jobs = length(param_combinations)
    
    # Setup directories
    base_project_dir = haskey(ENV, "CCR_PROJECT_DIR") ? ENV["CCR_PROJECT_DIR"] : "competition"
    base_output_dir = haskey(ENV, "CCR_OUTPUT_DIR") ? ENV["CCR_OUTPUT_DIR"] : "results"

    
    sweep_folder_name = "competition_sweep_pair_$(strategy_pair)_D_$(D_value)_grid$(grid_size)"
    submission_dir = joinpath(base_project_dir, sweep_folder_name)
    results_dir = joinpath(base_output_dir, "competition_pair_$(strategy_pair)_D_$(D_value)_grid$(grid_size)")

    if isdir(submission_dir)
        println("\nERROR: Directory already exists: $(submission_dir)")
        return
    end
    
    mkpath(submission_dir)
    mkpath(results_dir)
    mkpath(joinpath(submission_dir, "logs"))
    println("\nCreated directories...")

    # Copy necessary files
    needed_files = ["ccr_competition_chemostat_definitions.jl", "run_competition_chemostat_grid.jl"]
    for file in needed_files
        cp(file, joinpath(submission_dir, file), force=true)
    end
    println("Copied necessary files.")

    # Create parameter file mapping job indices to T1,T2 values
    param_file = joinpath(submission_dir, "parameter_combinations.txt")
    open(param_file, "w") do f
        for (i, (T1, T2)) in enumerate(param_combinations)
            println(f, "$i $T1 $T2")
        end
    end
    println("Created parameter file with $(total_jobs) combinations")

    
    # Computational requirements
    time_limit = "24:00:00"  # 24 hours per T1,T2 combination
    concurrent_limit = 80   # Allow concurrent jobs
    memory_per_job = "8G"
    cpus_per_task = 4
    
    
    println("\n=== SLURM Configuration ===")
    println("Time limit: $(time_limit)")
    println("Concurrent job limit: $(concurrent_limit)")
    println("Memory per job: $(memory_per_job)")
    println("CPUs per job: $(cpus_per_task)")
    println("Total jobs: $(total_jobs)")
    println("Volume fractions per job: $(grid_size)")
    println("Volume range: $(volume_min) to $(1-volume_min)")

    # Generate SLURM Script with updated environment variables
    slurm_script_name = "submit_large_scale_volume_sweep.sh"
    slurm_script = """#!/bin/bash
#SBATCH --job-name=CCR_Volume_Grid
#SBATCH --account=uoa04294
#SBATCH --time=$(time_limit)
#SBATCH --mem=$(memory_per_job)
#SBATCH --cpus-per-task=$(cpus_per_task)
#SBATCH --array=1-$(total_jobs)%$(concurrent_limit)
#SBATCH --partition=milan
#SBATCH --output=logs/slurm_%A_%a.out
#SBATCH --error=logs/slurm_%A_%a.err

module load Julia/1.11.3-GCC-12.3.0-VTune

# Prevent BLAS from over-threading (leave CPU cores to Julia threads)
export OPENBLAS_NUM_THREADS=1
    
# Export all necessary parameters
export CCR_D_VALUE=$(D_value)
export CCR_STRATEGY_PAIR="$(strategy_pair)"
export CCR_OUTPUT_DIR=$(results_dir)
export CCR_GRID_SIZE=$(grid_size)
export CCR_VOLUME_MIN=$(volume_min)
export CCR_NE=$(Ne)
export CCR_TAU=$(target_tau)
export CCR_ANALYSIS_CYCLES=$(analysis_cycles)
export CCR_EXTINCTION_THRESHOLD=$(extinction_threshold)
export CCR_REGRESSION_CYCLES=$(n_cycles_for_regression)
export CCR_INITIAL_CONDITIONS_FILE="$(initial_conditions_file)"
export CCR_USE_FILE_CONDITIONS=$(use_file_initial_conditions)

mkdir -p logs

# Get T1,T2 values for this job
PARAM_LINE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" parameter_combinations.txt)
read -r INDEX T1 T2 <<< "\$PARAM_LINE"

echo "=============================================="
echo "SLURM_JOB_ID: \$SLURM_JOB_ID"
echo "SLURM_ARRAY_TASK_ID: \$SLURM_ARRAY_TASK_ID"
echo "Running Volume Fraction Grid Sweep"
echo "Strategy pair: $(strategy_pair)"
echo "Parameters: T1=\$T1, T2=\$T2"
echo "Volume grid: $(grid_size) fractions from $(volume_min) to $(1-volume_min)"
echo "Using \$SLURM_CPUS_PER_TASK threads"
echo "Initial conditions: $(use_file_initial_conditions ? "From file" : "Default")"
echo "=============================================="

julia --threads=\$SLURM_CPUS_PER_TASK run_competition_chemostat_grid.jl \$T1 \$T2

echo "Job finished with exit code \$?"
"""
    
    script_file = joinpath(submission_dir, slurm_script_name)
    open(script_file, "w") do f
        write(f, slurm_script)
    end
    chmod(script_file, 0o755)
    println("\nCreated SLURM array script: $(script_file)")

    println("\n=== PARAMETER SWEEP SUMMARY ===")
    println("T1 range: $(T1_values[1]) to $(T1_values[end]) ($(length(T1_values)) values)")
    println("T2 range: $(T2_values[1]) to $(T2_values[end]) ($(length(T2_values)) values)")
    println("Strategy pair: $(strategy_pair)")
    println("Total parameter combinations: $(total_jobs)")
    println("Volume fractions per job: $(grid_size)")
    println("Volume range: $(volume_min) to $(1-volume_min)")
    println("Total analysis points: $(total_jobs * grid_size)")
    println("Expected runtime: $(total_jobs) jobs × up to $(time_limit) each")

    # Submit the job array
    try
        original_dir = pwd()
        cd(submission_dir)
        result = read(`sbatch $(slurm_script_name)`, String)
        println("\n" * "="^70)
        println("✅ VOLUME FRACTION GRID SWEEP JOB ARRAY SUBMITTED!")
        println("="^70)
        println(strip(result))
        println("\nJob array details:")
        println("• Submission directory: $(submission_dir)")
        println("• Results directory: $(results_dir)")
        println("• Log files: $(joinpath(submission_dir, "logs"))")
        println("• Expected output: $(total_jobs) pairs of CSV files")
        cd(original_dir)
    catch e
        println("\n" * "="^70)
        println("❌ ERROR: Failed to submit job array")
        println("="^70)
        println("Error message: $(e)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end