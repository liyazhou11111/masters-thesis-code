# Competition Chemostat

Chemostat competition simulations comparing the comparative advantages of transcriptional and post-transcriptional carbon catabolite repression (CCR) in a two-carbon-source environment with periodic substrate switching.

Three strategies are defined that differ solely in the level of CCR regulation:

| Index | Strategy | Description |
|-------|----------|-------------|
| 0 | No CCR (N) | No CCR regulation  |
| 1 | Transcriptional CCR (T) | CCR regulation on transcriptional level |
| 2 | Post-transcriptional CCR (P) | CCR regulation on post-transcriptional level |

Two competition experiments are performed:
- **No CCR vs Transcriptional CCR (0-1)**: Verifies that CCR confers a comparative advantage under tested conditions, using the initial parameter set.
- **Transcriptional CCR vs Post-transcriptional CCR (1-2)**: Compares whether different regulatory levels exhibit different comparative advantages across environments, using the default parameter set.

## Files

### Core model

- **`ccr_competition_chemostat_definitions.jl`** — Model definitions: parameter struct, 18-dimensional ODE system (2 substrates + 2×8 per-strategy variables), Ne×s neutral evolution analysis, and competition outcome classification.

### Model variants

- **`competition_model_constant_resource.jl`** — Fixed intracellular resource pool (c⁻ = 0.1). Eliminates the resource conservation effect.
- **`competition_model_accelerated_mRNA.jl`** — Accelerated mRNA dynamics (K_tau = 100). Eliminates the response speed difference between transcriptional and post-transcriptional CCR.

### HPC workflow (NeSI/SLURM)

- **`submit_competition_chemostat_grid.jl`** — Generates and submits SLURM job arrays. Each job runs one (T1, T2) parameter combination with a grid of initial volume fractions.
- **`run_competition_chemostat_grid.jl`** — Executed by each SLURM job. Runs competition simulations for a single (T1, T2) with `grid_size` initial conditions and outputs summary CSV files.
- **`collect_competition_chemostat_grid.jl`** — Collects individual job summary CSVs into a single combined file.

### Visualization

- **`Competition_chemostat_visualization.jl`** — Generates phase diagrams (T1 vs T2 plot) from combined summary CSV files. Supports multi-subplot layouts for comparing across D values or parameter variations.

## Environment Variables

All scripts read parameters from environment variables with defaults:

```bash
# Example configuration
export CCR_D_VALUE=0.04                   # Dilution rate
export CCR_STRATEGY_PAIR="1-2"            # Strategy pair: "0-1", "1-2"
export CCR_GRID_SIZE=9                    # Number of initial volume fractions
export CCR_VOLUME_MIN=1.0e-5              # Minimum volume fraction
export CCR_NE=1.0e8                       # Effective population size
export CCR_TAU=1.0e7                      # Target simulation time
export CCR_ANALYSIS_CYCLES=50             # Cycles for competition analysis
export CCR_EXTINCTION_THRESHOLD=1.0e-100  # Extinction threshold
export CCR_REGRESSION_CYCLES=50           # Cycles for fitness regression
export CCR_INITIAL_CONDITIONS_FILE="path/to/initial_conditions.csv"
export CCR_USE_FILE_CONDITIONS=true       # Use file-based initial conditions

# Directory configuration
export CCR_PROJECT_DIR="competition"      # Base project directory
export CCR_OUTPUT_DIR="results"           # Output directory for results
```

## Workflow

1. **Submit jobs**: `julia submit_competition_chemostat_grid.jl`
2. **Jobs run on cluster**: Each SLURM job executes `run_competition_chemostat_grid.jl`
3. **Collect results**: `julia collect_competition_chemostat_grid.jl <input_dir> <output_file>`
4. **Visualize**: Edit `csv_files` in `Competition_chemostat_visualization.jl` and run it

## Dependencies

- Julia 1.11+
- DifferentialEquations, DiffEqCallbacks, CSV, DataFrames, Statistics, LinearAlgebra, GLM, Printf, Dates
- Plots, LaTeXStrings, ColorSchemes (visualization only)
