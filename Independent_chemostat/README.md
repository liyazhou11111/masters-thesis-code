# Independent Chemostat

Independent cultivation simulations of CCR strategies in separate chemostats under identical environmental conditions.

This module serves two purposes:
1. **Providing initial conditions for co-culture competition**: Each strategy is run independently for a long period, and the endpoint state variables are extracted as initial conditions for the competition chemostat simulations (used for all competition sets).
2. **Biomass comparison**: Compares the average biomass over the last complete cycle between strategies grown independently, to test whether co-culture competition outcomes can be predicted from independent cultivation (performed for No CCR vs Transcriptional CCR only, using the initial parameter set).

## Files

### Model

- **`independent_chemostat_definitions.jl`** — Model definitions: parameter struct, 10-dimensional ODE system per strategy (S1, S2, X1, M1, E1, X2, M2, E2, c⁻, c), simulation runner, biomass analysis, and endpoint state extraction.

### HPC workflow (NeSI/SLURM)

- **`submit_independent_chemostat_grid.jl`** — Generates and submits SLURM job arrays. 
- **`run_independent_chemostat_grid.jl`** — Executed by each SLURM job. 
- **`collect_independent_chemostat_summaries.jl`** — Collects individual job summary CSVs into a single combined file.

### Visualization

- **`Independent_chemostat_visualization.jl`** — Generates comparison diagram showing which strategy achieves higher biomass. Supports multi-subplot layouts for comparing across D values.



## Workflow

1. **Submit jobs**: `julia submit_independent_chemostat_grid.jl`
2. **Jobs run on cluster**: Each SLURM job executes `run_independent_chemostat_grid.jl`
3. **Collect results**: `julia collect_independent_chemostat_summaries.jl <input_dir>`
4. **Visualize**: Edit `csv_files` in `Independent_chemostat_visualization.jl` and run it

## Environment Variables

```bash
# Example configuration
export CCR_D_VALUE=0.02                   # Dilution rate
export CCR_STRATEGY_PAIR="0-1"            # Strategy pair: "0-1" or "1-2"
export CCR_GRID_SIZE=5                    # Grid size for initial biomass scan
export CCR_BIOMASS_MIN_EXP=-5.0           # Min initial biomass exponent (10^x)
export CCR_BIOMASS_MAX_EXP=0.0            # Max initial biomass exponent (10^x)
export CCR_TAU=1.0e6                      # Target simulation time (generations)
export CCR_ANALYSIS_CYCLES=50             # Cycles for biomass averaging
export CCR_EXTINCTION_THRESHOLD=1.0e-100  # Extinction threshold

# Directory configuration
export CCR_PROJECT_DIR="independent_grid" # Base project directory
export CCR_OUTPUT_DIR="results"           # Output directory for results
```


