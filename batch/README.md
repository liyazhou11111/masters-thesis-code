# Batch Culture

Batch culture simulations for the three CCR strategies (No CCR, Transcriptional CCR, and Post-transcriptional CCR) to verify that the model reproduces the diauxic growth phenotype of carbon catabolite repression. The chemostat model is converted to a closed batch system by setting D=0 and eliminating carbon source inflow.

## Plot

- **Growth rate**: Biomass growth rate (dc/dt) and carbon source concentrations over time for each strategy, showing whether diauxic growth occurs.

## Files

- **`batch.jl`** — Self-contained module including: parameter struct, ODE systems for each strategy (N, T, P), batch simulation runner, and visualization functions.

## Usage

Run in Julia/VSCode: `include("batch.jl")`

Plots are displayed and saved as PDF to `plot_results/`.

## Dependencies

- Julia 1.11+
- DifferentialEquations, Printf
- Plots, Plots.PlotMeasures (visualization)
