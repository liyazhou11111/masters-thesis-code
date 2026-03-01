# Transcriptional and Post-transcriptional CCR Competition Model

Simulation code for the study: *Transcriptional and post-transcriptional carbon catabolite repression differentially manage the trade-off between high growth rate and fast adaptive capacity*.

This project investigates whether transcriptional and post-transcriptional carbon catabolite reression (CCR) differentially manage the trade-off between high growth rate and fast adaptive capacity through their respective properties of resource conservation and rapid response, thereby conferring comparative advantages in different environments. To this end, we developed a chemostat-based mathematical model that includes three CCR regulatory strategies:

| Index | Strategy | Description |
|-------|----------|-------------|
| 0 | No CCR (N) | No CCR regulation |
| 1 | Transcriptional CCR (T) | CCR regulation on transcriptional level |
| 2 | Post-transcriptional CCR (P) | CCR regulation on post-transcriptional level |

The model covers gene expression, nutrient uptake and metabolism, biosynthesis and degradation, gene induction, and cell growth dilution, adapted from [Narang and Pilyugin (2007)](https://doi.org/10.1016/j.jtbi.2006.08.007).

## Project Structure

### Batch/

Batch culture simulations to verify that the model reproduces the diauxic growth phenotype of CCR. Converts the chemostat model to a closed batch system (D=0, no inflow).

### Independent_chemostat/

Independent cultivation of each strategy in separate chemostats under identical environmental conditions. Provides initial conditions for all co-culture competition simulations and compares biomass between strategies grown independently.

### Competition_chemostat/

Co-culture competition simulations. Two strategies are placed in the same chemostat with periodic carbon source switching (S1 pulsed, S2 constant), and competition outcomes are classified across a grid of environmental conditions (T1, T2) to determine whether CCR confers a comparative advantage.

### Chemostat_competition/

Time-series analysis and visualization of co-culture competition between Transcriptional CCR and Post-transcriptional CCR, primarily used to determine whether the comparative advantage differences arise from the differential management of the trade-off between high growth rate and fast adaptive capacity through resource conservation and rapid response properties. Generates detailed plots (P1-P5) of biomass dynamics, state variables, variable differences, specific growth rate decomposition, and growth rate differences for the last complete cycle, and QuadGK integration analysis of synthesis rates and resource investment ratio.
