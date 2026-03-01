# Model variant: Fixed intracellular resource pool (c⁻ = 0.1)
#
# Eliminates the resource conservation effect by replacing the dynamic c⁻
# in mRNA and enzyme synthesis rates (Eqs. 10, 11, 14, 15) with a constant 0.1.
# This removes resource reallocation as a mechanism for competitive advantage,
# isolating other factors that contribute to fitness differences between strategies.
# All other equations remain identical to the default model.

function competition_model!(du, u, p, t)
    params, strategy_indices = p

    S1 = u[1]
    S2 = u[2]

    # Substrate dynamics: inflow/dilution (consumption added per strategy below)
    dS1dt = (S1_in(t, params) - S1) * params.D
    dS2dt = (S2_in(t, params) - S2) * params.D

    for (i, strategy_idx) in enumerate(strategy_indices)
        offset = 2 + (i - 1) * 8
        X1      = u[offset + 1]
        M1      = u[offset + 2]
        E1      = u[offset + 3]
        X2      = u[offset + 4]
        M2      = u[offset + 5]
        E2      = u[offset + 6]
        c_minus = u[offset + 7]
        c       = u[offset + 8]

        # Growth rate (Eq. 20)
        rg = (params.k1 * E1 * S1 / (params.KS1 + S1) +
              params.k2 * E2 * S2 / (params.KS2 + S2)) -
             ((1 - params.Y1) * params.alpha1 * X1 +
              (1 - params.Y2) * params.alpha2 * X2)

        # CCR regulation functions (Eq. 12-13, 16-17)
        #   kappa_M: transcriptional CCR modifier for M2 synthesis
        #   kappa_E: post-transcriptional CCR modifier for E2 synthesis
        if strategy_idx == 0       # Strategy N: No CCR
            kappa_M = 1.0
            kappa_E = 1.0
        elseif strategy_idx == 1   # Strategy T: Transcriptional CCR (Eq. 13)
            kappa_M = params.K_CCR1 / (params.K_CCR1 + X1)
            kappa_E = 1.0
        elseif strategy_idx == 2   # Strategy P: Post-transcriptional CCR (Eq. 17)
            kappa_M = 1.0
            kappa_E = params.K_CCR2 / (params.K_CCR2 + X1)
        else
            error("Unknown strategy_idx: $strategy_idx. Supported: 0 (N), 1 (T), 2 (P).")
        end

        # Substrate consumption by this strategy
        dS1dt -= params.k1 * E1 * S1 / (params.KS1 + S1) * c
        dS2dt -= params.k2 * E2 * S2 / (params.KS2 + S2) * c

        # ODEs
        # dX1/dt (Eq. 8)
        du[offset + 1] = params.k1 * E1 * S1 / (params.KS1 + S1) -
                          params.alpha1 * X1 - rg * X1
        # dM1/dt (Eq. 10)
        du[offset + 2] = params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * 0.1 -
                          params.delta1 * M1 - rg * M1
        # dE1/dt (Eq. 14)
        du[offset + 3] = params.kE * M1 * 0.1 -
                          params.delta2 * E1 - rg * E1
        # dX2/dt (Eq. 9)
        du[offset + 4] = params.k2 * E2 * S2 / (params.KS2 + S2) -
                          params.alpha2 * X2 - rg * X2
        # dM2/dt (Eq. 11) — includes kappa_M for transcriptional CCR
        du[offset + 5] = params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * 0.1 * kappa_M -
                          params.delta1 * M2 - rg * M2
        # dE2/dt (Eq. 15) — includes kappa_E for post-transcriptional CCR
        du[offset + 6] = params.kE * M2 * 0.1 * kappa_E -
                          params.delta2 * E2 - rg * E2
        # dc_minus/dt (Eq. 18)
        du[offset + 7] = (params.Y1 * params.alpha1 * X1 + params.Y2 * params.alpha2 * X2) +
                          (params.delta1 * M1 + params.delta2 * E1 +
                           params.delta1 * M2 + params.delta2 * E2) -
                          (params.kM * ((X1 + params.b * params.KX1) / (params.KX1 + X1)) * c_minus +
                           params.kE * M1 * c_minus +
                           params.kM * ((X2 + params.b * params.KX2) / (params.KX2 + X2)) * c_minus * kappa_M +
                           params.kE * M2 * c_minus * kappa_E) -
                          rg * c_minus
        # dc/dt (Eq. 21)
        du[offset + 8] = (rg - params.D) * c

        # Non-negativity clipping
        for j in (offset+1):(offset+8)
            if u[j] <= 0.0 && du[j] < 0.0
                du[j] = 0.0
            end
        end
    end

    # Non-negativity for substrates
    if u[1] <= 0.0 && dS1dt < 0.0; dS1dt = 0.0; end
    if u[2] <= 0.0 && dS2dt < 0.0; dS2dt = 0.0; end

    du[1] = dS1dt
    du[2] = dS2dt

    return nothing
end