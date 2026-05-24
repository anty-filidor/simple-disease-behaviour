function cond(x, I)
    if x >= 0 && x <= 1 && I >= 0 && I <= 1 && I > 0
        return true
    else
        return false
    end
end

function dynamics_rule_si(x0, param, t)
    R, p, c = param  # R: basic reproduction no; p: imitation rate; c: spontaneous reversion rate
    x, I = x0        # x: fraction with risky behaviour ∈ [0,1]; I: infected fraction ∈ [0,1]

    # Behaviour dynamics: two-sample pairwise imitation. An agent copies the opposite strategy only
    # if BOTH randomly sampled opponents hold it.
    risky_to_protective = p * (1+I) * x * (1-x)^2  # (1+I) weight: infection raises urgency; (1-x)^2: both sampled agents are protective
    protective_to_risky = p * (1-I) * (1-x) * x^2  # (1-I) weight: low prevalence breeds complacency; x^2: both sampled agents are risky
    spontaneous = c * (1-I) * (1-x)                # spontaneous reversion to risky at rate c, suppressed by high prevalence
    dxdt = -risky_to_protective + protective_to_risky + spontaneous

    # Epidemic dynamics: standard SIS (no permanent immunity).
    new_infections = R * x * I * (1-I)  # only risky agents (fraction x) transmit; R scales overall infectivity
    recoveries = I                      # recovery at normalised rate 1
    dIdt = new_infections - recoveries

    r = SVector{2}([dxdt; dIdt])
    return r
end
