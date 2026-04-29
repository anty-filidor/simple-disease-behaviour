function my_trajectory(x0, I, r, p, c)
    if cond(x0, I)
        X0 = [x0 I]
        pars = [r p c]
        dynamics = ContinuousDynamicalSystem(dynamics_rule_si, X0, pars; diffeq)
        SOL, t = trajectory(dynamics, T; Ttr=Ttr, Δt=sampling_time)
        t0 = Int((Ttr + 1) / sampling_time)
        x = SOL[t0:end, 1]
        I = SOL[t0:end, 2]
        t = t[t0:end]
        println(size(I))
        println(size(t))
        return [x, I]
    end
end

function my_trajectory_no_condition(x0, I, r, p, c)
    X0 = [x0 I]
    pars = [r p c]
    dynamics = ContinuousDynamicalSystem(dynamics_rule_si, X0, pars; diffeq)
    SOL, t = trajectory(dynamics, T; Δt=sampling_time, Ttr=Ttr)
    x = SOL[:, 1]
    I = SOL[:, 2]
    println(size(I))
    return [x, I]
end
