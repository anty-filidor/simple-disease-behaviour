function cond(x, I)
    if x >= 0 && x <= 1 && I >= 0 && I <= 1 && I > 0
        return true
    else
        return false
    end
end

function dynamics_rule_si(x0, param, t)
    R, p, c = param
    x, I = x0
    dxdt = -p*(1+I)*x*(1-x)^2 + p*(1-I)*(1-x)*x^2 + c*(1-I)*(1-x)
    dIdt = R*x*I*(1-I) - I
    r = SVector{2}([dxdt; dIdt])
    return r
end
