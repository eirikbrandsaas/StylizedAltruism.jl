function util(c::Real,γ::Real)
    util = (c^(1.0-γ) - 1.0)/(1.0-γ)
end

function xkn(xk,yk,ck,rf)
    xkn = (xk + yk - ck)*(1.0 + rf)
end

function xpn(xp,cp,tp,rf)
    xpn = (xp - cp - tp)*(1.0+rf)
end

function ck_bc(xk,yk,xkn,rf)
    ck = xk + yk - xkn/(1.0 + rf)
end

function cp_bc(xp,xpn,tp,rf)
    cp = xp - tp - xpn/(1.0+rf)
end

function BorrConstr()
    0.0 # Benchmark model has no-borrowing
end
