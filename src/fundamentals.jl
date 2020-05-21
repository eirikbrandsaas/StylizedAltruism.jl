function utilk(c::Real,h::Real,γ::Real,ξ::Real)
    util = ((c^(1.0-ξ)*h^(ξ))^(1.0-γ) - 1.0)/(1.0-γ)
end

function utilp(c::Real,γ::Real)
    util = (c^(1.0-γ) - 1.0)/(1.0-γ)
end
function xkn(xk,yk,ck,rf)
    xkn = (xk + yk - ck)*(1.0 + rf)
end

function xpn(xp,cp,tp,rf)
    xpn = (xp - cp - tp)*(1.0+rf)
end

function ck_bc(xk,wk,xkn,rf)
    ck = xk + wk - xkn/(1.0 + rf)
end

function cp_bc(xp,xpn,tp,rf)
    cp = xp - tp - xpn/(1.0+rf)
end

function BorrConstr()
    0.0 # Benchmark model has no-borrowing
end

function income(yk::Real,ia::Int,np::NumPar)
    wage = yk*np.inc_grd[ia]
end

function adj(h1::Real,h2::Real,κ::Real)
    if h1 != h2
        cost = κ*h1
    else
        cost = 0.0
    end

    return cost
end
