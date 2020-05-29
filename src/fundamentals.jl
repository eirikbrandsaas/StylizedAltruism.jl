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

function ck_bc(xk,wk,xkn,hk,hkn,rf,κ)
    ck = xk + wk - xkn/(1.0 + rf) - adj(hk,hkn,κ) - hkn
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

function gke_opt(np::StylizedAltruism.NumPar,gk::StylizedAltruism.Polk)
  gke = Polk_eq(np)
  gke.h = gk.h
  for ia in 1:np.na
    for ixk in 1:np.nx, iy in 1:np.ny
      if ia < np.na
        for ihki = 1:np.nhi, ixp in 1:np.nx
          gke.c[ia][ixk,ixp,iy,ihki] = gk.c[ia][ixk,ixp,iy,ihki,gk.h[ia][ixk,ixp,iy]]
          gke.x′[ia][ixk,ixp,iy,ihki] = gk.x′[ia][ixk,ixp,iy,ihki,gk.h[ia][ixk,ixp,iy]]
        end
      else
        for ihk = 1:np.nh
          gke.c[ia][ixk,iy,ihk] = gk.c[ia][ixk,iy,ihk,gk.h[ia][ixk,iy,ihk]]
          gke.x′[ia][ixk,iy,ihk] = gk.x′[ia][ixk,iy,ihk,gk.h[ia][ixk,iy,ihk]]
        end
      end
    end
  end

  return gke
end
