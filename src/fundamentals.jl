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

# function ck_bc(xk,wk,xkn,rf)
#     ck = xk + wk - xkn/(1.0 + rf)
# end

function ck_bc(xk,wk,xkn,hk,hkn,rf,κ,p,own::Bool)
    if own == false
        ck = xk + wk - xkn/(1.0 + rf) - adj(hk,hkn,κ) - hkn
    else
        ck = xk + wk - xkn/(1.0 + rf) - adj(hk,hkn,κ) - p*hkn
    end
end

function ck_bc2(xk,wk,xkn,hk,hkn,rf,κ,os,rs,s,ok::Bool)
    if ok == false
        ck = xk + wk - xkn/(1.0 + rf) - adj2(hk,hkn,κ,ok) - (1.0 + rs*s)*hkn
    else
        ck = xk + wk - xkn/(1.0 + rf) - adj2(hk,hkn,κ,ok) + (os*s)*hk
    end
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

function adj2(h1::Real,h2::Real,κ::Real,ok::Bool)
    if h1 != h2
        cost = κ*h1 + Inf64*ok
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
          ihkn = gk.h[ia][ixk,ixp,iy]
          iokn = gk.o[ia][ixk,ixp,iy]
          gke.c[ia][ixk,ixp,iy,ihki] = gk.c[ia][ixk,ixp,iy,ihki,ihkn,iokn]
          gke.x′[ia][ixk,ixp,iy,ihki] = gk.x′[ia][ixk,ixp,iy,ihki,ihkn,iokn]
        end
      else
        for ihk = 1:np.nh, io = 1:np.no, is = 1:np.ns
          gke.c[ia][ixk,iy,ihk,io,is] = gk.c[ia][ixk,iy,ihk,io,is,gk.h[ia][ixk,iy,ihk,io,is]]
          gke.x′[ia][ixk,iy,ihk,io,is] = gk.x′[ia][ixk,iy,ihk,io,is,gk.h[ia][ixk,iy,ihk,io,is]]
        end
      end
    end
  end

  return gke
end
