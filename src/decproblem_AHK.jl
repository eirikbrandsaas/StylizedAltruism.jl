
function SolveVk2!(M)
  np = M.np
  mp = M.mp
  ia = np.na
  for (ixk,xk) in enumerate(np.x_grd)
    for (iyk,yk) in enumerate(np.y_grd)
      for (ihk,hk) in enumerate(np.h_grd)
        wk = income(yk,ia,np)
        xkn = 0.0 # Don't save ion the last period
        ck = ck_bc(xk,wk,xkn,mp.rf)
        M.gk.c[ia][ixk,iyk,ihk] = ck
        M.gk.x′[ia][ixk,iyk,ihk] = xkn
        M.Vk[ia][ixk,iyk,ihk] = utilk(ck,0.0,mp.γ,mp.ξ)
      end
    end
  end
end

function SolveVp2!(M)
  np = M.np
  mp = M.mp
  gck_itp::Array{Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}},Gridded{Linear},Line{Nothing}},2} =
    [LinearInterpolation((np.x_grd,),M.gk.c[M.np.na][:,iyk,ihk],extrapolation_bc=Line()) for iyk = 1:np.ny,ihk = 1:np.nh]
  vtmp = fill(-Inf64,np.ntc)
  ia = M.np.na
  for (ixp,xp) in enumerate(np.x_grd)
    for (ixk,xk) in enumerate(np.x_grd)
      for (iyk,yk) in enumerate(np.y_grd)
        for (ihk,hk) in enumerate(np.h_grd)
        ## Loop over possible choices
        vtmp .= -Inf64
        xpn = 0.0
        for (itp,tp) in enumerate(M.np.tc_grd)
          cp = cp_bc(xp,xpn,tp,mp.rf)
          if cp > 0.0 && xpn >= BorrConstr()
            ck = gck_itp[iyk](xk + tp)
            vtmp[itp] = utilp(cp,mp.γ) + mp.η*utilk(ck,0.0,mp.γ,mp.ξ)
          end
        end
        imax = argmax(vtmp)
        tp = M.np.tc_grd[imax]
        M.gp.t[ia][ixp,ixk,iyk,ihk] = tp
        M.gp.c[ia][ixp,ixk,iyk,ihk] = cp_bc(xp,xpn,tp,mp.rf)
        M.Vp[ia][ixp,ixk,iyk,ihk] = vtmp[imax]
        end
      end
    end
  end
end

function SolveVk1!(M)
  np = M.np
  mp = M.mp
  vtmp = fill(-Inf64,np.nxc)
  ia = M.np.na-1

  gtp_itp::Array{Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}},Gridded{Linear},Line{Nothing}},2} =
    [LinearInterpolation((np.x_grd,),M.gp.t[ia+1][ixpn,:,iyk],extrapolation_bc=Line()) for ixpn = 1:np.nx, iyk = 1:np.ny]
  Vk_itp::Array{Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}},Gridded{Linear},Line{Nothing}},2} =
    [LinearInterpolation((np.x_grd,),M.Vk[ia+1][:,iyk],extrapolation_bc=Line()) for ixpn = 1:np.nx, iyk = 1:np.ny]

  for (ixpn,xp) in enumerate(np.x_grd)
    for (ixk,xk) in enumerate(np.x_grd)
      for (iyk,yk) in enumerate(np.y_grd)
        wk = income(yk,ia,np)
        ## Loop over possible choices
        vtmp .= -Inf64
        for (ixkn,xkn) in enumerate(M.np.xc_grd)
          ck = ck_bc(xk,wk,xkn,mp.rf)
          if ck > 0.0 && xkn >= BorrConstr()
            vtmp[ixkn] = utilk(ck,0.0,mp.γ,mp.ξ)
            for iykn = 1:np.ny
              tpn = gtp_itp[ixpn,iykn](xkn)
              vtmp[ixkn] += mp.β*np.Πy[iykn,iyk]*Vk_itp[iykn](xkn + tpn)
            end
          end
        end
        imax = argmax(vtmp)
        xkn = M.np.xc_grd[imax]
        M.gk.x′[ia][ixpn,ixk,iyk] = xkn
        M.gk.c[ia][ixpn,ixk,iyk] = ck_bc(xk,wk,xkn,mp.rf)
        M.Vk[ia][ixpn,ixk,iyk] = vtmp[imax]
      end
    end
  end
end

function SolveVp1!(M)
  np = M.np
  mp = M.mp
  vtmp = fill(-Inf64,np.nxc,np.ntc)
  ia = M.np.na-1

  # gck_itp = [LinearInterpolation((np.x_grd,np.x_grd),M.gk.c[ia][:,:,iyk],extrapolation_bc=Line()) for iyk = 1:np.ny]
  gxkn_itp::Array{Interpolations.Extrapolation{Float64,2,Interpolations.GriddedInterpolation{Float64,2,Float64,Gridded{Linear},Tuple{Array{Float64,1},Array{Float64,1}}},Gridded{Linear},Line{Nothing}},1} =
     [LinearInterpolation((np.x_grd,np.x_grd),M.gk.x′[ia][:,:,iyk],extrapolation_bc=Line()) for iyk = 1:np.ny]
  Vp_itp::Array{Interpolations.Extrapolation{Float64,2,Interpolations.GriddedInterpolation{Float64,2,Float64,Gridded{Linear},Tuple{Array{Float64,1},Array{Float64,1}}},Gridded{Linear},Line{Nothing}},1} =
     [LinearInterpolation((np.x_grd,np.x_grd),M.Vp[ia+1][:,:,iyk],extrapolation_bc=Line()) for iyk = 1:np.ny]

  for (ixp,xp) in enumerate(M.np.x_grd)
    for (ixk,xk) in enumerate(M.np.x_grd)
      for (iyk,yk) in enumerate(M.np.y_grd)
        wk = income(yk,ia,np)
        vtmp .= -Inf64

        ## Loop over possible choices
        for (ixpn,xpn) in enumerate(M.np.xc_grd)
          for (itp,tp) in enumerate(M.np.tc_grd)
          cp = cp_bc(xp,xpn,tp,mp.rf) #xp - xpn/(1.0 + mp.rf) - tp
          xkn = gxkn_itp[iyk](xpn,xk+tp)
          ck = ck_bc(xk+tp,wk,xkn,mp.rf)
            if cp > 0.0 && xpn >= BorrConstr()
              vtmp[ixpn,itp] = utilp(cp,mp.γ) +  mp.η*utilk(ck,0.0,mp.γ,mp.ξ)
              for iykn = 1:np.ny
                vtmp[ixpn,itp] += mp.β*np.Πy[iykn,iyk]*Vp_itp[iykn](xpn,xkn)
              end
            end
          end
        end
        imax = argmax(vtmp)
        xpn = M.np.xc_grd[imax[1]]
        tp = M.np.tc_grd[imax[2]]
        M.gp.x′[ia][ixp,ixk,iyk] = xpn
        M.gp.c[ia][ixp,ixk,iyk] = cp_bc(xp,xpn,tp,mp.rf)
        M.gp.t[ia][ixp,ixk,iyk] = tp
        M.Vp[ia][ixp,ixk,iyk] = vtmp[imax]
      end
    end
  end
end


function SolveAHK!(M)
  SolveVk2!(M)
  SolveVp2!(M)
  SolveVk1!(M)
  SolveVp1!(M)

  TestNumericalSolution(M) # Function that tests for numerical issues
end
