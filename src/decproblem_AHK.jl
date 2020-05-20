
function SolveVk2!(M)
  for (ixk,xk) in enumerate(M.np.x_grd)
    for (iyk,yk) in enumerate(M.np.y_grd)
      xkn = 0.0 # Don't save ion the last period
      ck = ck_bc(xk,yk,xkn,M.mp.rf)
      M.gk.c[M.np.na][ixk,iyk] = ck
      M.gk.x′[M.np.na][ixk,iyk] = xkn
      M.Vk[M.np.na][ixk,iyk] = util(ck,M.mp.γ)
    end
  end
end

function SolveVp2!(M)
  # gc_itp = []
  np = M.np
  mp = M.mp
  gck_itp = [LinearInterpolation((np.x_grd,),M.gk.c[M.np.na][:,iyk],extrapolation_bc=Line()) for iyk = 1:np.ny]
  vtmp = fill(-Inf64,np.ntc)
  ia = M.np.na
  for (ixp,xp) in enumerate(M.np.x_grd)
    for (ixk,xk) in enumerate(M.np.x_grd)
      for (iyk,yk) in enumerate(M.np.y_grd)
        ## Loop over possible choices
        vtmp .= -Inf64
        xpn = 0.0
        for (itp,tp) in enumerate(M.np.tc_grd)
          cp = cp_bc(xp,xpn,tp,mp.rf)
          if cp > 0.0
            ck = gck_itp[iyk](xk + tp)
            vtmp[itp] = util(cp,mp.γ) + mp.η*util(ck,mp.γ)
          end
        end
        imax = argmax(vtmp)
        tp = M.np.tc_grd[imax]
        M.gp.t[ia][ixp,ixk,iyk] = tp
        M.gp.c[ia][ixp,ixk,iyk] = cp_bc(xp,xpn,tp,mp.rf)
        M.Vp[ia][ixp,ixk,iyk] = vtmp[imax]
      end
    end
  end
end

function SolveVk1!(M)
  np = M.np
  mp = M.mp
  vtmp = fill(-Inf64,np.nxc)
  ia = M.np.na-1

  gtp_itp = [LinearInterpolation((np.x_grd,),M.gp.t[ia+1][ixpn,:,iyk],extrapolation_bc=Line()) for ixpn = 1:np.nx, iyk = 1:np.ny]
  Vk_itp =  [LinearInterpolation((np.x_grd,),M.Vk[ia+1][:,iyk],extrapolation_bc=Line()) for ixpn = 1:np.nx, iyk = 1:np.ny]

  for (ixpn,xp) in enumerate(M.np.x_grd)
    for (ixk,xk) in enumerate(M.np.x_grd)
      for (iyk,yk) in enumerate(M.np.y_grd)
        ## Loop over possible choices
        vtmp .= -Inf64
        for (ixkn,xkn) in enumerate(M.np.xc_grd)
          ck = ck_bc(xk,yk,xkn,mp.rf)
          if ck > 0.0 && xkn >= 0.0
            vtmp[ixkn] = util(ck,mp.γ)
            for iykn = 1:np.ny
              tpn = gtp_itp[ixpn,iykn](xkn)
              vtmp[ixkn] += mp.β*np.Πy[iykn,iyk]*Vk_itp[iykn](xkn + tpn)
            end
          end
        end
        imax = argmax(vtmp)
        xkn = M.np.xc_grd[imax]
        M.gk.x′[ia][ixpn,ixk,iyk] = xkn
        M.gk.c[ia][ixpn,ixk,iyk] = ck_bc(xk,yk,xkn,mp.rf)
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
  gxkn_itp = [LinearInterpolation((np.x_grd,np.x_grd),M.gk.x′[ia][:,:,iyk],extrapolation_bc=Line()) for iyk = 1:np.ny]
  Vp_itp = [LinearInterpolation((np.x_grd,np.x_grd),M.Vp[ia+1][:,:,iyk],extrapolation_bc=Line()) for iyk = 1:np.ny]

  for (ixp,xp) in enumerate(M.np.x_grd)
    for (ixk,xk) in enumerate(M.np.x_grd)
      for (iyk,yk) in enumerate(M.np.y_grd)
        ## Loop over possible choices
        vtmp .= -Inf64
        for (ixpn,xpn) in enumerate(M.np.xc_grd)
          for (itp,tp) in enumerate(M.np.tc_grd)
          cp = cp_bc(xp,xpn,tp,mp.rf) #xp - xpn/(1.0 + mp.rf) - tp
          xkn = gxkn_itp[iyk](xpn,xk+tp)
          ck = xk + yk - xkn/(1.0 + mp.rf)

            if cp > 0.0 && xpn >= 0.0
              vtmp[ixpn,itp] = util(cp,mp.γ) +  mp.η*util(ck,mp.γ)
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
end
