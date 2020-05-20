
function SolveVk2!(M)
  for (ixk,xk) in enumerate(M.np.x_grd)
    for (iyk,yk) in enumerate(M.np.y_grd)
      ck = xk + yk
      M.gk.c[M.np.na][ixk,iyk] = ck
      M.gk.x′[M.np.na][ixk,iyk] = 0.0
      M.Vk[M.np.na][ixk,iyk] = util(ck,M.mp.γ)
    end
  end
end

function SolveVp2!(M)
  # gc_itp = []
  np = M.np
  gck_itp = [LinearInterpolation((np.x_grd,),M.gk.c[M.np.na][:,iyk],extrapolation_bc=Line()) for iyk = 1:np.ny]
  vtmp = fill(-Inf64,np.ntc)
  ia = M.np.na
  for (ixp,xp) in enumerate(M.np.x_grd)
    for (ixk,xk) in enumerate(M.np.x_grd)
      for (iyk,yk) in enumerate(M.np.y_grd)
        ## Loop over possible choices
        vtmp .= -Inf64
        for (itp,tp) in enumerate(M.np.tc_grd)
          cp = xp - tp
          if cp > 0.0
            ck = gck_itp[iyk](xk + tp)
            vtmp[itp] = util(cp,M.mp.γ) + util(ck,M.mp.γ)
          end
        end
        imax = argmax(vtmp)
        tp = M.np.tc_grd[imax]
        M.gp.t[ia][ixp,ixk,iyk] = tp
        M.gp.c[ia][ixp,ixk,iyk] = xp - tp
        M.Vp[ia][ixp,ixk,iyk] = vtmp[imax]
      end
    end
  end
end

function SolveVk1!(M)
  np = M.np
  mp = M.mp
  vtmp = fill(-Inf64,np.ntc)
  ia = M.np.na-1

  gtp_itp = [LinearInterpolation((np.x_grd,),M.gp.t[ia+1][ixpn,:,iyk],extrapolation_bc=Line()) for ixpn = 1:np.nx, iyk = 1:np.ny]
  Vk_itp =  [LinearInterpolation((np.x_grd,),M.Vk[ia+1][:,iyk],extrapolation_bc=Line()) for ixpn = 1:np.nx, iyk = 1:np.ny]

  for (ixpn,xp) in enumerate(M.np.x_grd)
    for (ixk,xk) in enumerate(M.np.x_grd)
      for (iyk,yk) in enumerate(M.np.y_grd)
        ## Loop over possible choices
        vtmp .= -Inf64
        for (ixkn,xkn) in enumerate(M.np.xc_grd)
          ck = xk + yk - xkn/(1.0 + mp.rf)
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
        M.gk.c[ia][ixpn,ixk,iyk] = xk + yk - xkn/(1.0 + mp.rf)
        M.Vk[ia][ixpn,ixk,iyk] = vtmp[imax]
      end
    end
  end
end



function SolveAHK!(M)
  SolveVk2!(M)
  SolveVp2!(M)
  SolveVk1!(M)
end
