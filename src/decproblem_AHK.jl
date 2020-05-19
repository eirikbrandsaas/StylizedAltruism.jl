
function SolveVk2!(M)
  for (ixk,xk) in enumerate(M.np.x_grd)
    for (iyk,yk) in enumerate(M.np.y_grd)
      M.gk.c[M.np.na][ixk,iyk] = xk + yk
      M.gk.x′[M.np.na][ixk,iyk] = 0.0
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

function SolveAHK!(M)
  SolveVk2!(M)
  SolveVp2!(M)
end
