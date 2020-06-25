
function SolveVk2!(M)
  np = M.np
  mp = M.mp
  ia = np.na
  vtmp = fill(-Inf64,np.nh)
  for (ixk,xk) in enumerate(np.x_grd)
    for (iyk,yk) in enumerate(np.y_grd)
      for (ihk,hk) in enumerate(np.h_grd)
        wk = income(yk,ia,np)
        xkn = 0.0 # Don't save in the last period
        vtmp .= -Inf64
        for (ihkn,hkn) in enumerate(np.h_grd)
          ck = ck_bc(xk,wk,xkn,hk,hkn,mp.rf,mp.κ)

          M.gk.c[ia][ixk,iyk,ihk,ihkn] = ck
          M.gk.x′[ia][ixk,iyk,ihk,ihkn] = xkn
          if ck > 0.0 # If the transfer pushes the kid into infeasible choices!
            vtmp[ihkn] = utilk(ck,hkn,mp.γ,mp.ξ)
          end
        end
        imax = argmax(vtmp)
        M.gk.disc[ia][ixk,iyk,ihk,imax] = 1
        M.gk.h[ia][ixk,iyk,ihk] = imax
        M.Vk[ia][ixk,iyk,ihk] = vtmp[imax]
      end
    end
  end
end

function SolveVp2!(M)
  np = M.np
  mp = M.mp
  gck_itp::Array{Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}},Gridded{Linear},Line{Nothing}},3} =
    [LinearInterpolation((np.x_grd,),M.gk.c[M.np.na][:,iyk,ihk,ihkn],extrapolation_bc=Line()) for iyk = 1:np.ny,ihk = 1:np.nh, ihkn = 1:np.nh]
  gdisck_itp::Array{Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}},Gridded{Linear},Flat{Nothing}},3} =
    [LinearInterpolation((np.x_grd,),Float64.(M.gk.disc[M.np.na][:,iyk,ihk,ihkn]),extrapolation_bc=Flat()) for iyk = 1:np.ny,ihk = 1:np.nh, ihkn = 1:np.nh]
  vtmp = fill(-Inf64,np.ntc)
  ia = M.np.na
  for (ixp,xp) in enumerate(np.x_grd)
    for (ixk,xk) in enumerate(np.x_grd)
      for (iyk,yk) in enumerate(np.y_grd)
        for (ihk,hk) in enumerate(np.h_grd)
        ## Loop over possible choices
        vtmp .= -Inf64
        xpn = 0.0
        for (itp,tp) in enumerate(np.tc_grd)
        # Threads.@threads for (itp,tp) in collect(enumerate(np.tc_grd)) # Works, but is slower!?
          vtmp[itp] = evalVp2(xp,xk,iyk,ihk,xpn,tp,gdisck_itp,gck_itp,mp,np) # ,mp,np
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

function evalVp2(xp::Real,xk::Real,iyk::Int,ihk::Int,xpn::Real,tp::Real,
    gdisck_itp::Array{Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}},Gridded{Linear},Flat{Nothing}},3},
    gck_itp::Array{Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}},Gridded{Linear},Line{Nothing}},3},
    mp::ModPar,np::NumPar) # ,,np::NumPar

  cp = cp_bc(xp,xpn,tp,mp.rf)
  if cp > 0.0
    val = utilp(cp,mp.γ)
    for (ihkn,hkn) in enumerate(np.h_grd) # We loop over all the possible discrete housing choices..
      prob = gdisck_itp[iyk,ihk,ihkn](xk + tp) # And then find the probability that this choice is made!
      @assert prob >= 0.0
      @assert prob <= 1.0
      if xpn >= BorrConstr() && prob > 0
        ck = gck_itp[iyk,ihk,ihkn](xk + tp)
        if ck > 0.0
          val += prob*(mp.η*utilk(ck,hkn,mp.γ,mp.ξ))
        end
      end
    end
  else
    val = -Inf64
  end

  return val

end

function SolveVk1!(M)
  np = M.np
  mp = M.mp
  vtmp = fill(-Inf64,np.nxc,np.nh)
  ia = M.np.na-1

  gtp_itp::Array{Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}},Gridded{Linear},Line{Nothing}},3} =
    [LinearInterpolation((np.x_grd,),M.gp.t[ia+1][ixpn,:,iyk,ihkn],extrapolation_bc=Line()) for ixpn = 1:np.nx, iyk = 1:np.ny, ihkn = 1:np.nh]
  Vk_itp::Array{Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}},Gridded{Linear},Line{Nothing}},2} =
    [LinearInterpolation((np.x_grd,),M.Vk[ia+1][:,iyk,ihkn],extrapolation_bc=Line()) for iyk = 1:np.ny, ihkn = 1:np.nh]

  for (ixpn,xp) in enumerate(np.x_grd)
    for (ixk,xk) in enumerate(np.x_grd)
      for (iyk,yk) in enumerate(np.y_grd)
        for (ihki,hki) in enumerate(np.hi_grd)
        wk = income(yk,ia,np)
        ## Loop over possible choices
        vtmp .= -Inf64
        for (ihkn,hkn) in enumerate(np.h_grd)
          for (ixkn,xkn) in enumerate(np.xc_grd)
            hkn = np.h_grd[ihkn]
            ck = ck_bc(xk,wk,xkn,hki,hkn,mp.rf,mp.κ)
            if ck > 0.0 && xkn >= BorrConstr()
              vtmp[ixkn,ihkn] = utilk(ck,hkn,mp.γ,mp.ξ)
              for iykn = 1:np.ny
                tpn = gtp_itp[ixpn,iykn,ihkn](xkn)
                vtmp[ixkn,ihkn] += mp.β*np.Πy[iyk,iykn]*Vk_itp[iykn,ihkn](xkn + tpn)
              end
            end
          end
          imax = argmax(vtmp[:,ihkn])
          xkn = M.np.xc_grd[imax]
          M.gk.x′[ia][ixpn,ixk,iyk,ihki,ihkn] = xkn
          M.gk.c[ia][ixpn,ixk,iyk,ihki,ihkn] = ck_bc(xk,wk,xkn,0.0,hkn,mp.rf,mp.κ)
        end
        imax = argmax(vtmp) # Only need the optimal housing choice (second dimensio)
        ihmax = imax[2]
        M.gk.disc[ia][ixpn,ixk,iyk,ihki,ihmax] = 1
        M.gk.h[ia][ixpn,ixk,iyk,ihki] = ihmax
        M.Vk[ia][ixpn,ixk,iyk,ihki] = vtmp[imax]
      end
      end
    end
  end
end

##
function SolveVp1!(M)
  np = M.np
  mp = M.mp
  vtmp = fill(-Inf64,np.nxc,np.ntc)
  ia = M.np.na-1

  # gck_itp = [LinearInterpolation((np.x_grd,np.x_grd),M.gk.c[ia][:,:,iyk],extrapolation_bc=Line()) for iyk = 1:np.ny]
  gxkn_itp::Array{Interpolations.Extrapolation{Float64,2,Interpolations.GriddedInterpolation{Float64,2,Float64,Interpolations.Gridded{Interpolations.Linear},Tuple{Array{Float64,1},Array{Float64,1}}},Interpolations.Gridded{Interpolations.Linear},Interpolations.Line{Nothing}},3} =
     [LinearInterpolation((np.x_grd,np.x_grd),M.gk.x′[ia][:,:,iyk,ihki,ihkn],extrapolation_bc=Line()) for iyk = 1:np.ny, ihki = 1:np.nhi, ihkn = 1:np.nh]
  gdisck_itp::Array{Interpolations.Extrapolation{Float64,2,Interpolations.GriddedInterpolation{Float64,2,Float64,Interpolations.Gridded{Interpolations.Linear},Tuple{Array{Float64,1},Array{Float64,1}}},Interpolations.Gridded{Interpolations.Linear},Interpolations.Flat{Nothing}},3} =
     [LinearInterpolation((np.x_grd,np.x_grd),Float64.(M.gk.disc[ia][:,:,iyk,ihki,ihkn]),extrapolation_bc=Flat()) for iyk = 1:np.ny, ihki = 1:np.nhi, ihkn = 1:np.nh]
  Vp_itp::Array{Interpolations.Extrapolation{Float64,2,Interpolations.GriddedInterpolation{Float64,2,Float64,Gridded{Linear},Tuple{Array{Float64,1},Array{Float64,1}}},Gridded{Linear},Line{Nothing}},2} =
     [LinearInterpolation((np.x_grd,np.x_grd),M.Vp[ia+1][:,:,iyk,ihkn],extrapolation_bc=Line()) for iyk = 1:np.ny, ihkn = 1:np.nh]

  for (ixp,xp) in enumerate(M.np.x_grd)
    for (ixk,xk) in enumerate(M.np.x_grd)
      for (iyk,yk) in enumerate(M.np.y_grd)
        for (ihki,hki) in enumerate(np.hi_grd)
        wk = income(yk,ia,np)
        vtmp .= -Inf64

        ## Loop over possible choices
        Threads.@threads for (ixpn,xpn) in collect(enumerate(M.np.xc_grd))
           for (itp,tp) in enumerate(M.np.tc_grd)
            vtmp[ixpn,itp] = evalVp1(ixp,ixk,iyk,ihki,hki,wk,xk,xp,xpn,tp,gdisck_itp,gxkn_itp,Vp_itp,mp,np)
          end #tp
        end #xp
        imax = argmax(vtmp)
        xpn = M.np.xc_grd[imax[1]]
        tp = M.np.tc_grd[imax[2]]
        M.gp.x′[ia][ixp,ixk,iyk,ihki] = xpn
        M.gp.c[ia][ixp,ixk,iyk,ihki] = cp_bc(xp,xpn,tp,mp.rf)
        M.gp.t[ia][ixp,ixk,iyk,ihki] = tp
        M.Vp[ia][ixp,ixk,iyk,ihki] = vtmp[imax]
      end
      end
    end
  end
end

function evalVp1(ixp::Int,ixk::Int,iyk::Int,ihki::Int,hki::Real,wk::Real,xk::Real,xp::Real,xpn::Real,tp::Real,gdisck_itp,gxkn_itp,Vp_itp,mp::ModPar,np::NumPar)

  cp = cp_bc(xp,xpn,tp,mp.rf) #xp - xpn/(1.0 + mp.rf) - tp
  if cp > 0.0 && xpn >= BorrConstr()
    val  = utilp(cp,mp.γ)
    for (ihkn,hkn) in enumerate(np.h_grd)
      prob = gdisck_itp[iyk,ihki,ihkn](xpn,xk+tp)
      @assert prob >= 0.0
      @assert prob <= 1.0
      if prob > 0.0
        xkn = gxkn_itp[iyk,ihki,ihkn](xpn,xk+tp)
        ck = ck_bc(xk+tp,wk,xkn,hki,hkn,mp.rf,mp.κ)
        if ck > 0.0
          val += prob*(mp.η*utilk(ck,hkn,mp.γ,mp.ξ))
          for iykn = 1:np.ny
            val += prob*mp.β*np.Πy[iyk,iykn]*Vp_itp[iykn,ihkn](xpn,xkn)
          end
        end
      end # prob
    end
  else
    val = -Inf64
  end


  return val
end

function SolveAHK!(M)
  SolveVk2!(M)
  # @assert minimum(M.Vk[2])  > -Inf64
  SolveVp2!(M)
  # @assert minimum(M.Vp[2])  > -Inf64
  SolveVk1!(M)
  # @assert minimum(M.Vk[1])  > -Inf64
  SolveVp1!(M)
  # @assert minimum(M.Vp[1])  > -Inf64

  M.gke = gke_opt(M.np,M.gk)

  TestNumericalSolution(M) # Function that tests for numerical issues
end
