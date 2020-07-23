function SolveVc2(np,mp)
    cp_pol = fill(-1.0,np.nxf,np.ny,np.nh,np.no,np.ns,np.nθ)
    ck_pol = fill(-1.0,np.nxf,np.ny,np.nh,np.no,np.ns,np.nθ)
    xf_pol = fill(-1.0,np.nxf,np.ny,np.nh,np.no,np.ns,np.nθ)
    ihkn_pol = fill(-1,np.nxf,np.ny,np.nh,np.no,np.ns,np.nθ)
    val = fill(-1.0,np.nxf,np.ny,np.nh,np.no,np.ns,np.nθ)
    vtmp = fill(-Inf64,np.nh,np.ntc)
    ia = np.na
    for (ixf,xf) in enumerate(np.xf_grd)
      for (iyk,yk) in enumerate(np.y_grd)
        for (ihk,hk) in enumerate(np.h_grd)
          for (io,ok) in enumerate(np.o_grd)
            for (is,sk) in enumerate(np.s_grd)
              for (iθ,θ) in enumerate(np.θ_grd)
                vtmp .= -1e64
                xpn = 0.0
                wk = StylizedAltruism.income(yk,ia,np)
                for (ihkn,hkn) in enumerate(np.h_grd)
                  for (ick,ck) in enumerate(np.xc_grd)
                    cp = StylizedAltruism.ck_bc2(xf-ck,wk,xpn,hk,hkn,mp.rf,mp.κ,mp.os,mp.rs,sk,ok)
                    if ck > 0.0 && cp >0.0
                      vtmp[ihkn,ick] =  (1.0 - θ)*StylizedAltruism.utilp(cp,mp.γ) + (θ + (1.0 - θ)*mp.η)*StylizedAltruism.utilk(ck,hkn,mp.γ,mp.ξ)
                    end
                  end
                end
                imax = argmax(vtmp)
                ihkn = imax[1]
                ick = imax[2]
                hkn = np.h_grd[ihkn]
                ck = np.xc_grd[ick]
                cp_pol[ixf,iyk,ihk,io,is,iθ] = StylizedAltruism.ck_bc2(xf-ck,wk,xpn,hk,hkn,mp.rf,mp.κ,mp.os,mp.rs,sk,ok)
                ck_pol[ixf,iyk,ihk,io,is,iθ] = ck
                xf_pol[ixf,iyk,ihk,io,is,iθ] = xpn
                ihkn_pol[ixf,iyk,ihk,io,is,iθ] = ihkn
                val[ixf,iyk,ihk,io,is,iθ] = vtmp[imax]
              end
            end
          end
        end
      end
    end

  return cp_pol,ck_pol,ihkn_pol,xf_pol,val
end

function SolveVc1(np,mp,valnxt)
  cp_pol = fill(-2.0,np.nxf,np.ny,np.nhi,np.nθ)
  ck_pol = fill(-1.0,np.nxf,np.ny,np.nhi,np.nθ)
  xp′_pol = fill(-1.0,np.nxf,np.ny,np.nhi,np.nθ)
  ihkn_pol = fill(-1,np.nxf,np.ny,np.nhi,np.nθ)
  val = fill(-1.0,np.nxf,np.ny,np.nhi,np.nθ)
  vtmp = fill(-Inf64,np.nh,np.ntc,np.nxc)
  ia = np.na

  @assert np.no == 1 # Only works with one ownership status

  vnxt = [LinearInterpolation((np.xf_grd,),valnxt[:,iyk,ihkn,io,is,iθ],extrapolation_bc=Line()) for iyk = 1:np.ny, ihkn = 1:np.nh, io = 1:np.no, is=1:np.ns, iθ = 1:np.nθ]

  for (ixf,xf) in enumerate(np.xf_grd)
    for (iyk,yk) in enumerate(np.y_grd)
      for (ihk,hk) in enumerate(np.hi_grd)
        for (iθ,θ) in enumerate(np.θ_grd)
          ## Loop over possible choices
          vtmp .= -1e64
          wk = StylizedAltruism.income(yk,ia,np)
          #=Threads.@threads=# for (ihkn,hkn) in collect(enumerate(np.h_grd))
            for (ick,ck) in enumerate(np.xc_grd)
              for (ixpn,xpn) in enumerate(np.xc_grd)
                 vtmp[ihkn,ick,ixpn] = Vc1_eval(θ,xf+wk,ck,hkn,xpn,iyk,ihkn,mp,np,vnxt[:,:,:,:,iθ])
              end
            end
          end
          imax = argmax(vtmp)
          ihkn = imax[1]
          ick = imax[2]

          ixpn = imax[3]
          hkn = np.h_grd[ihkn]
          ck = np.xc_grd[ick]
          xpn = np.xc_grd[ixpn]
          ck_pol[ixf,iyk,ihk,iθ] = ck
          ihkn_pol[ixf,iyk,ihk,iθ] = ihkn
          cp_pol[ixf,iyk,ihk,iθ] = StylizedAltruism.cp_bc(xf+wk-ck-hkn,xpn,0.0,mp.rf)
          xp′_pol[ixf,iyk,ihk,iθ] = xpn
          val[ixf,iyk,ihk,iθ] = vtmp[imax]
        end
      end
    end
  end


  return cp_pol, ck_pol,ihkn_pol, xp′_pol, val
end


function Vc1_eval(θ,xf,ck,hkn,xpn,iyk,ihkn,mp,np,vnxt)
  cp = StylizedAltruism.cp_bc(xf-ck-hkn,xpn,0.0,mp.rf)
  iokn = 1

  vtmp = 0.0
  if ck > 0.0 && cp >0.0 && xpn > 0.0
    vtmp =  (1.0 - θ)*StylizedAltruism.utilp(cp,mp.γ) + (θ + (1.0 - θ)*mp.η)*StylizedAltruism.utilk(ck,hkn,mp.γ,mp.ξ)
    for iykn = 1:np.ny
      for is = 1:np.ns
          vtmp += mp.β*np.Πy[iyk,iykn]*np.Πs[is]*vnxt[iykn,ihkn,iokn,is](xpn)
      end
    end
  else
    vtmp = -1e60
  end

  return   vtmp
end
#


function SolveCommitment(np,mp,switch)
  Mc = ModelCommitment(np,mp,switch)
  Mc.cp[2], Mc.ck[2], Mc.ihkn[2], Mc.xf′[2], Mc.Vd[2] = SolveVc2(np,mp)
  Mc.cp[1], Mc.ck[1], Mc.ihkn[1], Mc.xf′[1], Mc.Vd[1] = SolveVc1(np,mp,Mc.Vd[2])

  return Mc
end

mutable struct ModelCommitment
  Vd :: Vector{Array{Float64}}
  Vk :: Vector{Array{Float64}}
  Vp :: Vector{Array{Float64}}

  cp :: Vector{Array{Float64}}
  ck :: Vector{Array{Float64}}
  xf′ :: Vector{Array{Float64}}
  ihkn :: Vector{Array{Int64}}
  tp :: Vector{Array{Float64}}

  np :: StylizedAltruism.NumPar
  mp :: StylizedAltruism.ModPar

  switch :: StylizedAltruism.switches

  function ModelCommitment(np::StylizedAltruism.NumPar,mp::StylizedAltruism.ModPar,switch::StylizedAltruism.switches)
    Vd = [[fill(-Inf64,np.nx,np.ny,np.nhi,np.nθ) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.ny,np.nh,np.no,np.ns,np.nθ)]]
    cp = [[fill(-Inf64,np.nx,np.ny,np.nhi,np.nθ) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.ny,np.nh,np.no,np.ns,np.nθ)]]
    ck = [[fill(-Inf64,np.nx,np.ny,np.nhi,np.nθ) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.ny,np.nh,np.no,np.ns,np.nθ)]]
    xd′ = [[fill(-Inf64,np.nx,np.ny,np.nhi,np.nθ) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.ny,np.nh,np.no,np.ns,np.nθ)]]
    ihkn = [[fill(-1,np.nx,np.ny,np.nhi,np.nθ) for ia = 1:np.na-1]; [fill(-1,np.nx,np.ny,np.nh,np.no,np.ns,np.nθ)]]
    tp = [[fill(-Inf64,np.nx,np.ny,np.nhi,np.nθ) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.ny,np.nh,np.no,np.ns,np.nθ)]]

    Vk =  [[fill(NaN64,np.nx,np.nx,np.ny,np.nhi,np.nθ) for ia = 1:np.na-1]; [fill(NaN64,np.nx,np.ny,np.nh,np.no,np.ns,np.nθ)]]
    Vp =  [[fill(NaN64,np.nx,np.nx,np.ny,np.nhi,np.nθ) for ia = 1:np.na-1]; [fill(NaN64,np.nx,np.nx,np.ny,np.nh,np.no,np.ns,np.nθ)]]

    new(Vd, Vk, Vp, cp, ck, xd′, ihkn, tp, np, mp, switch)
  end
end

"""
    function Vc_to_VkVp(Mc::ModelCommitment)

Converts the first-period commitment model solution to first-period first-stage
value functions for the kids and parents.
"""
function Vc_to_VkVp(Mc::ModelCommitment)
  np = Mc.np
  mp = Mc.mp

  @assert np.ny == 1 # Doesnt work with uncertainty
  @assert np.ns == 1 # Doesnt work with uncertainty
  @assert np.no == 1 # Doesnt work with owning
  ## Period 1
  for (ixp,xp) in enumerate(np.x_grd)
    for (ixk,xk) in enumerate(np.x_grd)
      for (iyk,yk) in enumerate(np.y_grd)
        for (ihki,hki) in enumerate(np.hi_grd)
          for (iθ,θ) in enumerate(np.θ_grd)
              # Find Policies
              ck1 = Mc.ck[1][ixk+ixp,iyk,ihki,iθ]
              cp1 = Mc.cp[1][ixk+ixp,iyk,ihki,iθ]
              ihkn = Mc.ihkn[1][ixk+ixp,iyk,ihki,iθ]
              hk1 = np.h_grd[ihkn]
              xf′ = Mc.xf′[1][ixk+ixp,iyk,ihki,iθ]
              ixfn = argmin(abs.(xf′ .- np.x_grd)) # Find next period starting position

              # Use policies to find first-period values
              valk = StylizedAltruism.utilk(ck1,hk1,mp.γ,mp.ξ)
              valp = StylizedAltruism.utilp(cp1,mp.γ)

              # Evaluate second period
              iokn = 1
              isn = 1
              ck2 = Mc.ck[2][ixfn,iyk,ihkn,iokn,isn,iθ]
              cp2 = Mc.cp[2][ixfn,iyk,ihkn,iokn,isn,iθ]
              hk2 = np.h_grd[Mc.ihkn[2][ixfn,iyk,ihkn,iokn,isn,iθ]]

              valk += mp.β*StylizedAltruism.utilk(ck2,hk2,mp.γ,mp.ξ)
              valp += mp.β*StylizedAltruism.utilp(cp2,mp.γ)
              valp += mp.η*valk

              Mc.Vk[1][ixk,ixp,iyk,ihki,iθ] = valk
              Mc.Vp[1][ixk,ixp,iyk,ihki,iθ] = valp

          end
        end
      end
    end
  end

  return Mc.Vk[1], Mc.Vp[1]
end

"""
    function Vk_to_Vkini(M::Model)

Converts the first-period second-stage value function of the kid to the
first-period first-stage value function (i.e. before the parent moves and
taking into accounts the parent's gift and savings decisions).
"""
function Vk_to_Vkini(M)
  np = M.np
  mp = M.mp

  Vk_itp = [LinearInterpolation((np.x_grd,np.x_grd),M.Vk[1][:,:,iyk,ihki],extrapolation_bc=Line()) for iyk = 1:np.ny, ihki = 1:np.nhi]
  hk_itp = [LinearInterpolation((np.x_grd,np.x_grd),np.h_grd[M.gke.h[1][:,:,iyk,ihki]],extrapolation_bc=Line()) for iyk = 1:np.ny, ihki = 1:np.nhi]
  ck_itp = [LinearInterpolation((np.x_grd,np.x_grd),M.gke.c[1][:,:,iyk,ihki],extrapolation_bc=Flat()) for iyk = 1:np.ny, ihki = 1:np.nhi]
  xk_itp = [LinearInterpolation((np.x_grd,np.x_grd),M.gke.x′[1][:,:,iyk,ihki],extrapolation_bc=Line()) for iyk = 1:np.ny, ihki = 1:np.nhi]
  Vk_ini = fill(NaN64,size(M.Vk[1]))
  hk_ini = fill(NaN64,size(M.gke.h[1]))
  ck_ini = fill(NaN64,size(M.gke.c[1]))
  xk_ini = fill(NaN64,size(M.gke.x′[1]))
  for (ixp,xp) in enumerate(np.x_grd)
    for (ixk,xk) in enumerate(np.x_grd)
      for (iyk,yk) in enumerate(np.y_grd)
        for (ihki,hki) in enumerate(np.hi_grd)
          tp = M.gp.t[1][ixp,ixk,iyk,ihki]
          xp′ = M.gp.x′[1][ixp,ixk,iyk,ihki]
          Vk_ini[ixp,ixk,iyk,ihki] = Vk_itp[iyk,ihki](xp′,xk+tp)
          ck_ini[ixp,ixk,iyk,ihki] = ck_itp[iyk,ihki](xp′,xk+tp)
          xk_ini[ixp,ixk,iyk,ihki] = xk_itp[iyk,ihki](xp′,xk+tp)
          hk_ini[ixp,ixk,iyk,ihki] = hk_itp[iyk,ihki](xp′,xk+tp)

        end
      end
    end
  end

  return Vk_ini, ck_ini, xk_ini, hk_ini
end

"""
    function lowest_pareto(Mcom,Mnocom)

finds the lowest pareto weights the kids and parent is willing to accept.
"""
function lowest_pareto(Mcom,Mnocom)
    np = Mnocom.np
    Vk_ini, _, _, _ = StylizedAltruism.Vk_to_Vkini(Mnocom)
    Vk,Vp = StylizedAltruism.Vc_to_VkVp(Mcom)

    ikpref = fill(-1,np.nx,np.nx,np.ny,np.nhi)
    ippref = fill(-1,np.nx,np.nx,np.ny,np.nhi)
    for ixk = 1:np.nx
        for ixp = 1:np.nx
            for iy = 1:np.ny,ihi=1:np.nhi
                ind = findfirst(Vk[ixk,ixp,iy,ihi,:] .>= Vk_ini[ixk,ixp,iy,ihi])
                if ind !=nothing
                  ikpref[ixk,ixp] = ind
                end

                ind = findlast(Vp[ixk,ixp,iy,ihi,:] .>= Mnocom.Vp[1][ixk,ixp,iy,ihi])
                if ind !=nothing
                 ippref[ixk,ixp] = ind
                end
            end
        end
    end

    return ikpref,ippref
end

"""
    function housing_commitment_kid(ikpreff,Mcom)

finds the housing choice for a kid using the lowest pareto-weight he is willing to accept
"""
function housing_commitment_kid(ikpref::Array{Int64},Mcom)
    np = Mcom.np
    hk = fill(NaN64,(np.nx,np.nx,np.ny,np.nhi))

    for (ixp,xp) in enumerate(np.x_grd)
      for (ixk,xk) in enumerate(np.x_grd)
        for (iyk,yk) in enumerate(np.y_grd)
          for (ihki,hki) in enumerate(np.hi_grd)
              iθ = ikpref[ixp,ixk,iyk,ihki]
              if iθ > 0
                  hk[ixp,ixk,iyk,ihki] = np.h_grd[Mcom.ihkn[1][ixk+ixp,iyk,ihki,iθ]]
              end
          end
        end
      end
    end

    return hk
end
