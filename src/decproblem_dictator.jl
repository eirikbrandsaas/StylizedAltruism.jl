function SolveVd2(np,mp)
  cp_pol = fill(-1.0,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)
  ck_pol = fill(-1.0,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)
  ihkn_pol = fill(-1,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)
  val = fill(-1.0,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)
  vtmp = fill(-Inf64,np.nh,np.ntc)
  ia = np.na
  for (ixp,xp) in enumerate(np.x_grd)
    for (ixk,xk) in enumerate(np.x_grd)
      for (iyk,yk) in enumerate(np.y_grd)
        for (ihk,hk) in enumerate(np.h_grd)
          for (io,ok) in enumerate(np.o_grd)
            for (is,sk) in enumerate(np.s_grd)
              ## Loop over possible choices
              vtmp .= -1e64
              xpn = 0.0
              wk = StylizedAltruism.income(yk,ia,np)
              for (ihkn,hkn) in enumerate(np.h_grd)
                for (ick,ck) in enumerate(np.tc_grd)
                  cp = StylizedAltruism.ck_bc2(xk+xp-ck,wk,xpn,hk,hkn,mp.rf,mp.κ,mp.os,mp.rs,sk,ok)
                  if ck > 0.0 && cp >0.0
                    vtmp[ihkn,ick] =  StylizedAltruism.utilp(cp,mp.γ) + mp.η*StylizedAltruism.utilk(ck,hkn,mp.γ,mp.ξ)
                  end
                end
              end
              imax = argmax(vtmp)
              ihkn = imax[1]
              ick = imax[2]
              hkn = np.h_grd[ihkn]
              ck = np.tc_grd[ick]
              cp_pol[ixp,ixk,iyk,ihk,io,is] = StylizedAltruism.ck_bc2(xk+xp-ck,wk,xpn,hk,hkn,mp.rf,mp.κ,mp.os,mp.rs,sk,ok)
              ck_pol[ixp,ixk,iyk,ihk,io,is] = ck
              ihkn_pol[ixp,ixk,iyk,ihk,io,is] = ihkn
              val[ixp,ixk,iyk,ihk,io,is] = vtmp[imax]
            end
          end
        end
      end
    end
  end

  return cp_pol, ck_pol,ihkn_pol,val
end

function SolveVd1(np,mp,valnxt)
  cp_pol = fill(-2.0,np.nx,np.nx,np.ny,np.nhi)
  ck_pol = fill(-1.0,np.nx,np.nx,np.ny,np.nhi)
  xp′_pol = fill(-1.0,np.nx,np.nx,np.ny,np.nhi)
  ihkn_pol = fill(-1,np.nx,np.nx,np.ny,np.nhi)
  val = fill(-1.0,np.nx,np.nx,np.ny,np.nhi)
  vtmp = fill(-Inf64,np.nh,np.ntc,np.nxc)
  ia = np.na

  @assert np.no == 1 # Only works with one ownership status

  vnxt = [LinearInterpolation((np.x_grd,np.x_grd),valnxt[:,:,iyk,ihkn,io,is],extrapolation_bc=Line()) for iyk = 1:np.ny, ihkn = 1:np.nh, io = 1:np.no, is=1:np.ns]

  for (ixp,xp) in enumerate(np.x_grd)
    for (ixk,xk) in enumerate(np.x_grd)
      for (iyk,yk) in enumerate(np.y_grd)
        for (ihk,hk) in enumerate(np.hi_grd)
          ## Loop over possible choices
          vtmp .= -1e64
          wk = StylizedAltruism.income(yk,ia,np)
          Threads.@threads for (ihkn,hkn) in collect(enumerate(np.h_grd))
            for (ick,ck) in enumerate(np.tc_grd)
              for (ixpn,xpn) in enumerate(np.xc_grd)
                 vtmp[ihkn,ick,ixpn] = Vd1_eval(xp,xk,ck,hkn,xpn,iyk,ihkn,mp,np,vnxt)
              end
            end
          end
          imax = argmax(vtmp)
          ihkn = imax[1]
          ick = imax[2]

          ixpn = imax[3]
          hkn = np.h_grd[ihkn]
          ck = np.tc_grd[ick]
          xpn = np.xc_grd[ixpn]
          cp_pol[ixp,ixk,iyk,ihk] = StylizedAltruism.cp_bc(xp + xk,xpn,0.0,mp.rf)
          ck_pol[ixp,ixk,iyk,ihk] = ck
          ihkn_pol[ixp,ixk,iyk,ihk] = ihkn
          xp′_pol[ixp,ixk,iyk,ihk] = xpn
          val[ixp,ixk,iyk,ihk] = vtmp[imax]
        end
      end
    end
  end


  return cp_pol, ck_pol,ihkn_pol, xp′_pol, val
end

function Vd1_eval(xp,xk,ck,hkn,xpn,iyk,ihkn,mp,np,vnxt)
  cp = StylizedAltruism.cp_bc(xp + xk-ck-hkn,xpn,0.0,mp.rf)
  iokn = 1
  vtmp = 0.0
  if ck > 0.0 && cp >0.0 && xpn > 0.0
    vtmp =  StylizedAltruism.utilp(cp,mp.γ) + mp.η*StylizedAltruism.utilk(ck,hkn,mp.γ,mp.ξ)
    for iykn = 1:np.ny
      for is = 1:np.ns
          vtmp += mp.β*np.Πy[iyk,iykn]*np.Πs[is]*vnxt[iykn,ihkn,iokn,is](xpn,0.0)
      end
    end
  else
    vtmp = -1e60
  end

  return   vtmp
end


mutable struct ModelDictator
  Vd :: Vector{Array{Float64}}

  cp :: Vector{Array{Float64}}
  ck :: Vector{Array{Float64}}
  xd′ :: Vector{Array{Float64}}
  ihkn :: Vector{Array{Int64}}
  tp :: Vector{Array{Float64}}

  np :: StylizedAltruism.NumPar
  mp :: StylizedAltruism.ModPar

  switch :: StylizedAltruism.switches

  function ModelDictator(np::StylizedAltruism.NumPar,mp::StylizedAltruism.ModPar,switch::StylizedAltruism.switches)
    Vd = [[fill(-Inf64,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)]]
    cp = [[fill(-Inf64,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)]]
    ck = [[fill(-Inf64,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)]]
    xd′ = [[fill(-Inf64,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)]]
    ihkn = [[fill(-1,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(-1,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)]]
    tp = [[fill(-Inf64,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)]]

    new(Vd, cp, ck, xd′, ihkn, tp, np, mp, switch)
  end
end

function SolveDictator(np,mp,switch)
  Md = ModelDictator(np,mp,switch)

  Md.cp[2], Md.ck[2], Md.ihkn[2], Md.Vd[2] = SolveVd2(np,mp)
  Md.cp[1], Md.ck[1], Md.ihkn[1], Md.xd′[1], Md.Vd[1] = SolveVd1(np,mp,Md.Vd[2])

  return Md
end
