function BenchmarkParameters(switch::StylizedAltruism.switches;
    η = 0.5, rf = 0.0, γ = 2.0, β = 0.9, ξ = 0.2,
    nx = 30, ntc = 60, ymin = 0.5, nxc = 80, ny = 3, xmax = 10, nh = 2)

  if switch.family == false
    η = 0.0
    ntc = 2
  end
  if switch.incomerisk == false
     ny = 1
   end
  if switch.incomes == false
    ny = 1
  end
  if switch.housing == false
    ξ = 0.0
    nh = 1
  else
  end

  mp = StylizedAltruism.ModPar(η=η,rf=rf,γ=γ,β=β,ξ=ξ)
  np = StylizedAltruism.NumPar(nx = nx, ntc = ntc, ymin=0.5,nxc=nxc,ny=ny,xmax=xmax,nh = nh, endowhouse = switch.endowhouse)
  if switch.incomes == false
    np.y_grd .= 0.0
    np.inc_grd .= 0.0
  end
  if switch.housing == false
    np.h_grd .= 0.0
  end


  return np,mp
end

function Solve(np,mp,switch)
  M = StylizedAltruism.Model(np,mp,switch)
  StylizedAltruism.SolveAHK!(M)
  return M
end
