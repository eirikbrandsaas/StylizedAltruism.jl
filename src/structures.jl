## File to hold all structures

mutable struct ModPar
  β :: Float64
  γ :: Float64
  η :: Float64
  ξ :: Float64
  χ :: Float64
  rf :: Float64
  κ :: Float64
  p :: Float64
  os :: Float64
  rs :: Float64

  function ModPar(;β=0.92, γ=1.5, η=0.3, ξ = 0.10, χ = .0, rf = 0.0, κ = 0.0,p=2.0,os = 1.0, rs =1.0)
    new(β, γ, η, ξ, χ, rf, κ, p, os, rs)
  end

end

mutable struct NumPar
  na :: Int64 # Ages
  nx :: Int64
  ny :: Int64
  nxc :: Int64
  ntc :: Int64
  nh :: Int64
  nhi :: Int64 # Size of housing endowment
  ns :: Int64
  no :: Int64
  x_grd :: Vector{Float64} # Stategrid (common for kids and parents)
  y_grd :: Vector{Float64} # Stategrid (only for kids)
  xc_grd :: Vector{Float64} # Choicegrid (common for kids and parents (for simplicity))
  tc_grd :: Vector{Float64} # Choicegrid (only for parents)
  inc_grd :: Vector{Float64} # Income grid (kids)
  h_grd :: Vector{Float64} # Housing grids
  hi_grd :: Vector{Float64} # Housing Endowment grid
  s_grd :: Vector{Float64} # House price uncertainty
  o_grd :: Vector{Bool}

  Πy :: Array{Float64,2} # Probability for y' given y
  Πs :: Array{Float64,1}

  function NumPar(;nh=1,na=2, nx=11, ny=5 ,nxc=20, ntc=20, ns = 1, xmax = 10.0, ymin = 0.5, ymax = 2.5,hmin=0.2,hmax = 1.0, endowhouse = false,no=1)
    x_grd = range(1e-5,stop=xmax,length=nx)
    if ny > 1
      y_grd = range(ymin,stop=ymax,length=ny)
    else
      y_grd = range(1.0,stop=1.0,length=ny)
    end
    xc_grd = range(minimum(x_grd),stop=maximum(x_grd),length=nxc)
    tc_grd = range(0.0,stop=maximum(x_grd)/2,length=ntc)

    if ns == 1
      s_grd = range(0.0,stop=0.0,length=ns)
    else
      s_grd = range(-0.3,stop=0.3,length=ns)
    end

    if no == 1
      o_grd = [false]
    elseif no ==2
      o_grd = [false,true]
    else
      throw("no must be 1 or 2")
    end
    Πs = fill(1.0/ns,(ns))

    Πy = fill(1.0/ny,(ny,ny))
    if ny == 3
      Πy = [2/3 1/3 0; 1/3 1/3 1/3; 0 1/3 2/3]
    end

    inc_grd = range(1.,stop = 1.0,length=na)
    if nh > 1
      h_grd = range(hmin,stop=hmax,length=nh)
    else
      h_grd = range(minimum(x_grd),stop=minimum(x_grd),length=nh)
    end

    if endowhouse == false
      nhi = 1
      hi_grd = range(0.0,stop=0.0,length=nhi)
    else
      nhi = nh
      hi_grd = h_grd
    end

    new(na,nx, ny, nxc, ntc, nh, nhi, ns, no,x_grd, y_grd, xc_grd, tc_grd, inc_grd, h_grd, hi_grd, s_grd, o_grd, Πy, Πs)
  end
end

mutable struct Polp
  c :: Vector{Array{Float64}}
  x′ :: Vector{Array{Float64}}
  t :: Vector{Array{Float64}}

  function Polp(np::NumPar)
    c = [[fill(0.0,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(0.0,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)]]
    x′ = [[fill(0.0,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(0.0,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)]]
    t = [[fill(0.0,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(0.0,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)]]

    new(c,x′,t)
  end
end

mutable struct Polk # Policy functions for the kids, for each housing choice!
  c :: Vector{Array{Float64}}
  x′ :: Vector{Array{Float64}}
  h :: Vector{Array{Int64}}
  o :: Vector{Array{Int64}}
  disc :: Vector{Array{Int64}}

  function Polk(np::NumPar)
    c = [[fill(-Inf64,np.nx,np.nx,np.ny,np.nhi,np.nh,np.no) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.ny,np.nh,np.no,np.ns,np.nh)]]
    x′ = [[fill(0.0,np.nx,np.nx,np.ny,np.nhi,np.nh,np.no) for ia = 1:np.na-1]; [fill(0.0,np.nx,np.ny,np.nh,np.no,np.ns,np.nh)]]
    h = [[fill(0,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(0,np.nx,np.ny,np.nh,np.no,np.ns)]]
    o = [[fill(0,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(0,np.nx,np.ny,np.nh,np.no,np.ns)]]
    disc = [[fill(0,np.nx,np.nx,np.ny,np.nhi,np.nh,np.no) for ia = 1:np.na-1]; [fill(0,np.nx,np.ny,np.nh,np.no,np.ns,np.nh)]]
    new(c,x′,h,o,disc)
  end
end


mutable struct Polk_eq # Policy functinos for the kids
  c :: Vector{Array{Float64}}
  x′ :: Vector{Array{Float64}}
  h :: Vector{Array{Int64}}

  function Polk_eq(np::NumPar)
    c = [[fill(-Inf64,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.ny,np.nh,np.no,np.ns)]]
    x′ = [[fill(0.0,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(0.0,np.nx,np.ny,np.nh,np.no,np.ns)]]
    h = [[fill(0.0,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(0.0,np.nx,np.ny,np.nh,np.no,np.ns)]]
    new(c,x′,h)
  end
end

mutable struct switches

  family :: Bool
  housing :: Bool
  endowhouse :: Bool
  incomerisk :: Bool
  incomes :: Bool
  pricerisk :: Bool

  function switches(;family=false, endowhouse=false, housing=false, incomerisk=false,incomes=false,pricerisk=false)

  new(family, housing, endowhouse,incomerisk,incomes,pricerisk)
  end
end


mutable struct Model
  Vk :: Vector{Array{Float64}}
  Vp :: Vector{Array{Float64}}

  gk :: Polk
  gp :: Polp
  gke :: Polk_eq

  np :: NumPar
  mp :: ModPar

  switch ::  switches

  function Model(np::NumPar,mp::ModPar,switch::switches)

    Vk =  [[fill(-Inf64,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.ny,np.nh,np.no,np.ns)]]
    Vp =  [[fill(-Inf64,np.nx,np.nx,np.ny,np.nhi) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.nx,np.ny,np.nh,np.no,np.ns)]]

    gk = Polk(np)
    gp = Polp(np)
    gke = Polk_eq(np)

    switch = switch

    new(Vk, Vp, gk ,gp, gke, np, mp, switch)
  end

end
