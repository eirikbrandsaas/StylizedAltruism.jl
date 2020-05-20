## File to hold all structures

mutable struct ModPar
  β :: Float64
  γ :: Float64
  η :: Float64
  ξ :: Float64
  χ :: Float64
  rf :: Float64


  function ModPar(;β=0.92, γ=1.5, η=0.3, ξ = 0.10, χ = .0, rf = 0.0)
    new(β, γ, η, ξ, χ, rf)
  end

end

mutable struct NumPar
  na :: Int64 # Ages
  nx :: Int64
  ny :: Int64
  nxc :: Int64
  ntc :: Int64
  x_grd :: Vector{Float64} # Stategrid (common for kids and parents)
  y_grd :: Vector{Float64} # Stategrid (only for kids)
  xc_grd :: Vector{Float64} # Choicegrid (common for kids and parents (for simplicity))
  tc_grd :: Vector{Float64} # Choicegrid (only for parents)
  inc_grd :: Vector{Float64}

  Πy :: Array{Float64,2} # Probability for y' given y

  function NumPar(;na=2, nx=11, ny=5 ,nxc=20, ntc=20, xmax = 10.0, ymin = 0.5, ymax = 2.5)
    x_grd = range(1e-5,stop=xmax,length=nx)
    if ny > 1
      y_grd = range(ymin,stop=ymax,length=ny)
    else
      y_grd = range(1.0,stop=1.0,length=ny)
    end
    xc_grd = range(minimum(x_grd),stop=maximum(x_grd),length=nxc)
    tc_grd = range(0.0,stop=maximum(x_grd),length=ntc)

    Πy = fill(1.0/ny,(ny,ny))
    if ny == 3
      Πy = [2/3 1/3 0; 1/3 1/3 1/3; 0 1/3 2/3]
    end

    inc_grd = range(1.,stop = 1.0,length=na)

    new(na,nx, ny, nxc, ntc, x_grd, y_grd, xc_grd, tc_grd, inc_grd, Πy)
  end
end

mutable struct Polp
  c :: Vector{Array{Float64}}
  x′ :: Vector{Array{Float64}}
  t :: Vector{Array{Float64}}

  function Polp(np::NumPar)
    c = [fill(0.0,np.nx,np.nx,np.ny) for ia = 1:np.na]
    x′ = [fill(0.0,np.nx,np.nx,np.ny) for ia = 1:np.na]
    t = [fill(0.0,np.nx,np.nx,np.ny) for ia = 1:np.na]

    new(c,x′,t)
  end
end


mutable struct Polk
  c :: Vector{Array{Float64}}
  x′ :: Vector{Array{Float64}}

  function Polk(np::NumPar)
    c = [[fill(0.0,np.nx,np.nx,np.ny) for ia = 1:np.na-1]; [fill(0.0,np.nx,np.ny)]]
    x′ = [[fill(0.0,np.nx,np.nx,np.ny) for ia = 1:np.na-1]; [fill(0.0,np.nx,np.ny)]]

    new(c,x′)
  end
end

mutable struct Model
  Vk :: Vector{Array{Float64}}
  Vp :: Vector{Array{Float64}}

  gk :: Polk
  gp :: Polp

  np :: NumPar
  mp :: ModPar

  function Model(np::NumPar,mp::ModPar)

    Vk =  [[fill(-Inf64,np.nx,np.nx,np.ny) for ia = 1:np.na-1]; [fill(-Inf64,np.nx,np.ny)]]
    Vp =  [fill(-Inf64,np.nx,np.nx,np.ny) for ia = 1:np.na]

    gk = Polk(np)
    gp = Polp(np)

    new(Vk, Vp, gk ,gp, np, mp)
  end

end
