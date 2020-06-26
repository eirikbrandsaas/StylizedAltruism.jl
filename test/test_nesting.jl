function test_altr_vs_eta()
    nx = 25
    ntc = 31
    nxc = 22

    # Model turning of transfers/family stuff explicitly
    switch = StylizedAltruism.switches(family=false,incomes=true,incomerisk=true,housing=true)
    np, mp = BenchmarkParameters(switch, nx = nx, ntc = ntc, nxc = nxc)
    M1 = Solve(np,mp,switch)

    # Setting the altruism parameter = 0
    switch = StylizedAltruism.switches(family=true,incomes=true,incomerisk=true,housing=true)
    np, mp = BenchmarkParameters(switch, nx = nx, ntc = ntc, nxc = nxc)
    mp.η= 0.0
    M2 = Solve(np,mp,switch)

    sum_identical = sum([1.0 .- sum(M1.gp.t[ia] .== M2.gp.t[ia])/length(M2.gp.t[ia]) for ia = 1:np.na])/np.na

    return sum_identical
end

function test_houserisk()
    nx = 25
    nh = 1
    ntc = 50
    nxc = 50
    ny = 1

    save = false # Switch to save plots

    switch = StylizedAltruism.switches(family = true,incomes = true, incomerisk = true, housing = true,endowhouse=false,pricerisk=true)
    np,mp = BenchmarkParameters(switch, nx = nx, η = 0.5,nh = nh, nxc = nxc, ntc = ntc,ξ = 0.4, xmax=13, ny=ny,ns=5,no=2)
    np.s_grd .= 0
    M1 = Solve(np,mp,switch)

    switch = StylizedAltruism.switches(family = true,incomes = true, incomerisk = true, housing = true,endowhouse=false,pricerisk=true)
    np,mp = BenchmarkParameters(switch, nx = nx, η = 0.5,nh = nh, nxc = nxc, ntc = ntc,ξ = 0.4, xmax=13, ny=ny,ns=1,no=2)
    M2 = Solve(np,mp,switch)

    M1.Vk[1] ≈  M2.Vk[1]

end

test_houserisk()
