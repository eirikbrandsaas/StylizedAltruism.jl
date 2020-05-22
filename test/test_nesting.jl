function test_altr_vs_eta()
    nx = 25
    ntc = 31
    nxc = 22

    # Model turning of transfers/family stuff explicitly
    switch = switches(family=false,incomes=true,incomerisk=true,housing=true)
    np, mp = BenchmarkParameters(switch, nx = nx, ntc = ntc, nxc = nxc)
    M1 = Solve(np,mp,switch)

    # Setting the altruism parameter = 0
    switch = switches(family=true,incomes=true,incomerisk=true,housing=true)
    np, mp = BenchmarkParameters(switch, nx = nx, ntc = ntc, nxc = nxc)
    mp.η= 0.0
    M2 = Solve(np,mp,switch)

    sum_identical = sum([1.0 .- sum(M1.gp.t[ia] .== M2.gp.t[ia])/length(M2.gp.t[ia]) for ia = 1:np.na])/np.na

    return sum_identical
end
