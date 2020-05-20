function TestNumericalSolution(M)
    # Tests if the maximum transfer choice is too low
    @assert maximum(maximum.(M.gp.t)) <= maximum(M.np.tc_grd)

    # Tests that household savings is never bound by the choice grid
    @assert maximum(maximum.(M.gp.x′)) <= maximum(M.np.xc_grd)
    @assert maximum(maximum.(M.gk.x′)) <= maximum(M.np.xc_grd)
end
