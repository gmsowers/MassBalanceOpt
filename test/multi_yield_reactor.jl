tests = []
assert_tests = []

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A
    comps2 = @components B C
    in = Stream(:in, fs, comps1)
    out = Stream(:out, fs, comps2)

    r1 = MultiYieldReactor(m, :r1, fs, [in], [out], [:pureA], :reactor)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:r1_out_mass]) ≈ 1.0)
                (value(m[:r1_out_C_massfrac]) ≈ 1.0)
            ])
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A
    comps2 = @components B
    comps3 = @components C D
    eth = Stream(:eth, fs, comps1)
    pro = Stream(:pro, fs, comps2)
    cgeth = Stream(:cgeth, fs, comps3)
    cgpro = Stream(:cgpro, fs, comps3)

    r1 = MultiYieldReactor(m, :r1, fs, [eth, pro], [cgeth, cgpro], [:eth, :pro], :furn)
    fix(m[:r1_eth_mass], 100.0)
    fix(m[:r1_pro_mass], 50.0)
    fix(m[:r1_eth_rate], 25.0)
    fix(m[:r1_pro_rate], 25.0)
    fix(m[:r1_eth_y_C_from_A], 0.55)
    fix(m[:r1_pro_y_C_from_B], 0.4)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:r1_total_feed_mass]) ≈ 150.0)
                (value(m[:r1_cgeth_C_massfrac]) ≈ 0.55)
                (value(m[:r1_cgpro_C_massfrac]) ≈ 0.4)
                (value(m[:r1_eth_n_furn]) ≈ 4.0)
                (value(m[:r1_pro_n_furn]) ≈ 2.0)
            ])
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A
    comps2 = @components B
    comps3 = @components C D
    eth = Stream(:eth, fs, comps1)
    pro = Stream(:pro, fs, comps2)
    cgeth = Stream(:cgeth, fs, comps3)
    cgpro = Stream(:cgpro, fs, comps3, basis=FLOW)

    r1 = MultiYieldReactor(m, :r1, fs, [eth, pro], [cgeth, cgpro], [:eth, :pro], :furn)
    fix(m[:r1_eth_mass], 100.0)
    fix(m[:r1_pro_mass], 0.0)
    fix(m[:r1_eth_rate], 25.0)
    fix(m[:r1_pro_rate], 25.0)
    fix(m[:r1_eth_y_C_from_A], 0.55)
    fix(m[:r1_pro_y_C_from_B], 0.4)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:r1_total_feed_mass]) ≈ 100.0)
                (value(m[:r1_cgeth_C_massfrac]) ≈ 0.55)
                (value(m[:r1_eth_n_furn]) ≈ 4.0)
                (isapprox(value(m[:r1_pro_n_furn]), 0.0, atol=1.0e-8))
            ])
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components C2H6
    comps2 = @components C3H8
    comps3 = @components OFFGAS C2H4 C2H6 C3H6 C3H8 C4PLUS
    eth = Stream(:eth, fs, comps1)
    pro = Stream(:pro, fs, comps2)
    cgeth = Stream(:cgeth, fs, comps3, basis=FLOW)
    cgpro = Stream(:cgpro, fs, comps3, basis=FLOW)

    r1 = MultiYieldReactor(m, :r1, fs, [eth, pro], [cgeth, cgpro], [:eth, :pro], :furn)
    fix(m[:r1_eth_mass], 100.0)
    fix(m[:r1_pro_mass], 25.0)
    fix(m[:r1_eth_rate], 25.0)
    fix(m[:r1_pro_rate], 25.0)
    fix(m[:r1_eth_y_OFFGAS_from_C2H6], 0.10)
    fix(m[:r1_eth_y_C2H4_from_C2H6], 0.55)
    fix(m[:r1_eth_y_C3H6_from_C2H6], 0.01)
    fix(m[:r1_eth_y_C3H8_from_C2H6], 0.004)
    fix(m[:r1_eth_y_C4PLUS_from_C2H6], 0.01)
    fix(m[:r1_pro_y_OFFGAS_from_C3H8], 0.17)
    fix(m[:r1_pro_y_C2H4_from_C3H8], 0.38)
    fix(m[:r1_pro_y_C2H6_from_C3H8], 0.05)
    fix(m[:r1_pro_y_C3H6_from_C3H8], 0.19)
    fix(m[:r1_pro_y_C4PLUS_from_C3H8], 0.1)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:r1_total_feed_mass]) ≈ 125.0)
                (value(m[:r1_eth_y_C2H6_from_C2H6]) ≈ 0.326)
                (value(m[:r1_pro_y_C3H8_from_C3H8]) ≈ 0.11)
                (value(m[:r1_cgeth_C2H4_mass]) ≈ 55.0)
                (value(m[:r1_cgpro_C2H4_mass]) ≈ 9.5)
            ])
end)

@testset "MultiYieldReactor tests" begin
    for (e, t) in assert_tests
        @test_throws e t()
    end
    for t in tests
        @test t()
    end
end
nothing