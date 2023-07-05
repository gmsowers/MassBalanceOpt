function test1()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A
    comps2 = @components B
    in1 = Stream(:in1, fs, comps1)
    out1 = Stream(:out1, fs, comps2)

    stoic_coef = @stoic(A => B)
    mw = Dict(:A => 1.0, :B => 1.0)
    arx = StoicReactor(m, :arx, fs, in1, out1, stoic_coef, mw)
    set_start_value(m[:arx_extent_rx_1], 0.5)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:arx_out1_B_massfrac]) ≈ 1.0)
                ≈(value(m[:arx_out1_mass]), 1.0, rtol=1.0e-8)
                ≈(value(m[:arx_extent_rx_1]), 1.0, rtol=1.0e-8)
            ])
end

function test2()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A
    comps2 = @components A B
    in1 = Stream(:in1, fs, comps1)
    out1 = Stream(:out1, fs, comps2)

    stoic_coef = @stoic A => B
    mw = Dict(:A => 1.0, :B => 1.0)
    conv = OrderedDict(1 => (c=:A, X=0.0))
    arx = StoicReactor(m, :arx, fs, in1, out1, stoic_coef, mw, conv)
    set_start_value(m[:arx_out1_A_mass], 0.5)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:arx_out1_A_massfrac]) ≈ 1.0)
                ≈(value(m[:arx_out1_B_massfrac]), 0.0, rtol=1.0e-8)
                ≈(value(m[:arx_out1_mass]), 1.0, rtol=1.0e-8)
                ≈(value(m[:arx_extent_rx_1]), 0.0, rtol=1.0e-8)
            ])
end

function test3()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A
    comps2 = @components A B
    in1 = Stream(:in1, fs, comps1)
    out1 = Stream(:out1, fs, comps2)

    stoic_coef = @stoic A => B
    mw = Dict(:A => 1.0, :B => 1.0)
    conv = OrderedDict(1 => (c=:A, X=0.5))
    @block(arx, StoicReactor, in1, out1, stoic_coef, mw, conv)
    set_start_value(m[:arx_out1_A_mass], 0.25)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:arx_out1_A_massfrac]) ≈ 0.5)
                ≈(value(m[:arx_out1_B_massfrac]), 0.5, rtol=1.0e-8)
                ≈(value(m[:arx_out1_mass]), 1.0, rtol=1.0e-8)
                ≈(value(m[:arx_extent_rx_1]), 0.5, rtol=1.0e-8)
            ])
end

function test4()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A B
    comps2 = @components B C
    in1 = Stream(:in1, fs, comps1, basis=FLOW)
    out1 = Stream(:out1, fs, comps2, basis=FLOW)

    stoic_coef = @stoic A + B => C
    mw = Dict(:A => 1.0, :B => 1.0, :C => 2.0)
    arx = StoicReactor(m, :arx, fs, in1, out1, stoic_coef, mw)
    set_start_value(m[:arx_out1_mass], 0.5)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:arx_out1_C_mass]) ≈ 1.0)
                ≈(value(m[:arx_out1_B_mass]), 0.0, atol=1.0e-8)
                ≈(value(m[:arx_out1_mass]), 1.0, rtol=1.0e-8)
                ≈(value(m[:arx_extent_rx_1]), 0.5, rtol=1.0e-8)
            ])
end

function test5()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components C2H2 H2
    comps2 = @components C2H2 H2 C2H4
    in1 = Stream(:in1, fs, comps1)
    out1 = Stream(:out1, fs, comps2)

    stoic_coef = @stoic C2H2 + H2 => C2H4
    mw = Dict(:C2H2 => 26.03728, :H2 => 2.01588, :C2H4 => 28.05316)
    conv = OrderedDict(1 => (c=:C2H2, X=0.0))
    arx = StoicReactor(m, :arx, fs, in1, out1, stoic_coef, mw, conv)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:arx_out1_C2H2_massfrac]) ≈ 0.5)
                (value(m[:arx_out1_H2_massfrac]) ≈ 0.5)
            ])
end

function test6()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components C2H2 H2 N2
    comps2 = @components C2H2 H2 C2H4 N2
    in1 = Stream(:in1, fs, comps1)
    out1 = Stream(:out1, fs, comps2)

    stoic_coef = @stoic C2H2 + H2 => C2H4
    mw = Dict(:C2H2 => 26.03728, :H2 => 2.01588, :C2H4 => 28.05316)
    conv = OrderedDict(1 => (c=:C2H2, X=1.0))
    arx = StoicReactor(m, :arx, fs, in1, out1, stoic_coef, mw, conv)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                ≈(value(m[:arx_out1_C2H2_massfrac]), 0.0, atol=1.0e-8)
                ≈(value(m[:arx_out1_H2_moles]), 0.1525516013, rtol=1.0e-8)
                (value(m[:arx_out1_C2H4_moles]) ≈ 0.012802156497)
                ≈(value(m[:arx_out1_N2_mass]), 0.33333333333, rtol=1.0e-8)
            ])
end

function test7()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components C2H2 H2
    comps2 = @components C2H2 H2 C2H4 C2H6
    in1 = Stream(:in1, fs, comps1, basis=FLOW)
    out1 = Stream(:out1, fs, comps2, basis=FLOW)

    stoic_coef = @stoic begin
                            C2H2 + H2  => C2H4
                            C2H2 + 2H2 => C2H6
                         end
    mw = Dict(:C2H2 => 26.03728, :H2 => 2.01588, :C2H4 => 28.05316, :C2H6 => 30.06904)
    conv = OrderedDict(1 => (c=:C2H2, X=0.7),
                       2 => (c=:C2H2, X=0.3))
    arx = StoicReactor(m, :arx, fs, in1, out1, stoic_coef, mw, conv)

    @values begin
        arx_in1_C2H2_mass = 100.0
        arx_in1_H2_mass   = 100.0
        arx_in1_mass      = 200.0
    end
    
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                ≈(value(m[:arx_out1_C2H2_mass]), 0.0, atol=1.0e-8)
                ≈(value(m[:arx_out1_H2_moles]), 44.613, rtol=1.0e-4)
                ≈(value(m[:arx_out1_C2H4_moles]), 2.6885, rtol=1.0e-4)
            ])
end

@testset "StoicReactor tests" begin
    @test test1()
    @test test2()
    @test test3()
    @test test4()
    @test test5()
    @test test6()
    @test test7()
end
nothing