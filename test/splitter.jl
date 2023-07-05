tests = []
assert_tests = []

push!(assert_tests, (AssertionError, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A
    comps2 = @components B
    in = Stream(:in, fs, comps1, basis=FLOW)
    out = Stream(:out, fs, comps2)

    spl = Splitter(m, :spl, fs, in, [out])
    return nothing
end))

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A
    in = Stream(:in, fs, comps, basis=FLOW)
    out1 = Stream(:out1, fs, comps)
    out2 = Stream(:out2, fs, comps)

    spl = Splitter(m, :spl, fs, in, [out1, out2])

    set_silent(m)
    optimize!(m)
    return ( (termination_status(m) == LOCALLY_SOLVED) && 
                (value(m[:spl_out1_mass]) ≈ 0.5) &&
                (value(m[:spl_out2_mass]) ≈ 0.5) &&
                (value(m[:spl_out2_split_frac]) ≈ 0.5) &&
                (value(m[:spl_out1_A_massfrac]) ≈ 1.0) &&
                (value(m[:spl_out2_A_massfrac]) ≈ 1.0)
            )
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A
    in = Stream(:in, fs, comps, basis=FLOW)
    out1 = Stream(:out1, fs, comps)

    spl = Splitter(m, :spl, fs, in, [out1])
    set_silent(m)
    optimize!(m)
    
    return ( (termination_status(m) == LOCALLY_SOLVED) && 
                (value(m[:spl_out1_mass]) ≈ 1.0) &&
                (value(m[:spl_out1_A_massfrac]) ≈ 1.0) &&
                (value(m[:spl_out1_split_frac]) ≈ 1.0)
            )
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A B
    in = Stream(:in, fs, comps)
    out1 = Stream(:out1, fs, comps)
    out2 = Stream(:out2, fs, comps)

    spl = Splitter(m, :spl, fs, in, [out1, out2])
    fix(m[:spl_out1_split_frac], 0.0)
    
    set_silent(m)
    optimize!(m)
    return ( (termination_status(m) == LOCALLY_SOLVED) && 
                (value(m[:spl_out1_mass]) ≈ 0.0) &&
                (value(m[:spl_out1_A_massfrac]) ≈ 0.5) &&
                (value(m[:spl_out2_split_frac]) ≈ 1.0)
            )
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A B
    in = Stream(:in, fs, comps)
    out1 = Stream(:out1, fs, comps)
    out2 = Stream(:out2, fs, comps, basis=FLOW)

    spl = Splitter(m, :spl, fs, in, [out1, out2])
    fix(m[:spl_in_mass], 100)
    fix(m[:spl_in_A_massfrac], 0.6)
    fix(m[:spl_in_B_massfrac], 0.4)
    fix(m[:spl_out1_split_frac], 0.3)
    
    set_silent(m)
    optimize!(m)
    return ( (termination_status(m) == LOCALLY_SOLVED) && 
                (value(m[:spl_out1_mass]) ≈ 30.0) &&
                (value(m[:spl_out1_A_massfrac]) ≈ 0.6) &&
                (value(m[:spl_out1_B_massfrac]) ≈ 0.4) &&
                (value(m[:spl_out2_mass]) ≈ 70.0) &&
                (value(m[:spl_out2_A_mass]) ≈ 42.0) &&
                (value(m[:spl_out2_B_mass]) ≈ 28.0) &&
                (value(m[:spl_out2_split_frac]) ≈ 0.7)
            )
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A B
    in = Stream(:in, fs, comps, basis=FLOW)
    out1 = Stream(:out1, fs, comps, basis=FLOW)
    out2 = Stream(:out2, fs, comps, basis=FLOW)

    spl = Splitter(m, :spl, fs, in, [out1, out2])
    fix(m[:spl_in_mass], 100.0)
    fix(m[:spl_in_A_mass], 60.0)
    fix(m[:spl_in_B_mass], 40.0)
    fix(m[:spl_out1_split_frac], 0.3)
    
    set_silent(m)
    optimize!(m)
    return ( (termination_status(m) == LOCALLY_SOLVED) && 
                (value(m[:spl_out1_mass]) ≈ 30.0) &&
                (value(m[:spl_out1_A_mass]) ≈ 18.0) &&
                (value(m[:spl_out1_B_mass]) ≈ 12.0) &&
                (value(m[:spl_out2_mass]) ≈ 70.0) &&
                (value(m[:spl_out2_A_mass]) ≈ 42.0) &&
                (value(m[:spl_out2_B_mass]) ≈ 28.0) &&
                (value(m[:spl_out2_split_frac]) ≈ 0.7)
            )
end)

@testset "Splitter tests" begin
    for (e, t) in assert_tests
        @test_throws e t()
    end
    for t in tests
        @test t()
    end
end
nothing