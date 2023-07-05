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

    mix = Mixer(m, :mix, fs, [in], out)
    return nothing
end))

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A
    in = Stream(:in, fs, comps, basis=FLOW)
    out = Stream(:out, fs, comps)

    mix = Mixer(m, :mix, fs, [in], out)

    set_silent(m)
    optimize!(m)
    return ( (termination_status(m) == LOCALLY_SOLVED) && 
                (value(m[:mix_out_mass]) ≈ 1.0) &&
                (value(m[:mix_out_A_massfrac]) ≈ 1.0)
            )
end)

push!(tests,
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A
    in1 = Stream(:in1, fs, comps, basis=FLOW)
    in2 = Stream(:in2, fs, comps)
    out = Stream(:out, fs, comps)

    mix = Mixer(m, :mix, fs, [in1, in2], out)

    set_silent(m)
    optimize!(m)
    return ( (termination_status(m) == LOCALLY_SOLVED) && 
                (value(m[:mix_out_mass]) ≈ 2.0) &&
                (value(m[:mix_out_A_massfrac]) ≈ 1.0)
            )
end)

push!(tests,
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A B
    in1 = Stream(:in1, fs, comps, basis=FLOW)
    in2 = Stream(:in2, fs, comps)
    out = Stream(:out, fs, comps, basis=FLOW)

    mix = Mixer(m, :mix, fs, [in1, in2], out)

    fix(m[:mix_in1_mass],       300.0)
    fix(m[:mix_in1_A_mass],     100.0)
    fix(m[:mix_in1_B_mass],     200.0)
    fix(m[:mix_in2_mass],       200.0)
    fix(m[:mix_in2_A_massfrac], 0.3)
    fix(m[:mix_in2_B_massfrac], 0.7)
    
    set_silent(m)
    optimize!(m)
    # print_vars(m)
    return ( (termination_status(m) == LOCALLY_SOLVED) && 
                (value(m[:mix_out_mass]) ≈ 500.0) &&
                (value(m[:mix_out_A_mass]) ≈ 160.0) &&
                (value(m[:mix_out_B_mass]) ≈ 340.0)
             )
end)

push!(tests,
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A
    comps2 = @components B
    in1 = Stream(:in1, fs, comps1)
    in2 = Stream(:in2, fs, comps2)
    out = Stream(:out, fs, union(comps1, comps2))

    mix = Mixer(m, :mix, fs, [in1, in2], out)

    fix(m[:mix_in1_mass],       100.0)
    fix(m[:mix_in2_mass],       200.0)
    
    set_silent(m)
    optimize!(m)
 
    return ( (termination_status(m) == LOCALLY_SOLVED) && 
                (value(m[:mix_out_mass]) ≈ 300.0) &&
                (value(m[:mix_out_A_massfrac]) ≈ 0.333333333) &&
                (value(m[:mix_out_B_massfrac]) ≈ 0.666666667)
             )
end)

push!(tests,
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    in1 = Stream(:in1, fs, @components A)
    in2 = Stream(:in2, fs, @components B)
    in3 = Stream(:in3, fs, @components C)
    out = Stream(:out, fs, @components A B C)

    mix = Mixer(m, :mix, fs, [in1, in2, in3], out)

    fix(m[:mix_in1_mass],       100.0)
    fix(m[:mix_in2_mass],       100.0)
    fix(m[:mix_in3_mass],       100.0)
    
    set_silent(m)
    optimize!(m)
 
    return ( (termination_status(m) == LOCALLY_SOLVED) && 
                (value(m[:mix_out_mass]) ≈ 300.0) &&
                (value(m[:mix_out_A_massfrac]) ≈ 0.333333333) &&
                (value(m[:mix_out_B_massfrac]) ≈ 0.333333333) &&
                (value(m[:mix_out_C_massfrac]) ≈ 0.333333333)
             )
end)

push!(tests,
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A B C
    in1 = Stream(:in1, fs, comps)
    in2 = Stream(:in2, fs, comps)
    in3 = Stream(:in3, fs, comps, basis=FLOW)
    out = Stream(:out, fs, comps)

    mix = Mixer(m, :mix, fs, [in1, in2, in3], out)

    fix(m[:mix_in1_mass],       100.0)
    fix(m[:mix_in1_A_massfrac], 0.1)
    fix(m[:mix_in1_B_massfrac], 0.4)
    fix(m[:mix_in1_C_massfrac], 0.5)
    fix(m[:mix_in2_mass],       100.0)
    fix(m[:mix_in2_A_massfrac], 0.3)
    fix(m[:mix_in2_B_massfrac], 0.5)
    fix(m[:mix_in2_C_massfrac], 0.2)
    fix(m[:mix_in3_mass],       100.0)
    fix(m[:mix_in3_A_mass],     25.0)
    fix(m[:mix_in3_B_mass],     0.0)
    fix(m[:mix_in3_C_mass],     75.0)
    
    set_silent(m)
    optimize!(m)
 
    return ( (termination_status(m) == LOCALLY_SOLVED) && 
                (value(m[:mix_out_mass]) ≈ 300.0) &&
                (value(m[:mix_out_A_massfrac]) ≈ 0.216666667) &&
                (value(m[:mix_out_B_massfrac]) ≈ 0.3) &&
                (value(m[:mix_out_C_massfrac]) ≈ 0.483333333)
             )
end)

push!(tests,
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    in1 = Stream(:in1, fs, @components A)
    in2 = Stream(:in2, fs, @components B)
    out = Stream(:out, fs, @components A B)

    mix = Mixer(m, :mix, fs, [in1, in2], out)

    fix(m[:mix_in1_mass],       100.0)
    fix(m[:mix_in2_mass],       0.0)
    
    set_silent(m)
    optimize!(m)
 
    return ( (termination_status(m) == LOCALLY_SOLVED) && 
                (value(m[:mix_out_mass]) ≈ 100.0) &&
                (value(m[:mix_out_A_massfrac]) ≈ 1.0) &&
                (value(m[:mix_out_B_massfrac]) ≈ 0.0)
             )
end)

@testset "Mixer tests" begin
    for (e, t) in assert_tests
        @test_throws e t()
    end
    for t in tests
        @test t()
    end
end
nothing