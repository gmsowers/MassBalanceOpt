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

    r1 = YieldReactor(m, :r1, fs, in, out)
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
    comps2 = @components A B C
    in = Stream(:in, fs, comps1)
    out = Stream(:out, fs, comps2)

    r1 = YieldReactor(m, :r1, fs, in, out)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:r1_out_mass]) ≈ 1.0)
                (value(m[:r1_out_A_massfrac]) ≈ 1.0)
            ])
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A
    comps2 = @components A B C
    in = Stream(:in, fs, comps1)
    out = Stream(:out, fs, comps2)

    r1 = YieldReactor(m, :r1, fs, in, out)
    fix(m[:r1_y_B_from_A], 0.3)
    fix(m[:r1_y_C_from_A], 0.4)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:r1_out_mass]) ≈ 1.0)
                (value(m[:r1_out_A_massfrac]) ≈ 0.3)
                (value(m[:r1_out_B_massfrac]) ≈ 0.3)
                (value(m[:r1_out_C_massfrac]) ≈ 0.4)
            ])
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A B
    comps2 = @components A B C D
    in = Stream(:in, fs, comps1)
    out = Stream(:out, fs, comps2)

    r1 = YieldReactor(m, :r1, fs, in, out)
    fix(m[:r1_y_C_from_A], 0.3)
    fix(m[:r1_y_D_from_A], 0.4)
    fix(m[:r1_y_C_from_B], 0.6)
    fix(m[:r1_y_D_from_B], 0.1)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:r1_out_mass]) ≈ 1.0)
                (value(m[:r1_out_A_massfrac]) ≈ 0.15)
                (value(m[:r1_out_B_massfrac]) ≈ 0.15)
                (value(m[:r1_out_C_massfrac]) ≈ 0.45)
                (value(m[:r1_out_D_massfrac]) ≈ 0.25)
            ])
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A B
    comps2 = @components C D
    in = Stream(:in, fs, comps1)
    out = Stream(:out, fs, comps2)

    r1 = YieldReactor(m, :r1, fs, in, out)
    fix(m[:r1_y_C_from_A], 0.3)
    fix(m[:r1_y_C_from_B], 0.6)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:r1_out_mass]) ≈ 1.0)
                (value(m[:r1_out_C_massfrac]) ≈ 0.45)
                (value(m[:r1_out_D_massfrac]) ≈ 0.55)
            ])
end)

@testset "YieldReactor tests" begin
    for (e, t) in assert_tests
        @test_throws e t()
    end
    for t in tests
        @test t()
    end
end
nothing