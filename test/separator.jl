tests = []
assert_tests = []

push!(assert_tests, (AssertionError, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A
    comps2 = @components B
    in = Stream(:in, fs, comps1)
    out = Stream(:out, fs, comps2)

    sep = Separator(m, :sep, fs, in, [out])
    return nothing
end))

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps1 = @components A B
    comps2 = @components A
    comps3 = @components B
    in = Stream(:in, fs, comps1)
    out1 = Stream(:out1, fs, comps2)
    out2 = Stream(:out2, fs, comps3)

    sep = Separator(m, :sep, fs, in, [out1, out2])
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:sep_out1_mass]) ≈ 0.5)
                (value(m[:sep_out2_mass]) ≈ 0.5)
                (value(m[:sep_out1_A_massfrac]) ≈ 1.0)
                (value(m[:sep_out2_B_massfrac]) ≈ 1.0)
                (value(m[:sep_A_out1_split]) ≈ 1.0)
                (value(m[:sep_B_out2_split]) ≈ 1.0)
            ])
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A B
    in = Stream(:in, fs, comps)
    out1 = Stream(:out1, fs, comps)
    out2 = Stream(:out2, fs, comps)

    sep = Separator(m, :sep, fs, in, [out1, out2])
    fix(m[:sep_A_out1_split], 0.3)
    fix(m[:sep_B_out1_split], 0.6)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:sep_out1_mass]) ≈ 0.45)
                (value(m[:sep_out2_mass]) ≈ 0.55)
                (value(m[:sep_out1_A_massfrac]) ≈ 0.33333333)
                (value(m[:sep_out2_B_massfrac]) ≈ 0.36363636)
                (value(m[:sep_A_out2_split]) ≈ 0.7)
                (value(m[:sep_B_out2_split]) ≈ 0.4)
            ])
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A B C
    comps2 = @components B C

    in = Stream(:in, fs, comps)
    out1 = Stream(:out1, fs, comps, basis=FLOW)
    out2 = Stream(:out2, fs, comps2, basis=FLOW)

    sep = Separator(m, :sep, fs, in, [out1, out2])
    fix(m[:sep_B_out1_split], 0.6)
    fix(m[:sep_C_out1_split], 0.3)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:sep_out1_mass]) ≈ 0.633333333)
                (value(m[:sep_out2_mass]) ≈ 0.366666667)
                (value(m[:sep_out1_A_mass]) ≈ 0.33333333)
                (value(m[:sep_out1_B_mass]) ≈ 0.2)
                (value(m[:sep_out1_C_mass]) ≈ 0.1)
                (value(m[:sep_out2_B_mass]) ≈ 0.133333333333)
                (value(m[:sep_out2_C_mass]) ≈ 0.23333333)
                (value(m[:sep_B_out2_split]) ≈ 0.4)
                (value(m[:sep_C_out2_split]) ≈ 0.7)
            ])
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A B C

    in = Stream(:in, fs, comps)
    out1 = Stream(:out1, fs, comps)
    out2 = Stream(:out2, fs, comps)

    sep = Separator(m, :sep, fs, in, [out1, out2])
    fix(m[:sep_A_out1_split], 0.0)
    fix(m[:sep_B_out1_split], 0.6)
    fix(m[:sep_C_out1_split], 0.3)
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:sep_out1_mass]) ≈ 0.3)
                (value(m[:sep_out2_mass]) ≈ 0.7)
                isapprox(value(m[:sep_out1_A_massfrac]), 0.0; atol=1.0e-8)
                (value(m[:sep_out1_B_massfrac]) ≈ 0.6666666666667)
                (value(m[:sep_out1_C_massfrac]) ≈ 0.3333333333333)
                (value(m[:sep_out2_A_massfrac]) ≈ 0.47619048)
                (value(m[:sep_out2_B_massfrac]) ≈ 0.19047619)
                (value(m[:sep_out2_C_massfrac]) ≈ 0.3333333333333)
                (value(m[:sep_A_out2_split]) ≈ 1.0)
                (value(m[:sep_B_out2_split]) ≈ 0.4)
                (value(m[:sep_C_out2_split]) ≈ 0.7)
            ])
end)

push!(tests, 
function ()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A B C D
    comps1 = @components A B
    comps2 = @components C D

    in = Stream(:in, fs, comps)
    out1 = Stream(:out1, fs, comps1, basis=FLOW)
    out2 = Stream(:out2, fs, comps2, basis=FLOW)

    sep = Separator(m, :sep, fs, in, [out1, out2])
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:sep_out1_mass]) ≈ 0.5)
                (value(m[:sep_out2_mass]) ≈ 0.5)
                (value(m[:sep_out1_A_mass]) ≈ 0.25)
                (value(m[:sep_out1_B_mass]) ≈ 0.25)
                (value(m[:sep_out2_C_mass]) ≈ 0.25)
                (value(m[:sep_out2_D_mass]) ≈ 0.25)
            ])
end)

@testset "Separator tests" begin
    for (e, t) in assert_tests
        @test_throws e t()
    end
    for t in tests
        @test t()
    end
end
nothing