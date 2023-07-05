function test1()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    compsA = @components A
    compsB = @components B
    comps = @components A B

    in1 = Stream(:in1, fs, compsA)
    in2 = Stream(:in2, fs, compsB)
    mixed = Stream(:mixed, fs, comps)
    mix = Mixer(m, :mix, fs, [in1, in2], mixed)

    out1 = Stream(:out1, fs, comps)
    out2 = Stream(:out2, fs, comps)
    spl = Splitter(m, :spl, fs, mixed, [out1, out2])

    connect(fs)

    @set mix_in1_mass=100.0
    @set mix_in2_mass=200.0
    @set spl_out1_split_frac=0.3

    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:mix_mixed_mass]) ≈ 300.0)
                (value(m[:spl_out1_mass]) ≈ 90.0)
                (value(m[:spl_out1_B_massfrac]) ≈ 0.66666666667)
            ])
end

function test2()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A B

    in1 = Stream(:in1, fs, comps)
    out1 = Stream(:out1, fs, comps)
    out2 = Stream(:out2, fs, comps)
    spl = Splitter(m, :spl, fs, in1, [out1, out2])
    
    mixed = Stream(:mixed, fs, comps)
    mix = Mixer(m, :mix, fs, [out1, out2], mixed)

    connect(fs)

    @values begin
        spl_in1_mass        = 100.0
        spl_in1_A_massfrac  = 0.25
        spl_in1_B_massfrac  = 0.75
        spl_out1_split_frac = 0.7
    end

    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:mix_mixed_mass]) ≈ 100.0)
                (value(m[:spl_out1_mass]) ≈ 70.0)
                (value(m[:mix_mixed_A_massfrac]) ≈ 0.25)
            ])
end

function test3()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    compsA = @components A
    compsB = @components B
    comps = @components A B

    in1 = Stream(:in1, fs, compsA)
    in2 = Stream(:in2, fs, compsB)
    mixed = Stream(:mixed, fs, comps)
    mix = Mixer(m, :mix, fs, [in1, in2], mixed)

    out1 = Stream(:out1, fs, comps)
    out2 = Stream(:out2, fs, comps)
    spl = Splitter(m, :spl, fs, mixed, [out1, out2])

    mixed2 = Stream(:mixed2, fs, comps)
    mix2 = Mixer(m, :mix2, fs, [out1, out2], mixed2)

    connect(fs)

    @set mix_in1_mass=100.0
    @set mix_in2_mass=200.0
    @set spl_out1_split_frac=0.3

    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:mix_mixed_mass]) ≈ 300.0)
                (value(m[:spl_out1_mass]) ≈ 90.0)
                (value(m[:mix2_mixed2_A_massfrac]) ≈ value(m[:mix_mixed_A_massfrac]))
            ])
end

function test4()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    compsA = @components A
    compsB = @components B
    comps = @components A B

    in1 = Stream(:in1, fs, compsA, basis=FLOW)
    in2 = Stream(:in2, fs, compsB, basis=FLOW)
    mixed = Stream(:mixed, fs, comps)
    mix = Mixer(m, :mix, fs, [in1, in2], mixed)

    out1 = Stream(:out1, fs, comps)
    out2 = Stream(:out2, fs, comps)
    spl = Splitter(m, :spl, fs, mixed, [out1, out2])

    mixed2 = Stream(:mixed2, fs, comps, basis=FLOW)
    mix2 = Mixer(m, :mix2, fs, [out1, out2], mixed2)

    connect(fs)

    @values begin
        mix_in1_mass=300.0
        mix_in1_A_mass=300.0
        mix_in2_mass=150.0
        mix_in2_B_mass=150.0
        spl_out1_split_frac=0.3
    end
    
    set_silent(m)
    optimize!(m)
    return all([ (termination_status(m) == LOCALLY_SOLVED)
                (value(m[:mix_mixed_mass]) ≈ 450.0)
                (value(m[:spl_out1_mass]) ≈ 135.0)
                (value(m[:mix2_mixed2_A_mass]) ≈ value(m[:mix_in1_A_mass]))
                (value(m[:mix2_mixed2_B_mass]) ≈ value(m[:mix_in2_B_mass]))
            ])
end

@testset "Mixer+Splitter tests" begin
    @test test1()
    @test test2()
    @test test3()
    @test test4()
end
nothing