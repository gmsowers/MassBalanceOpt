function test1()
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

    set_start_values(mix)

    @test (start_value(m[:mix_out_mass]) ≈ 500.0 &&
           start_value(m[:mix_out_A_mass]) ≈ 160.0 &&
           start_value(m[:mix_out_B_mass]) ≈ 340.0)
end

function test2()
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
    
    set_start_values(mix)
 
    @test ((start_value(m[:mix_out_mass]) ≈ 300.0) &&
           (start_value(m[:mix_out_A_massfrac]) ≈ 0.333333333) &&
           (start_value(m[:mix_out_B_massfrac]) ≈ 0.333333333) &&
           (start_value(m[:mix_out_C_massfrac]) ≈ 0.333333333)
             )
end

function test3()
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

    set_start_values(spl)
    @test ((start_value(m[:spl_out1_mass]) ≈ 30.0) &&
           (start_value(m[:spl_out1_A_massfrac]) ≈ 0.6) &&
           (start_value(m[:spl_out1_B_massfrac]) ≈ 0.4) &&
           (start_value(m[:spl_out2_mass]) ≈ 70.0) &&
           (start_value(m[:spl_out2_A_mass]) ≈ 42.0) &&
           (start_value(m[:spl_out2_B_mass]) ≈ 28.0) &&
           (start_value(m[:spl_out2_split_frac]) ≈ 0.7))
end

function test4()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()
    comps = @components A B
    in = Stream(:in, fs, comps)
    out1 = Stream(:out1, fs, comps)
    out2 = Stream(:out2, fs, comps)

    sep = Separator(m, :sep, fs, in, [out1, out2])
    fix(m[:sep_A_out1_split], 0.3)
    fix(m[:sep_B_out1_split], 0.6)

    set_start_values(sep)
    @test all([(start_value(m[:sep_out1_mass]) ≈ 0.45)
               (start_value(m[:sep_out2_mass]) ≈ 0.55)
               (start_value(m[:sep_out1_A_massfrac]) ≈ 0.33333333)
               (start_value(m[:sep_out2_B_massfrac]) ≈ 0.36363636)
               (start_value(m[:sep_A_out2_split]) ≈ 0.7)
               (start_value(m[:sep_B_out2_split]) ≈ 0.4)
            ])
end

function test5()
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

    set_start_values(r1)
    @test all([ (start_value(m[:r1_out_mass]) ≈ 1.0)
                (start_value(m[:r1_out_A_massfrac]) ≈ 0.15)
                (start_value(m[:r1_out_B_massfrac]) ≈ 0.15)
                (start_value(m[:r1_out_C_massfrac]) ≈ 0.45)
                (start_value(m[:r1_out_D_massfrac]) ≈ 0.25)
            ])
end

function test6()
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

    set_start_values(r1)
    @test all([(start_value(m[:r1_total_feed_mass]) ≈ 150.0)
               (start_value(m[:r1_cgeth_C_massfrac]) ≈ 0.55)
               (start_value(m[:r1_cgpro_C_massfrac]) ≈ 0.4)
               (start_value(m[:r1_eth_n_furn]) ≈ 4.0)
               (start_value(m[:r1_pro_n_furn]) ≈ 2.0)
            ])
end

function test7()
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

    set_start_values(arx)
    @test all([ ≈(start_value(m[:arx_out1_C2H2_massfrac]), 0.0, atol=1.0e-8)
                ≈(start_value(m[:arx_out1_H2_moles]), 0.1525516013, rtol=1.0e-8)
                 (start_value(m[:arx_out1_C2H4_moles]) ≈ 0.012802156497)
                ≈(start_value(m[:arx_out1_N2_mass]), 0.33333333333, rtol=1.0e-8)
            ])
end

function test8()
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

    @values begin
        mix_in1_mass=300.0
        mix_in1_A_mass=300.0
        mix_in2_mass=150.0
        mix_in2_B_mass=150.0
        spl_out1_split_frac=0.3
    end

    connect(fs)
    set_start_values([mix, spl, mix2])

    @test all([ (start_value(m[:mix_mixed_mass]) ≈ 450.0)
                (start_value(m[:spl_out1_mass]) ≈ 135.0)
                (start_value(m[:mix2_mixed2_A_mass]) ≈ start_value(m[:mix_in1_A_mass]))
                (start_value(m[:mix2_mixed2_B_mass]) ≈ start_value(m[:mix_in2_B_mass]))
            ])
end

@testset "set_start_values tests" begin
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
    test7()
    test8()
end
nothing
