include("model_init.jl")

function model_test()
    m = Model(Ipopt.Optimizer)
    fs = Flowsheet()

    compsCG      = @components  h2 ch4 c2h2 c2h4 c2h6 c3h6 c3h8 mapd btd c4s benz tx bzraff foil
    compsEthFeed = @components  c2h6
    compsEthRec  = @components  c2h4 c2h6
    compsEthFurn = @components  c2h4 c2h6
    compsProFeed = @components  c3h8
    compsProRec  = @components  c3h6 c3h8 mapd
    compsProFurn = @components  c3h6 c3h8 mapd
    compsHpygas  = @components  benz tx bzraff foil
    compsCGCfeed = @components  h2 ch4 c2h2 c2h4 c2h6 c3h6 c3h8 mapd btd c4s benz tx bzraff
    compsLpygas  = @components  benz tx bzraff
    compsDC2     = @components  h2 ch4 c2h2 c2h4 c2h6 c3h6 c3h8 mapd btd c4s
    compsDC2oh   = @components  h2 ch4 c2h2 c2h4 c2h6
    compsDC2bt   = @components  c3h6 c3h8 mapd btd c4s
    compsARX     = @components  h2 ch4 c2h4 c2h6
    compsOG      = @components  h2 ch4 c2h4
    compsC2S     = @components  c2h4 c2h6
    compsDC3oh   = @components  c3h6 c3h8 mapd
    compsDC3bt   = @components  btd c4s
    compsC3Soh   = @components  c3h6 c3h8

    # Ethane feed header.
    (ethfd, ethre, ethf1f5) = @streams begin
        ethfd   , compsEthFeed
        ethre   , compsEthRec
        ethf1f5 , compsEthFurn
    end
    @block(ethhdr, Mixer, [ethfd, ethre], ethf1f5)

    # Propane feed header.
    (profd, prore, prof1f5) = @streams begin
        profd   , compsProFeed
        prore   , compsProRec
        prof1f5 , compsProFurn
    end
    @block(prohdr, Mixer, [profd, prore], prof1f5)

    # Furnaces, F1--F5.
    (cgethf1f5, cgprof1f5) = @streams begin
        cgethf1f5 , compsCG, (basis=FLOW)
        cgprof1f5 , compsCG, (basis=FLOW)
    end
    @block(f1f5, MultiYieldReactor, [ethf1f5, prof1f5], [cgethf1f5, cgprof1f5], [:eth, :pro], :furn)
    connect([ethf1f5, prof1f5])

    # Cracked gas header.
    cgas = @stream(cgas, compsCG)
    @block(cghdr, Mixer, [cgethf1f5, cgprof1f5], cgas)
    connect([cgethf1f5, cgprof1f5])

    # Quench section.
    (hpygas, cgcfeed) = @streams begin
        hpygas  , compsHpygas
        cgcfeed , compsCGCfeed
    end
    @block(quench, Separator, cgas, [hpygas, cgcfeed])
    connect(cgas)

    # Cracked gas compressor section.
    (lpygas, cgcout) = @streams begin
        lpygas , compsLpygas
        cgcout , compsDC2
    end
    @block(cgc, Separator, cgcfeed, [lpygas, cgcout])
    connect(cgcfeed)

    # Front-end deethanizer.
    (dc2oh, dc2bt) = @streams begin
        dc2oh , compsDC2oh
        dc2bt , compsDC2bt
    end
    @block(dc2, Separator, cgcout, [dc2oh, dc2bt])
    connect(cgcout)

    # Acetylene reactors.
    arxo = @stream(arxo, compsARX)
    arx_stoic = @stoic begin
        c2h2 + h2  => c2h4
        c2h2 + 2h2 => c2h6
    end
    arx_mw = Dict(:c2h2 => 26.03728,
                :h2   =>  2.01588,
                :c2h4 => 28.05316,
                :c2h6 => 30.06904)
    arx_conv = OrderedDict(1 => (c=:c2h2, X=0.7))
    @block(arx, StoicReactor, dc2oh, arxo, arx_stoic, arx_mw, arx_conv)
    connect(dc2oh)

    # Cold train/demethanizer.
    (offgas, dc1bt) = @streams begin
        offgas , compsOG
        dc1bt  , compsC2S
    end
    @block(ctdc1, Separator, arxo, [offgas, dc1bt])
    @specs(ctdc1_c2h4_offgas_split ~ ctdc1_offgas_c2h4_massfrac)
    connect(arxo)

    # C2 splitter.
    c2h4p = @stream(c2h4p, compsC2S)
    @block(c2s, Separator, dc1bt, [c2h4p, ethre])
    @specs begin
        c2s_c2h4p_c2h6_massfrac ~ c2s_c2h4_c2h4p_split
        c2s_ethre_c2h4_massfrac ~ c2s_c2h6_c2h4p_split
    end
    connect([dc1bt, ethre])

    # Depropanizer.
    (dc3oh, c4sp) = @streams begin
        dc3oh  , compsDC3oh
        c4sp   , compsDC3bt
    end
    @block(dc3, Separator, dc2bt, [dc3oh, c4sp])
    connect(dc2bt)

    # C3 splitter.
    c3h6p = @stream(c3h6p, compsC3Soh)
    @block(c3s, Separator, dc3oh, [c3h6p, prore])
    @specs begin
        c3s_prore_c3h6_massfrac ~ c3s_c3h6_prore_split
        c3s_c3h6p_c3h8_massfrac ~ c3s_c3h8_prore_split
    end
    connect([dc3oh, prore])

    # Calculate the plant mass balance.
    @variable(m, plant_mass_imbalance)
    @constraint(m, plant_mass_imbalance ==
        m[:ethhdr_ethfd_mass]  +
        m[:prohdr_profd_mass]  -
        m[:quench_hpygas_mass] -
        m[:cgc_lpygas_mass]    -
        m[:ctdc1_offgas_mass]  -
        m[:c2s_c2h4p_mass]     -
        m[:dc3_c4sp_mass]      -
        m[:c3s_c3h6p_mass], base_name="plant_mass_balance")

    # Start values for the variables.
    initialize(m)

    # Prices and objective function.
    prices = Dict(
                :ethane_feed  => 10.0,
                :propane_feed => 15.0,
                :c2h4_prod    => 40.0,
                :c3h6_prod    => 30.0,
                :c4s_prod     => 28.0,
                :pygas_prod   => 23.0,
                :offgas_prod  => 25.0
                )
    map!(p -> p/100.0, values(prices))

    costs = @expression(m, m[:ethhdr_ethfd_mass]  * prices[:ethane_feed] +
                           m[:prohdr_profd_mass]  * prices[:propane_feed])
    sales = @expression(m, m[:c2s_c2h4p_mass]     * prices[:c2h4_prod] +
                           m[:c3s_c3h6p_mass]     * prices[:c3h6_prod] +
                           m[:dc3_c4sp_mass]      * prices[:c4s_prod]  +
                           m[:quench_hpygas_mass] * prices[:pygas_prod] +
                           m[:cgc_lpygas_mass]    * prices[:pygas_prod] +
                           m[:ctdc1_offgas_mass]  * prices[:offgas_prod])
    @objective(m, Max, sales - costs)

    # Make the ethane and propane feed flow rates degrees of freedom.
    @specs begin
        -ethhdr_ethfd_mass
        -prohdr_profd_mass
    end

    # Some variable bounds to make the problem feasible.
    @bounds begin
        1.0 < ethhdr_ethfd_mass < 2.0e5
        1.0 < prohdr_profd_mass < 2.0e5
        f1f5_n_furn < 5.0
        0.0 < f1f5_pro_n_furn < 3.0
        1.0 < f1f5_eth_n_furn
    end

    set_silent(m)
    @solve

    return (value(m[:f1f5_eth_n_furn]) ≈ 2.0 && value(m[:f1f5_pro_n_furn]) ≈ 3.0 && value(m[:f1f5_n_furn]) ≈ 5.0)
end # model_test

@testset "Model test" begin
    @test model_test()
end
nothing