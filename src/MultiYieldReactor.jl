"""
    MultiYieldReactor <: AbstractBlock

Model a reactor with multiple inlet and outlet streams. This is used to represent a set of YieldReactors operating
in parallel. The number of inlet streams must equal the number of outlet
streams, and the length of `feed_names` must equal the number of streams.

Create a `MultiYieldReactor` block:
```julia
    MultiYieldReactor(m::Model, name::Symbol, fs::Flowsheet, inlets::Vector{Stream}, outlets::Vector{Stream},
         feed_names::Vector{Symbol}, reactor_name::Symbol)
```
or if the current scope contains `m` and `fs` bound to a `Model` and `Flowsheet`:
```
    @block(name, MultiYieldReactor, inlets, outlets, feed_names, reactor_name)
```

### Examples
```julia
julia> m = Model(); fs = Flowsheet();

julia> compsA = @components(A); compsB = @components(B); comps_out = @components A B C D;

julia> (feedA, feedB, outA, outB) = @streams begin
    feedA, compsA
    feedB, compsB
    outA , comps_out
    outB , comps_out
end;

julia> r1 = @block(r1, MultiYieldReactor, [feedA, feedB], [outA, outB], [:feedA, :feedB], :rx);

julia> print_vars(r1)
              Name               Fix      Value          Lower          Upper          Start     
--------------------------------|---|--------------|--------------|--------------|--------------|
r1_feedA_mass                     ==              1|              |              |             1|
r1_feedA_A_massfrac               ==              1|              |              |             1|
r1_feedB_mass                     ==              1|              |              |             1|
r1_feedB_B_massfrac               ==              1|              |              |             1|
r1_outA_mass                                       |              |              |              |
r1_outA_A_massfrac                                 |              |              |              |
r1_outA_B_massfrac                                 |              |              |              |
r1_outA_C_massfrac                                 |              |              |              |
r1_outA_D_massfrac                                 |              |              |              |
r1_outB_mass                                       |              |              |              |
r1_outB_A_massfrac                                 |              |              |              |
r1_outB_B_massfrac                                 |              |              |              |
r1_outB_C_massfrac                                 |              |              |              |
r1_outB_D_massfrac                                 |              |              |              |
r1_feedA_A_mass                                    |              |              |              |
r1_feedB_B_mass                                    |              |              |              |
r1_outA_A_mass                                     |              |              |              |
r1_outA_B_mass                                     |              |              |              |
r1_outA_C_mass                                     |              |              |              |
r1_outA_D_mass                                     |              |              |              |
r1_outB_A_mass                                     |              |              |              |
r1_outB_B_mass                                     |              |              |              |
r1_outB_C_mass                                     |              |              |              |
r1_outB_D_mass                                     |              |              |              |
r1_total_feed_mass                                 |              |              |              |
r1_feedA_n_rx                                      |              |              |              |
r1_feedA_rate                     ==              1|              |              |             1|
r1_feedB_n_rx                                      |              |              |              |
r1_feedB_rate                     ==              1|              |              |             1|
r1_n_rx                                            |              |              |              |
r1_feedA_y_A_from_A                                |              |              |             0|
r1_feedA_y_B_from_A               ==              0|              |              |             0|
r1_feedA_y_C_from_A               ==              0|              |              |             0|
r1_feedA_y_D_from_A               ==              0|              |              |             0|
r1_feedB_y_A_from_B               ==              0|              |              |             0|
r1_feedB_y_B_from_B                                |              |              |             0|
r1_feedB_y_C_from_B               ==              0|              |              |             0|
r1_feedB_y_D_from_B               ==              0|              |              |             0|
38 variables
```

See also [`Stream`](@ref), [`YieldReactor`](@ref), [`StoicReactor`](@ref)
"""
struct MultiYieldReactor <: AbstractBlock
    @Block_fields
    feed_names::Vector{Symbol}
    mass_flows::Dict{Symbol, Dict}
    v_nrx::Dict{Symbol, VariableRef}
    v_rates::Dict{Symbol, VariableRef}
    v_nrx_total::VariableRef
    yields::Dict{Symbol, Dict}
    total_feed_mass::VariableRef

    function MultiYieldReactor(m::Model,
                               name::Symbol,
                               fs::Flowsheet,
                               inlets::Vector{Stream},
                               outlets::Vector{Stream},
                               feed_names::Vector{Symbol},
                               reactor_name::Symbol)

        n_in = length(inlets)
        n_out = length(outlets)
        n_feeds = length(feed_names)
        @assert(n_in == n_out,
            """In MultiYieldReactor block $name, the number of inlet streams ($n_in) is not equal to
            the number of outlet streams ($n_out)""")
        @assert(n_feeds == n_in,
            """In MultiYieldReactor block $name, the number of feed names ($n_feeds) is not equal to
            the number of inlet streams ($n_in)""")
             
        @Block_init

        mass_flows = make_mass_flow_vars!(m, prefix, inlets, outlets, strm_vars, var_list, eq_list)
    
        # Make a variable for the total feed mass flow rate and an equation to calculate it.
        #    e.g., sum(rx_in_mass for all inlet streams) == rx_total_mass
        total_feed_mass = make_var!(m, prefix * "total_feed_mass", var_list)
        eq_name = prefix * "total_feed_mass_calc"
        make_eq!(m, eq_name,
            @constraint(m, 
                sum(strm_vars[s.name][:total_mass] for s in inlets) == total_feed_mass,
            base_name=eq_name), eq_list)

        # Make the n_reactor and feed_rate variables and the equations relating them.
        #   e.g., n_feed1_rx * feed1_rate == in1_mass
        v_nrx = Dict{Symbol, VariableRef}()
        v_rates = Dict{Symbol, VariableRef}()
        for (sin, fd) in zip(inlets, feed_names)
            v_nrx_i = v_nrx[fd] = make_var!(m, prefix * string(fd) * "_n_" * string(reactor_name), var_list)
            v_rates_i = v_rates[fd] = make_var!(m, prefix * string(fd) * "_rate", var_list)
            fix(v_rates_i, 1.0)

            eq_name = prefix * string(fd) * "_total_mass_calc"
            make_eq!(m, eq_name,
                @constraint(m, 
                    v_nrx_i * v_rates_i == strm_vars[sin.name][:total_mass],
                base_name=eq_name), eq_list)
        end

        # Make a variable and equation for totaling n_reactors,
        #  e.g., n_feed1_reactor + n_feed2_reactor == n_reactor
        v_nrx_total = make_var!(m, prefix * "n_" * string(reactor_name), var_list)
        eq_name = prefix * "total_" * string(reactor_name)
        make_eq!(m, eq_name,
            @constraint(m, 
                sum(v_nrx[fd] for fd in feed_names) == v_nrx_total,
            base_name=eq_name), eq_list)

        # Make the yield variables, e.g.,
        #    rx_feed1_y_cOut_from_cIn for cIn in inlet.comps, cOut in outlet.comps
        # yields[feed_name[i]][cIn][cOut] is the yield of outlet component cOut from inlet component cIn in feed_name[i].
        # All the yields are fixed at 0 except one. The free yield is yields[cIn][cIn] if cIn is in the
        #    outlet component list, otherwise the yield to the last outlet component is free.
        yields = Dict{Symbol, Dict}()
        for (sin, sout, fd) in zip(inlets, outlets, feed_names)
            yields[fd] = Dict{Symbol, OrderedDict}()
            for c_in in sin.comps
                y = yields[fd][c_in] = OrderedDict{Symbol, VariableRef}()
                for c_out in sout.comps
                    y[c_out] = make_var!(m, prefix * string(fd) * "_y_" * string(c_out) * "_from_" * string(c_in), var_list)
                    fix(y[c_out], 0.0)
                end
                c_in in keys(y) ? free(y[c_in]) : free(last_value(y))
            end
        end

        # Make the balance equations
        #   e.g., sum(rx_in1_cIn_mass * rx_feed1_y_cOut_from_cIn for cIn in inlet.comps) == rx_out1_cOut_mass
        # and equations to sum the yields to one for each inlet component, in each feed stream.
        #   e.g., sum(rx_feed1_y_cOut_from_cIn for cIn in inlet.comps) == 1.0 for cOut in outlet.comps
        for (sin, sout, fd) in zip(inlets, outlets, feed_names)
            for c_out in sout.comps
                eq_name = prefix * string(fd) * "_" * string(c_out) * "_mass_balance"
                make_eq!(m, eq_name,
                    @constraint(m, 
                        sum(mass_flows[sin.name][c_in] * yields[fd][c_in][c_out] for c_in in sin.comps) == mass_flows[sout.name][c_out],
                    base_name=eq_name), eq_list)
            end
            for c_in in sin.comps
                eq_name = prefix * string(fd) * "_" * string(c_in) * "_yield_sum"
                make_eq!(m, eq_name,
                    @constraint(m, 
                        sum(yields[fd][c_in][c_out] for c_out in sout.comps) == 1.0,
                    base_name=eq_name), eq_list)
            end
        end

        self = new(name, fs, inlets, outlets, strm_vars, var_list, eq_list, feed_names, mass_flows, v_nrx, v_rates, v_nrx_total, yields, total_feed_mass)
        @Block_finish
        return self
    end
end

function JuMP.set_start_values(blk::MultiYieldReactor; copy_inlets::Bool=true)
    copy_inlets && copy_streams(blk.inlets)
    inlets = blk.inlets
    outlets = blk.outlets
    strm_vars = blk.strm_vars
    mass_flows = blk.mass_flows
    feed_names = blk.feed_names

    # Values of the free yields and block inlet component mass flow rates.
    for (sin, fd) in zip(inlets, feed_names)
        for c in sin.comps
            set_start_values_of_free_fracs!(blk.yields[fd][c])
            is_frac(sin) && set_start_value(mass_flows[sin.name][c], get_value(strm_vars[sin.name][:fx][c]) * get_value(strm_vars[sin.name][:total_mass]))
        end
    end

    # # Values of the block outlet component mass flow rates.
    for (sin, sout, fd) in zip(inlets, outlets, feed_names)
        for c_out in sout.comps
            set_start_value(mass_flows[sout.name][c_out], sum(get_value(mass_flows[sin.name][c_in]) * 
                        get_value(blk.yields[fd][c_in][c_out]) for c_in in sin.comps)) 
        end
    end

    # Value of total_feed_mass variable.
    set_start_value(blk.total_feed_mass, sum(get_value(strm_vars[s.name][:total_mass]) for s in inlets))

    # Values of n_reactor variables.
    for (sin, fd) in zip(inlets, feed_names)
        set_start_value(blk.v_nrx[fd], get_value(strm_vars[sin.name][:total_mass]) / max(1.0e-8, get_value(blk.v_rates[fd])))
    end
    set_start_value(blk.v_nrx_total, sum(get_value(blk.v_nrx[fd]) for fd in feed_names))

    # Values of the outlet stream total mass flow rates, and mass fractions or component mass flows.
    for s in outlets
        set_start_value(strm_vars[s.name][:total_mass], sum(get_value(mass_flows[s.name][c]) for c in s.comps))
        out_is_frac = is_frac(s)
        for c in s.comps
            if out_is_frac
                set_start_value(strm_vars[s.name][:fx][c], get_value(mass_flows[s.name][c]) / max(1.0e-8, get_value(strm_vars[s.name][:total_mass])))
            else
                set_start_value(strm_vars[s.name][:fx][c], get_value(mass_flows[s.name][c]))
            end
        end
    end
end
