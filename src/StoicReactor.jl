"""
    StoicReactor <: AbstractBlock

Model a reactor in which a set of stoichiometric reactions with specified conversions takes place.

Create a `StoicReactor` block:
```julia
    StoicReactor(m::Model, name::Symbol, fs::Flowsheet, inlet::Stream, outlet::Stream,
        stoic_coef::Vector{OrderedDict{Symbol, Real}}, 
        mw::Dict{Symbol, Float64},
        conv::OrderedDict{Int, NamedTuple{(:c, :X), Tuple{Symbol, Float64}}}= OrderedDict{Int, NamedTuple{(:c, :X),
                               Tuple{Symbol, Float64}}}())
```
or if the current scope contains `m` and `fs` bound to a `Model` and `Flowsheet`:
```
    @block(name, StoicReactor, inlet, outlet, stoic_coef, mw, conv)
```

### Examples
```julia
julia> m = Model(); fs = Flowsheet(); comps1 = @components A B; comps2 = @components A B C D;

julia> (in1, out1) = @streams begin
    in1, comps1
    out1, comps2
end;

julia> mw = Dict(:A => 30.0, :B => 28.0, :C => 35.0, :D => 30.0);

julia> coef = @stoic A + B => C + D
1-element Vector{OrderedDict{Symbol, Real}}:
 OrderedDict(:A => -1, :B => -1, :C => 1, :D => 1)

 julia> conv = OrderedDict(1 => (c=:A, X=0.65))
 OrderedDict{Int64, NamedTuple{(:c, :X), Tuple{Symbol, Float64}}} with 1 entry:
   1 => (c = :A, X = 0.65)
   
julia> r1 = @block(r1, StoicReactor, in1, out1, coef, mw, conv);

julia> print_vars(r1)
              Name               Fix      Value          Lower          Upper          Start     
--------------------------------|---|--------------|--------------|--------------|--------------|
r1_in1_mass                       ==              1|              |              |             1|
r1_in1_A_massfrac                 ==            0.5|              |              |           0.5|
r1_in1_B_massfrac                 ==            0.5|              |              |           0.5|
r1_out1_mass                                       |              |              |              |
r1_out1_A_massfrac                                 |              |              |              |
r1_out1_B_massfrac                                 |              |              |              |
r1_out1_C_massfrac                                 |              |              |              |
r1_out1_D_massfrac                                 |              |              |              |
r1_in1_A_mass                                      |              |              |              |
r1_in1_B_mass                                      |              |              |              |
r1_out1_A_mass                                     |              |              |              |
r1_out1_B_mass                                     |              |              |              |
r1_out1_C_mass                                     |              |              |              |
r1_out1_D_mass                                     |              |              |              |
r1_in1_A_moles                                     |              |              |              |
r1_in1_B_moles                                     |              |              |              |
r1_out1_A_moles                                    |              |              |              |
r1_out1_B_moles                                    |              |              |              |
r1_out1_C_moles                                    |              |              |              |
r1_out1_D_moles                                    |              |              |              |
r1_extent_rx_1                                     |              |              |              |
r1_conv_A_rx_1                    ==           0.65|              |              |          0.65|
22 variables
```

See also [`@stoic`](@ref), [`Stream`](@ref), [`YieldReactor`](@ref), [`MultiYieldReactor`](@ref)
"""
struct StoicReactor <: AbstractBlock
    @Block_fields
    mass_flows::Dict{Symbol, Dict}
    mole_flows::Dict{Symbol, OrderedDict}
    mw::Dict{Symbol, Float64}
    stoic_coef::Vector{OrderedDict{Symbol, Real}}
    rx_comps::OrderedSet
    inert_comps::OrderedSet
    extents::Vector{VariableRef}
    conversions::OrderedDict
    conversion_vars::Vector{VariableRef}

    function StoicReactor(m::Model,
                          name::Symbol,
                          fs::Flowsheet,
                          inlet::Stream,
                          outlet::Stream,
                          stoic_coef::Vector{OrderedDict{Symbol, Real}},
                          mw::Dict{Symbol, Float64},
                          conversions::OrderedDict{Int, NamedTuple{(:c, :X), Tuple{Symbol, Float64}}}=
                                       OrderedDict{Int, NamedTuple{(:c, :X), Tuple{Symbol, Float64}}}())

        inlets = [inlet]
        outlets = [outlet]
        @Block_init

        mass_flows = make_mass_flow_vars!(m, prefix, inlets, outlets, strm_vars, var_list, eq_list)
        mole_flows = Dict{Symbol, OrderedDict}()

        rx_comps = OrderedSet(c for sc in stoic_coef for c in keys(sc))
        inert_comps = OrderedSet(c for c in outlet.comps if !(c in rx_comps))

        # Equations to calculate the molar flow rates of the reacting components, in both
        #    inlet and outlet streams.
        #    e.g., arx_in_c2h2_mass == arx_in_c2h2_moles * mw[:c2h2]
        for s in (inlet, outlet)
            s_mole_flows = mole_flows[s.name] = OrderedDict{Symbol, VariableRef}()
            prefix_ = prefix * string(s.name) * '_'
            for c in s.comps
                if c in rx_comps
                    s_mole_flows[c] = make_var!(m, prefix_ * string(c) * "_moles", var_list)

                    eq_name = prefix_ * string(c) * "_moles_calc"
                    make_eq!(m, eq_name,
                        @constraint(m, 
                            mass_flows[s.name][c] == s_mole_flows[c] * mw[c],
                        base_name=eq_name), eq_list)
                end
            end
        end

        # Equations relating the extents of reaction and the component molar flow rates.
        #    e.g., arx_out_h2_moles == arx_in_h2_moles + sum(stoic_coef[i][:h2] * arx_extent_rx_i)
        rx_idxs = collect(eachindex(stoic_coef))
        extents = [make_var!(m, prefix * "extent_rx_" * string(i), var_list) for i in rx_idxs]
        for c in rx_comps
            inlet_mole_flow = get(mole_flows[inlet.name], c, 0.0)
            outlet_mole_flow = get(mole_flows[outlet.name], c, 0.0)
            eq_name = prefix * string(c) * "_balance"
            make_eq!(m, eq_name,
                @constraint(m, 
                    outlet_mole_flow == inlet_mole_flow + sum(extents[i] * get(stoic_coef[i], c, 0.0) for i in rx_idxs),
                base_name=eq_name), eq_list)
        end

        # Equations relating the extents of reaction and the conversions.
        #    e.g., arx_in_c2h2_moles * arx_conv_c2h2_rx_1 == -stoic_coef[1][:c2h2] * arx.extent_rx_1
        conversion_vars = Vector{VariableRef}(undef, length(conversions))
        for (i, cx) in conversions
            i_str = string(i)
            c_str = string(cx.c)
            conversion_vars[i] = make_var!(m, prefix * "conv_" * c_str * "_rx_" * i_str, var_list)
            fix(conversion_vars[i], cx.X)
            eq_name = prefix * "conv_" * c_str * "_rx_" * i_str * "_def"
            make_eq!(m, eq_name,
                @constraint(m, 
                    mole_flows[inlet.name][cx.c] * conversion_vars[i] == -stoic_coef[i][cx.c] * extents[i],
                base_name=eq_name), eq_list)
        end

        # For the inert components, copy the inlet component mass flow to the outlet.
        #    e.g., arx_out_co2_mass == arx_in_co2_mass
        for c in inert_comps
            eq_name = prefix * string(c) * "_balance"
            make_eq!(m, eq_name,
                @constraint(m, 
                    mass_flows[outlet.name][c] == mass_flows[inlet.name][c],
                base_name=eq_name), eq_list)
        end

        self = new(name, fs, inlets, outlets, strm_vars, var_list, eq_list, mass_flows, mole_flows, mw, stoic_coef,
                     rx_comps, inert_comps, extents, conversions, conversion_vars)
        @Block_finish
        return self
    end
end

function JuMP.set_start_values(blk::StoicReactor; copy_inlets::Bool=true)
    copy_inlets && copy_streams(blk.inlets)
    inlet = blk.inlets[1]
    outlet = blk.outlets[1]
    in_strm_vars = blk.strm_vars[inlet.name]
    in_mass_flows = blk.mass_flows[inlet.name]
    out_mass_flows = blk.mass_flows[outlet.name]
    in_mole_flows = blk.mole_flows[inlet.name]
    out_mole_flows = blk.mole_flows[outlet.name]
    out_strm_vars = blk.strm_vars[outlet.name]
    
    # Values of block inlet component mass flow rates and component molar flow rates of the inlet reacting components.
    for c in inlet.comps
        is_frac(inlet) && set_start_value(in_mass_flows[c], get_value(in_strm_vars[:fx][c]) * get_value(in_strm_vars[:total_mass]))
        (c in blk.rx_comps) && set_start_value(in_mole_flows[c], get_value(in_mass_flows[c]) / blk.mw[c])
    end

    # Values of the outlet component mass flow rates of the inert components.
    for c in blk.inert_comps
        set_start_value(out_mass_flows[c], get_value(in_mass_flows[c]))
    end

    # Values of the extents.
    for (i, cx) in blk.conversions
        set_start_value(blk.extents[i], -get_value(in_mole_flows[cx.c]) * get_value(blk.conversion_vars[i]) / blk.stoic_coef[i][cx.c])
    end

    # Values of the outlet molar and mass flow rates.
    rx_idxs = collect(eachindex(blk.stoic_coef))
    inlet_comps = keys(in_mole_flows)
    outlet_comps = keys(out_mole_flows)
    for c in blk.rx_comps
        if c in outlet_comps
            inlet_mole_flow = (c in inlet_comps ? get_value(in_mole_flows[c]) : 0.0)
            outlet_mole_flow = inlet_mole_flow + sum(get_value(blk.extents[i]) * get(blk.stoic_coef[i], c, 0.0) for i in rx_idxs)
            set_start_value(out_mole_flows[c], outlet_mole_flow)
            set_start_value(out_mass_flows[c], outlet_mole_flow * blk.mw[c])
        end
    end

    # Values of outlet total mass flow rate and outlet mass fractions or component mass flows.
    set_start_value(out_strm_vars[:total_mass], sum(get_value(out_mass_flows[c]) for c in outlet.comps))
    out_is_frac = is_frac(outlet)
    for c in outlet.comps
        if out_is_frac
            set_start_value(out_strm_vars[:fx][c], get_value(out_mass_flows[c]) / max(1.0e-8, get_value(out_strm_vars[:total_mass])))
        else
            set_start_value(out_strm_vars[:fx][c], get_value(out_mass_flows[c]))
        end
    end
end
