"""
    YieldReactor <: AbstractBlock

Model a reactor with one inlet and one outlet stream, in which each component in the inlet stream has a specified
yield to the components in the outlet stream.

Create a `YieldReactor` block:
```julia
    YieldReactor(m::Model, name::Symbol, fs::Flowsheet, inlet::Stream, outlet::Stream)
```julia
or if the current scope contains `m` and `fs` bound to a `Model` and `Flowsheet`:
```
    @block(name, YieldReactor, inlet, outlet)
```

### Examples
```julia
julia> m = Model(); fs = Flowsheet(); comps1 = @components A B; comps2 = @components A B C D;

julia> (in1, out1) = @streams begin
    in1, comps1
    out1, comps2
end;

julia> r1 = @block(r1, YieldReactor, in1, out1);

julia> @set begin
    r1_y_C_from_A = 0.3
    r1_y_D_from_A = 0.4
    r1_y_C_from_B = 0.6
    r1_y_D_from_B = 0.1
end

julia> set_start_values(r1)

julia> print_vars(r1)
Name                             Fix      Value          Lower          Upper          Start     
--------------------------------|---|--------------|--------------|--------------|--------------|
r1_in1_mass                       ==              1|              |              |             1|
r1_in1_A_massfrac                 ==            0.5|              |              |           0.5|
r1_in1_B_massfrac                 ==            0.5|              |              |           0.5|
r1_out1_mass                                       |              |              |             1|
r1_out1_A_massfrac                                 |              |              |          0.15|
r1_out1_B_massfrac                                 |              |              |          0.15|
r1_out1_C_massfrac                                 |              |              |          0.45|
r1_out1_D_massfrac                                 |              |              |          0.25|
r1_in1_A_mass                                      |              |              |           0.5|
r1_in1_B_mass                                      |              |              |           0.5|
r1_out1_A_mass                                     |              |              |          0.15|
r1_out1_B_mass                                     |              |              |          0.15|
r1_out1_C_mass                                     |              |              |          0.45|
r1_out1_D_mass                                     |              |              |          0.25|
r1_y_A_from_A                                      |              |              |           0.3|
r1_y_B_from_A                     ==              0|              |              |             0|
r1_y_C_from_A                     ==            0.3|              |              |           0.3|
r1_y_D_from_A                     ==            0.4|              |              |           0.4|
r1_y_A_from_B                     ==              0|              |              |             0|
r1_y_B_from_B                                      |              |              |           0.3|
r1_y_C_from_B                     ==            0.6|              |              |           0.6|
r1_y_D_from_B                     ==            0.1|              |              |           0.1|
22 variables
```

See also [`Stream`](@ref), [`MultiYieldReactor`](@ref), [`StoicReactor`](@ref)
"""
struct YieldReactor <: AbstractBlock
    @Block_fields
    mass_flows::Dict{Symbol, Dict}
    yields::Dict{Symbol, OrderedDict}

    function YieldReactor(m::Model,
                          name::Symbol,
                          fs::Flowsheet,
                          inlet::Stream,
                          outlet::Stream)

        inlets = [inlet]
        outlets = [outlet]
        @Block_init

        mass_flows = make_mass_flow_vars!(m, prefix, inlets, outlets, strm_vars, var_list, eq_list)

        # Make the yield variables, e.g.,
        #    rx_y_cOut_from_cIn for cIn in inlet.comps, cOut in outlet.comps
        # yields[cIn][cOut] is the yield of outlet component cOut from inlet component cIn.
        # All the yields are fixed at 0 except one. The free yield is yields[cIn][cIn] if cIn is in the
        #    outlet component list, otherwise the yield to the last outlet component is free.
        yields = Dict{Symbol, OrderedDict}()
        for c_in in inlet.comps
            y = yields[c_in] = OrderedDict{Symbol, VariableRef}()
            for c_out in outlet.comps
                y[c_out] = make_var!(m, prefix * "y_" * string(c_out) * "_from_" * string(c_in), var_list)
                fix(y[c_out], 0.0)
            end
            c_in in keys(y) ? free(y[c_in]) : free(last_value(y))
        end
        
        # Make the balance equations
        #   e.g., sum(rx_in_cIn_mass * rx_y_cOut_from_cIn for cIn in inlet.comps) == rx_out_cOut_mass
        # and equations to sum the yields to one for each inlet component.
        #   e.g., sum(rx_y_cOut_from_cIn for cIn in inlet.comps) == 1.0
        for c_out in outlet.comps
            eq_name = prefix * string(c_out) * "_mass_balance"
            make_eq!(m, eq_name,
                @constraint(m, 
                    sum(mass_flows[inlet.name][c_in] * yields[c_in][c_out] for c_in in inlet.comps) == mass_flows[outlet.name][c_out],
                base_name=eq_name), eq_list)
        end
        for c_in in inlet.comps
            eq_name = prefix * string(c_in) * "_yield_sum"
            make_eq!(m, eq_name,
                @constraint(m, 
                    sum(yields[c_in][c_out] for c_out in outlet.comps) == 1.0,
                base_name=eq_name), eq_list)
        end

        self = new(name, fs, inlets, outlets, strm_vars, var_list, eq_list, mass_flows, yields)
        @Block_finish
        return self
    end
end

function JuMP.set_start_values(blk::YieldReactor; copy_inlets::Bool=true)
    copy_inlets && copy_streams(blk.inlets)
    inlet = blk.inlets[1]
    outlet = blk.outlets[1]
    in_vars = blk.strm_vars[inlet.name]
    in_mass_flows = blk.mass_flows[inlet.name]
    out_vars = blk.strm_vars[outlet.name]
    out_mass_flows = blk.mass_flows[outlet.name]

    # Values of the free yields and the block inlet component mass flow rates.
    for c in inlet.comps
        set_start_values_of_free_fracs!(blk.yields[c])
        is_frac(inlet) && set_start_value(in_mass_flows[c], get_value(in_vars[:fx][c]) * get_value(in_vars[:total_mass]))
    end

    # Values of the block outlet component mass flow rates.
    for c_out in outlet.comps
        set_start_value(out_mass_flows[c_out], sum(get_value(in_mass_flows[c_in]) * get_value(blk.yields[c_in][c_out]) for c_in in inlet.comps))
    end

    # Values of the outlet stream total mass flow rates, and mass fractions or component mass flows.
    set_start_value(out_vars[:total_mass], sum(get_value(out_mass_flows[c]) for c in outlet.comps))
    out_is_frac = is_frac(outlet)
    for c in outlet.comps
        if out_is_frac
            set_start_value(out_vars[:fx][c], get_value(out_mass_flows[c]) / max(1.0e-8, get_value(out_vars[:total_mass])))
        else
            set_start_value(out_vars[:fx][c], get_value(out_mass_flows[c]))
        end
    end
end
