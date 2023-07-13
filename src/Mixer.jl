"""
    Mixer <: AbstractBlock

Mix two or more inlet streams into one outlet stream.

Create a `Mixer` block:
```julia
    Mixer(m::Model, name::Symbol, fs::Flowsheet, inlets::Vector{Stream}, outlet::Stream)
```
or if the current scope contains `m` and `fs` bound to a `Model` and `Flowsheet`:
```julia
    @block(name, Mixer, inlets, outlet)
```

### Examples

```julia
julia> m = Model(); fs = Flowsheet(); comps = @components A B;

julia> (in1, in2, out) = @streams begin
           in1, comps
           in2, comps
           out, comps
       end;

julia> mix1 = @block(mix1, Mixer, [in1, in2], out);

julia> print_vars(mix1)
              Name               Fix      Value          Lower          Upper          Start     
--------------------------------|---|--------------|--------------|--------------|--------------|
mix1_in1_mass                     ==              1|              |              |             1|
mix1_in1_A_massfrac               ==            0.5|              |              |           0.5|
mix1_in1_B_massfrac               ==            0.5|              |              |           0.5|
mix1_in2_mass                     ==              1|              |              |             1|
mix1_in2_A_massfrac               ==            0.5|              |              |           0.5|
mix1_in2_B_massfrac               ==            0.5|              |              |           0.5|
mix1_out_mass                                      |              |              |              |
mix1_out_A_massfrac                                |              |              |              |
mix1_out_B_massfrac                                |              |              |              |
9 variables
```

See also [`Stream`](@ref), [`Splitter`](@ref), [`Separator`](@ref), [`YieldReactor`](@ref), [`MultiYieldReactor`](@ref), [`StoicReactor`](@ref)
"""
struct Mixer <: AbstractBlock
    @Block_fields

    function Mixer(m::Model,
                   name::Symbol,
                   fs::Flowsheet,
                   inlets::Vector{Stream},
                   outlet::Stream)
        @assert(union([strm.comps for strm in inlets]...) == outlet.comps,
            """In Mixer $name, the outlet stream component set is not the same as
               the union of the inlet stream component sets.""")

        outlets = [outlet]
        @Block_init

        # Total mass balance, e.g., in1_mass + in2_mass == out_mass
        sout = strm_vars[outlet.name]
        eq_name = prefix * "total_mass_balance"
        make_eq!(m, eq_name,
            @constraint(m, 
                sum(strm_vars[s.name][:total_mass] for s in inlets) == sout[:total_mass],
            base_name=eq_name), eq_list)

        # Component mass balances
        #    e.g., in1_h2_mass + in2_h2_mass == out_h2_mass
        #      or  in1_h2_massfrac * in1_mass + in2_h2_massfrac * in2_mass == out_h2_massfrac * out_mass
        mass_out = (is_frac(outlet) ? sout[:total_mass] : 1.0)
        for c in outlet.comps
            eq_name = prefix * string(c) * "_mass_balance"
            make_eq!(m, eq_name,
                @constraint(m,
                    sum(strm_vars[s.name][:fx][c] * (is_frac(s) ? strm_vars[s.name][:total_mass] : 1.0) for s in inlets if c in s.comps) == 
                        sout[:fx][c] * mass_out,
                base_name=eq_name), eq_list)
        end

        self = new(name, fs, inlets, outlets, strm_vars, var_list, eq_list)
        @Block_finish
        return self
    end
end

function JuMP.set_start_values(blk::Mixer; copy_inlets::Bool=true)
    copy_inlets && copy_streams(blk.inlets)
    outlet = blk.outlets[1]
    out_vars = blk.strm_vars[outlet.name]

    # Outlet stream total mass flow rate.
    out_mass_val = max(1.0e-8, sum(get_value(blk.strm_vars[s.name][:total_mass]) for s in blk.inlets))
    set_start_value(out_vars[:total_mass], out_mass_val)

    # Outlet stream mass fractions or component mass flow rates.
    mass_out = (is_frac(outlet) ? out_mass_val : 1.0)
    for c in outlet.comps
        fx_out_val = sum(get_value(blk.strm_vars[s.name][:fx][c]) *
            (is_frac(s) ? get_value(blk.strm_vars[s.name][:total_mass]) : 1.0) for s in blk.inlets if c in s.comps) / mass_out
        set_start_value(out_vars[:fx][c], fx_out_val)
    end
    return
end
