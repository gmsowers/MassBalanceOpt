"""
    Splitter <: AbstractBlock

Split one inlet stream into two or more outlet streams.

Create a `Splitter` block:
```julia
    Splitter(m::Model, name::Symbol, fs::Flowsheet, inlet::Stream, outlets::Vector{Stream})
```
or if the current scope contains `m` and `fs` bound to a `Model` and `Flowsheet`:
```julia
    @block(name, Splitter, inlet, outlets)
```

### Examples
    
```julia
julia> m = Model(); fs = Flowsheet(); comps = @components A B;

julia> (in1, out1, out2) = @streams begin
           in1, comps
           out1, comps
           out2, comps
       end;

julia> spl = @block(spl, Splitter, in1, [out1, out2]);

julia> print_vars(spl)
              Name               Fix      Value          Lower          Upper          Start     
--------------------------------|---|--------------|--------------|--------------|--------------|
spl_in1_mass                      ==              1|              |              |             1|
spl_in1_A_massfrac                ==            0.5|              |              |           0.5|
spl_in1_B_massfrac                ==            0.5|              |              |           0.5|
spl_out1_mass                                      |              |              |              |
spl_out1_A_massfrac                                |              |              |              |
spl_out1_B_massfrac                                |              |              |              |
spl_out2_mass                                      |              |              |              |
spl_out2_A_massfrac                                |              |              |              |
spl_out2_B_massfrac                                |              |              |              |
spl_out1_split_frac               ==            0.5|              |              |           0.5|
spl_out2_split_frac                                |              |              |           0.5|
11 variables
```

See also [`Stream`](@ref), [`Mixer`](@ref), [`Separator`](@ref), [`YieldReactor`](@ref), [`MultiYieldReactor`](@ref), [`StoicReactor`](@ref)
"""
struct Splitter <: AbstractBlock
    @Block_fields
    split_fracs::Dict{Symbol, VariableRef}

    function Splitter(m::Model,
                      name::Symbol,
                      fs::Flowsheet,
                      inlet::Stream,
                      outlets::Vector{Stream})
        @assert(all(strm.comps == inlet.comps for strm in outlets),
            """In Splitter $name, all the outlet streams must have the same component set
                as the inlet stream.""")

        inlets = [inlet]
        @Block_init

        # Make split frac variables.
        split_fracs = Dict{Symbol, VariableRef}()
        f = 1.0/length(outlets)
        for s in outlets
            v = make_var!(m, prefix * string(s.name) * "_split_frac", var_list)
            fix(v, f)
            split_fracs[s.name] = v
        end
        free(split_fracs[outlets[end].name])

        # Total mass balance, e.g., spl_out1_mass + spl_out2_mass == spl_in_mass
        sin = inlet
        eq_name = prefix * "total_mass_balance"
        make_eq!(m, eq_name,
            @constraint(m, 
                sum(strm_vars[s.name][:total_mass] for s in outlets) == strm_vars[sin.name][:total_mass],
            base_name=eq_name), eq_list)

        for sout in outlets
            # Split frac definition, e.g., spl_in_mass * spl_out1_split_frac == out1_mass
            eq_name = prefix * string(sout.name) * "_split_frac_def"
            make_eq!(m, eq_name,
                @constraint(m, 
                    strm_vars[sin.name][:total_mass] * split_fracs[sout.name] == strm_vars[sout.name][:total_mass],
                base_name=eq_name), eq_list)
        
            # Outlet composition copy, e.g., spl_in_h2_massfrac == spl_out1_h2_massfrac
            #                            or  spl_in_h2_mass * spl_out1_mass == spl_out1_h2_mass * spl_in_mass
            for c in sin.comps
                eq_name = prefix * string(sout.name) * "_" * string(c) * "_copy"
                make_eq!(m, eq_name,
                    @constraint(m, 
                        strm_vars[sin.name][:fx][c] * (is_flow(sout) ? strm_vars[sout.name][:total_mass] : 1.0) ==
                            strm_vars[sout.name][:fx][c] * (is_flow(sin) ? strm_vars[sin.name][:total_mass] : 1.0),
                    base_name=eq_name), eq_list)
            end
        end

        self = new(name, fs, inlets, outlets, strm_vars, var_list, eq_list, split_fracs)
        @Block_finish
        return self
    end
end

function JuMP.set_start_values(blk::Splitter; copy_inlets::Bool=true)
    copy_inlets && copy_streams(blk.inlets)
    inlet = blk.inlets[1]
    in_vars = blk.strm_vars[inlet.name]

    # Values of free split fracs. Can't assume that only the last split_frac is free.
    set_start_values_of_free_fracs!(blk.split_fracs)

    in_mass_val = get_value(in_vars[:total_mass])
    in_fx = in_vars[:fx]
    for s in blk.outlets
        out_vars = blk.strm_vars[s.name]
        
        # Outlet stream total mass flow rates.
        out_mass_val = in_mass_val * get_value(blk.split_fracs[s.name])
        set_start_value(out_vars[:total_mass], out_mass_val)

        # Outlet stream mass fracs or component mass flow rates.
        for c in inlet.comps
            in_mass_c_val = max(1.0e-8, get_value(in_fx[c]))
            out_fx_val = in_mass_c_val * (is_flow(s) ? out_mass_val : 1.0) / (is_flow(inlet) ? in_mass_val : 1.0)
            set_start_value(out_vars[:fx][c], out_fx_val)
        end
    end
end
