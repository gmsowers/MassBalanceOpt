"""
    Separator <: AbstractBlock

Separate the components in the inlet stream into two or more outlet streams.

Create a `Separator` block:
```julia
    Separator(m::Model, name::Symbol, fs::Flowsheet, inlet::Stream, outlets::Vector{Stream})
```
or if the current scope contains `m` and `fs` bound to a `Model` and `Flowsheet`:
```julia
    @block(name, Separator, inlet, outlets)
```

### Examples
```julia
julia> m = Model(); fs = Flowsheet(); comps = @components A B;

julia> (in1, out1, out2) = @streams begin
    in1, comps
    out1, comps
    out2, comps
end;

julia> sep = @block(sep, Separator, in1, [out1, out2]);

julia> @set sep_A_out1_split = 0.3

julia> @set sep_B_out1_split = 0.6

julia> set_start_values(sep)

julia> print_vars(sep)
              Name               Fix      Value          Lower          Upper          Start     
--------------------------------|---|--------------|--------------|--------------|--------------|
sep_in1_mass                      ==              1|              |              |             1|
sep_in1_A_massfrac                ==            0.5|              |              |           0.5|
sep_in1_B_massfrac                ==            0.5|              |              |           0.5|
sep_out1_mass                                      |              |              |          0.45|
sep_out1_A_massfrac                                |              |              |     0.3333333|
sep_out1_B_massfrac                                |              |              |     0.6666667|
sep_out2_mass                                      |              |              |          0.55|
sep_out2_A_massfrac                                |              |              |     0.6363636|
sep_out2_B_massfrac                                |              |              |     0.3636364|
sep_in1_A_mass                                     |              |              |           0.5|
sep_in1_B_mass                                     |              |              |           0.5|
sep_out1_A_mass                                    |              |              |          0.15|
sep_out1_B_mass                                    |              |              |           0.3|
sep_out2_A_mass                                    |              |              |          0.35|
sep_out2_B_mass                                    |              |              |           0.2|
sep_A_out1_split                  ==            0.3|              |              |           0.3|
sep_A_out2_split                                   |              |              |           0.7|
sep_B_out1_split                  ==            0.6|              |              |           0.6|
sep_B_out2_split                                   |              |              |           0.4|
19 variables
```

See also [`Stream`](@ref), [`Mixer`](@ref), [`Splitter`](@ref), [`YieldReactor`](@ref), [`MultiYieldReactor`](@ref), [`StoicReactor`](@ref)
"""
struct Separator <: AbstractBlock
    @Block_fields
    mass_flows::Dict{Symbol, Dict}
    splits::OrderedDict{Symbol, Dict}

    function Separator(m::Model,
                       name::Symbol,
                       fs::Flowsheet,
                       inlet::Stream,
                       outlets::Vector{Stream})
        @assert(union([strm.comps for strm in outlets]...) == inlet.comps,
        """In Separator $name, the inlet stream component set is not the same as
            the union of the outlet stream component sets.""")
           
        inlets = [inlet]
        @Block_init

        mass_flows = make_mass_flow_vars!(m, prefix, inlets, outlets, strm_vars, var_list, eq_list)

        # If a component appears in an outlet stream, make a split variable for that component
        #   and outlet stream along with an equation that defines it.
        splits = OrderedDict{Symbol, Dict}()
        sin = inlet
        for c in sin.comps
            splits[c] = Dict{Symbol, VariableRef}()
            for s in outlets
                if c in s.comps
                    prefix_ = prefix * string(c) * '_' * string(s.name)
                    splits[c][s.name] = make_var!(m, prefix_ * "_split", var_list)

                    # e.g., sep_in_h2_mass * sep_h2_out1_split == sep.out1_h2_mass
                    eq_name = prefix_ * "_split_calc"
                    make_eq!(m, eq_name,
                        @constraint(m, 
                            mass_flows[sin.name][c] * splits[c][s.name] == mass_flows[s.name][c],
                        base_name=eq_name), eq_list)
                end
            end

            # Add an equation to sum the splits to one for each component.
            #   e.g., sep_h2_out1_split + sep_h2_out2_split == 1.0
            if !isempty(splits[c])
                eq_name = prefix * string(c) * "_sum_splits"
                make_eq!(m, eq_name,
                    @constraint(m, 
                        sum(splits[c][s] for s in keys(splits[c])) == 1.0,
                    base_name=eq_name), eq_list)
            end
        end

        # Fix all but the last split variable for each component in each outlet stream.
        for (c, c_splits) in splits
            n_splits = length(c_splits)
            if n_splits > 0
                split = 1.0/n_splits
                for split_var in values(c_splits)
                    fix(split_var, split)
                end
                free(last_value(splits[c]))
            end
        end

        self = new(name, fs, inlets, outlets, strm_vars, var_list, eq_list, mass_flows, splits)
        @Block_finish
        return self
    end
end

function JuMP.set_start_values(blk::Separator; copy_inlets::Bool=true)
    copy_inlets && copy_streams(blk.inlets)
    inlet = blk.inlets[1]
    in_vars = blk.strm_vars[inlet.name]
    in_mass_flows = blk.mass_flows[inlet.name]

    # Values of free split fracs. Can't assume that only the last split_frac is free.
    foreach(set_start_values_of_free_fracs!, values(blk.splits))
    
    # Values of the block inlet component mass flow rates. If the inlet stream is FLOW based, these will already exist.
    if is_frac(inlet)
        for c in inlet.comps
            set_start_value(in_mass_flows[c], get_value(in_vars[:fx][c]) * get_value(in_vars[:total_mass]))
        end
    end

    # Values of the block outlet component mass flow rates.
    for c in inlet.comps
        for s in blk.outlets
            if c in s.comps
                set_start_value(blk.mass_flows[s.name][c], get_value(in_mass_flows[c]) * get_value(blk.splits[c][s.name]))
            end
        end
    end

    # Values of the outlet stream total mass flow rates, and mass fractions or component mass flows.
    for s in blk.outlets
        s_mass_flows = blk.mass_flows[s.name]
        s_strm_vars = blk.strm_vars[s.name]
        set_start_value(s_strm_vars[:total_mass], sum(get_value(s_mass_flows[c]) for c in s.comps))
        s_is_frac = is_frac(s)
        for c in s.comps
            if s_is_frac
                set_start_value(s_strm_vars[:fx][c], get_value(s_mass_flows[c]) / max(1.0e-8, get_value(s_strm_vars[:total_mass])))
            else
                set_start_value(s_strm_vars[:fx][c], get_value(s_mass_flows[c]))
            end
        end
    end
end

