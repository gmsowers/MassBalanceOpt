"""
    @Block_fields

Create the fields typically used in a block. They are:
```julia
    name::Symbol
    fs::Flowsheet
    inlets::Vector{Stream}
    outlets::Vector{Stream}
    strm_vars::Dict{Symbol, Dict}
    var_list::Vector{VariableRef}
    eq_list::Vector{ConstraintRef}
```

See also [`AbstractBlock`](@ref), [`@Block_init`](@ref), [`@Block_finish`](@ref)
"""
macro Block_fields()
    esc(quote
        name::Symbol
        fs::Flowsheet
        inlets::Vector{Stream}
        outlets::Vector{Stream}
        strm_vars::Dict{Symbol, Dict}
        var_list::Vector{VariableRef}
        eq_list::Vector{ConstraintRef}
    end)
end

"""
    @Block_init

Create the stream variables used in a block and set their specs.

See also [`AbstractBlock`](@ref), [`@Block_fields`](@ref), [`@Block_finish`](@ref)
"""
macro Block_init()
    esc(quote
        var_list = Vector{VariableRef}()
        eq_list = Vector{ConstraintRef}()
        prefix = fs.prefix * string(name) * "_"

        strm_vars = make_stream_vars!(m, [inlets; outlets], prefix, var_list)
        set_stream_var_specs!(inlets, strm_vars)
    end)
end

"""
    @Block_finish

Assign a newly created block referenced by `self` to its container flowsheet, and
record the block in the `to` and `from` fields of the inlet and outlet streams.

See also [`AbstractBlock`](@ref), [`@Block_fields`](@ref), [`@Block_init`](@ref)
"""
macro Block_finish()
    esc(quote
        fs.blocks[name] = self
        for s in inlets s.to = self end
        for s in outlets s.from = self end
    end)
end

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

"""
    set_start_values(blk::AbstractBlock; copy_inlets::Bool=true)

Set the start values of the block variables in `blk`. If `copy_inlets=false`, don't copy the inlet stream variable start values
from the upstream block.
"""
function JuMP.set_start_values(blk::AbstractBlock; copy_inlets::Bool=true)
end

"""
    set_start_values(blks::Vector{<:AbstractBlock})

Set the start values of the block variables in all the blocks in `blks`.
"""
JuMP.set_start_values(blks::Vector{<:AbstractBlock}) = foreach(set_start_values, blks)


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

"""
    set_start_values_of_free_fracs!(fracs)

Set the start values of all the free variables in `fracs`.
"""
function set_start_values_of_free_fracs!(fracs)
    fracs_c = collect(values(fracs))
    fixed_fracs = filter(is_fixed, fracs_c)
    (length(fixed_fracs) == length(fracs_c)) && return  # all the fracs are fixed
    fixed_frac_sum = (isempty(fixed_fracs) ? 0.0 : sum(map(get_value, fixed_fracs)))
    free_fracs = filter(!is_fixed, fracs_c)
    free_frac_val = (1.0 - fixed_frac_sum)/length(free_fracs)
    foreach(free_fracs) do var set_start_value(var, free_frac_val) end
    return
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

"""
    make_mass_flow_vars!(m::Model,
                         prefix::AbstractString,
                         inlets::Vector{Stream},
                         outlets::Vector{Stream},
                         strm_vars::Dict{Symbol, Dict},
                         var_list::Vector{VariableRef},
                         eq_list::Vector{ConstraintRef})::Dict{Symbol, Dict}

Conditionally create component mass flow rate variables for the streams in
`inlets` and `outlets`. If the stream has a `FLOW` basis, component mass flow rate variables already exist.
If the stream has a `FRAC` basis, make new component mass flow variables and equations
to calculate them. Add equations to calculate the outlet stream total mass flow rates
from the component mass flow rates.
"""
function make_mass_flow_vars!(m::Model,
                              prefix::AbstractString,
                              inlets::Vector{Stream},
                              outlets::Vector{Stream},
                              strm_vars::Dict{Symbol, Dict},
                              var_list::Vector{VariableRef},
                              eq_list::Vector{ConstraintRef})::Dict{Symbol, Dict}
    mass_flows = Dict{Symbol, Dict}()
    for s in union(inlets, outlets)
        s_mass_flows = mass_flows[s.name] = Dict{Symbol, VariableRef}()
        s_is_frac = is_frac(s)
        for c in s.comps
            if s_is_frac
                prefix_ = prefix * string(s.name) * '_' * string(c)
                s_mass_flows[c] = make_var!(m, prefix_ * "_mass", var_list)

                # e.g., blk_in_ch4_mass == blk_in_ch4_massfrac * blk_in_mass
                eq_name = prefix_ * "_mass_flow_calc"
                make_eq!(m, eq_name,
                    @constraint(m, 
                        s_mass_flows[c] == strm_vars[s.name][:fx][c] * strm_vars[s.name][:total_mass],
                    base_name=eq_name), eq_list)
            else
                s_mass_flows[c] = strm_vars[s.name][:fx][c]
            end
        end
    end

    # e.g., blk_out1_h2_mass + blk_out1_ch4_mass == blk_out1_mass
    for s in outlets
        eq_name = prefix * string(s.name) * "_total_mass_calc"
        make_eq!(m, eq_name,
            @constraint(m, 
                sum(mass_flows[s.name][c] for c in s.comps) == strm_vars[s.name][:total_mass],
            base_name=eq_name), eq_list)
    end

    return mass_flows
end

last_value = last ∘ collect ∘ values

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