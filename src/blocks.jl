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

builtin_block_implementations = (
                                 "Mixer.jl",
                                 "Splitter.jl",
                                 "Separator.jl",
                                 "YieldReactor.jl",
                                 "MultiYieldReactor.jl",
                                 "StoicReactor.jl",
                                 )
foreach(include, builtin_block_implementations)
