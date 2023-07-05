module MassBalanceOpt

using JuMP, Ipopt
using OrderedCollections
using Printf

export StreamBasis,
       FRAC,
       FLOW,
       AbstractBlock,
       @components,
       make_var!,
       print_vars,
       print_fixed,
       print_free,
       print_bounds,
       print_active,
       eval_obj,
       write_vars,
       free,
       set_lower,
       set_upper,
       delete_lower,
       delete_upper,
       flip,
       connect,
       Stream,
       is_frac,
       is_flow,
       Flowsheet,
       make_eq!,
       print_eqs,
       print_model,
       copy_stream,
       copy_streams,
       make_stream_vars!,
       set_stream_var_specs!,
       @Block_fields,
       @Block_init,
       @Block_finish,
       Mixer,
       Splitter,
       Separator,
       YieldReactor,
       MultiYieldReactor,
       make_mass_flow_vars!,
       @stoic,
       StoicReactor,
       @solve,
       @set,
       @values,
       @specs,
       @bounds,
       @stream,
       @streams,
       @block
    
"""
    @enum StreamBasis

`FRAC`: Generate mass fraction variables (default).
`FLOW`: Generate component mass flow rate variables.

See also [`Stream`](@ref)
"""
@enum StreamBasis FRAC FLOW

"""
    AbstractBlock

Abstract supertype for all blocks (Mixer, Splitter, etc.)

See also [`Mixer`](@ref), [`Splitter`](@ref), [`Separator`](@ref), [`YieldReactor`](@ref),
    [`MultiYieldReactor`](@ref), [`StoicReactor`](@ref), [`@Block_fields`](@ref),
    [`@Block_init`](@ref), [`@Block_finish`](@ref)
"""
abstract type AbstractBlock end

"""
    make_var!(m::Model, name::AbstractString, var_list::Vector{VariableRef})

Create a new JuMP variable named `name`, register it in model `m` as `m[:name]`, and add it to `var_list`.
"""
function make_var!(m::Model,
                   name::AbstractString,
                   var_list::Vector{VariableRef})
    vref = @variable(m, base_name=name)
    m[Symbol(name)] = vref
    push!(var_list, vref)
    return vref
end

function print_var(io::IO, v::VariableRef)
    fixed = is_fixed(v)
    _name = name(v)
    _value = if has_values(owner_model(v))
                @sprintf("%14.7g", value(v))
            else
                (fixed ? @sprintf("%14.7g", fix_value(v)) : "")
            end
    lower = (has_lower_bound(v) ? @sprintf("%14.7g", lower_bound(v)) : "")
    upper = (has_upper_bound(v) ? @sprintf("%14.7g", upper_bound(v)) : "")
    start = (has_start_value(v) ? @sprintf("%14.7g", start_value(v)) : "")
    spec = (fixed ? "  == " : "")
    @printf(io, "%-32s%5s%14s|%14s|%14s|%14s|\n", _name, spec, _value, lower, upper, start)
    return nothing
end
print_var(v::VariableRef) = print_var(stdout, v)

"""
    print_vars([io::IO], vars::Vector{VariableRef})

Print the variables in `vars` in a table.

See also [`print_fixed`](@ref), [`print_free`](@ref), [`print_bounds`](@ref), [`print_active`](@ref), [`print_eqs`](@ref),
    [`write_vars`](@ref)
"""
function print_vars(io::IO, vars::Vector{VariableRef})
    n = length(vars)
    n == 0 && return
    println(io, "              Name               Fix      Value          Lower          Upper          Start     ")
    println(io, "--------------------------------|---|--------------|--------------|--------------|--------------|")
    foreach(vars) do var print_var(io, var) end
    suffix = (n > 1 ? 's' : "")
    println(io, "$n variable$suffix")
    return
end
print_vars(vars::Vector{VariableRef}) = print_vars(stdout, vars)

"""
    print_vars([io::IO], var::VariableRef)

Print a single variable `var`.
"""
print_vars(io::IO, var::VariableRef) = print_vars(io, [var])
print_vars(var::VariableRef) = print_vars(stdout, var)

"""
    print_vars([io::IO], m::Model)

Print all the variables in model `m`.
"""
print_vars(io::IO, m::Model) = print_vars(io, all_variables(m))
print_vars(m::Model) = print_vars(stdout, m)

"""
    print_vars([io::IO], blk::AbstractBlock)

Print all the variables in block `blk`.
"""
print_vars(io::IO, blk::AbstractBlock) = print_vars(io, blk.var_list)
print_vars(blk::AbstractBlock) = print_vars(stdout, blk)

"""
    print_vars([io::IO], m::Model, glob::AbstractString)

Print the variables in model `m` that match `glob`, where `glob` is a string
containing one or more `*` wildcard characters.

### Examples
```julia
julia> print_vars(m, "*mass")

julia> print_vars(m, "*cgcfeed*")
```
"""
function print_vars(io::IO, m::Model, glob::AbstractString)
    if occursin("*", glob)
        print_vars(io, m, Regex(replace(glob, "*" => ".*")))
    else
        print_vars(io, m, Regex(glob))
    end
end
print_vars(m::Model, glob::AbstractString) = print_vars(stdout, m, glob)

"""
    print_vars([io::IO], m::Model, pattern::Regex)

Print the variables in model `m` that match the regular expression `pattern`.

### Examples
```julia
julia> print_vars(m, r"^.*mass\$")
```
"""
function print_vars(io::IO, m::Model, pattern::Regex)
    vars = Vector{VariableRef}()
    for var in all_variables(m)
        occursin(pattern, name(var)) && push!(vars, var)
    end
    print_vars(io, vars)
end
print_vars(m::Model, pattern::Regex) = print_vars(stdout, m, pattern)


"""
    print_fixed([io::IO], m::Model)

Print all the fixed variables in model `m`.

See also [`print_free`](@ref), [`print_bounds`](@ref), [`print_active`](@ref)
"""
print_fixed(io::IO, m::Model) = print_vars(io, filter(is_fixed, all_variables(m)))
print_fixed(m::Model) = print_fixed(stdout, m)

"""
    print_fixed([io::IO], blk::AbstractBlock)

Print all the fixed variables in block `blk`.
"""
print_fixed(io::IO, blk::AbstractBlock) = print_vars(io, filter(is_fixed, blk.var_list))
print_fixed(blk::AbstractBlock) = print_fixed(stdout, blk)

"""
    print_free([io::IO], m::Model)

Print all the free variables in model `m`.

See also [`print_fixed`](@ref), [`print_bounds`](@ref), [`print_active`](@ref)
"""
print_free(io::IO, m::Model) = print_vars(io, filter(!is_fixed, all_variables(m)))
print_free(m::Model) = print_free(stdout, m)

"""
    print_free([op::IO], blk::AbstractBlock)

Print all the free variables in block `blk`.
"""
print_free(io::IO, blk::AbstractBlock) = print_vars(io, filter(!is_fixed, blk.var_list))
print_free(blk::AbstractBlock) = print_free(stdout, blk)

bounded_vars(vars::Vector{VariableRef}) = filter(vars) do var
                                                           (has_lower_bound(var) || has_upper_bound(var))
                                                       end

"""
    print_bounds([io::IO], m::Model)

Print the variables in `m` with lower or upper bounds.

See also [`print_fixed`](@ref), [`print_free`](@ref), [`print_active`](@ref)
"""
print_bounds(io::IO, m::Model) = print_vars(io, bounded_vars(all_variables(m)))
print_bounds(m::Model) = print_bounds(stdout, m)

"""
    print_bounds([io::IO], blk::AbstractBlock)

Print the variables in `blk` with lower or upper bounds.
"""
print_bounds(io::IO, blk::AbstractBlock) = print_vars(io, bounded_vars(blk.var_list))
print_bounds(blk::AbstractBlock) = print_bounds(stdout, blk)

const tol = 2.0 * sqrt(eps())

"""
    print_active([io::IO], m::Model)

Print the variables in `m` whose values equal their lower or upper bounds.

See also [`print_fixed`](@ref), [`print_free`](@ref), [`print_bounds`](@ref)
"""
function print_active(io::IO, m::Model)
    has_values(m) || return
    print_vars(io, filter(bounded_vars(all_variables(m))) do var
                            (has_lower_bound(var) && (isapprox(value(var), lower_bound(var); atol=tol, rtol=tol))) ||
                            (has_upper_bound(var) && (isapprox(value(var), upper_bound(var); atol=tol, rtol=tol)))
                        end)
end
print_active(m::Model) = print_active(stdout, m)

"""
    eval_obj(m::Model)

Evalulate the objective function in `m`, using the solution values if a solution exists, otherwise
using the start values.
"""
eval_obj(m::Model) = value(get_value, objective_function(m))

"""
    write_vars([io::IO], m::Model)

Write the values of all the variables in `m` to `io` in the format:
```julia
@values begin
    var_1 = 1.0
    var_2 = 2.0
     ...
end
```
"""
function write_vars(io::IO, m::Model)
    max_len = maximum(length.(name.(all_variables(m))))
    println(io, "@values begin")
    foreach(all_variables(m)) do var
        Printf.format(io, Printf.Format("    %-$(max_len)s = %s\n"), name(var), Base.sprint(Base.show, get_value(var)))
    end
    print(io, "end")
end
write_vars(m::Model) = write_vars(stdout, m)

"""
    get_value(var::VariableRef)

Return `var`'s:
    - solution value, if a solution exists
    - fixed value, if `var` is fixed
    - start value, if `var` has a start value
    - 0.0
in that order.
"""
function get_value(var::VariableRef)
    return (has_values(owner_model(var)) ? value(var) :
           (is_fixed(var)                ? fix_value(var) :
           (has_start_value(var)         ? start_value(var) : 0.0)))
end

"""
    fix(var::VariableRef, [val::Real])

If `val` is passed, fix the variable `var` equal to `val` and set its start value to `val`.
If no `val` is passed, fix the variable equal to its value in the most recent
solution. If no solution is available, fix the variable equal to its
start value. If no start value is available, fix the variable equal to 0.

See also [`free`](@ref), [`@specs`](@ref)
"""
JuMP.fix(var::VariableRef) = fix(var, get_value(var))
function JuMP.fix(var::VariableRef, val::Real)
    fix(var, val, force=true)
    set_start_value(var, val)
end

"""
    free(var::VariableRef)

Delete the fixing constraint for variable `var`.

See also [`fix`](@ref), [`@specs`](@ref)
"""
free(var::VariableRef) = is_fixed(var) && unfix(var)

"""
    set_value(var::VariableRef, val::Real)

If `var` is fixed, reset its fixed value to `val`. If `var` is free, set
its start value to `val`.

See also [`fix`](@ref), [`@set`](@ref), [`@values`](@ref)
"""
function JuMP.set_value(var::VariableRef, val::Real)
    is_fixed(var) && fix(var, val)
    set_start_value(var, val)
    return
end

"""
    flip(var1::VariableRef, [var2::VariableRef])

If `var2` is passed, invert the specs of `var1` and `var2`, i.e., if `var1` is fixed and
`var2` is free, make `var1` free and `var2` fixed; if `var1` is free and
`var2` is fixed, make `var1` fixed and `var2` free. If both vars have the same spec,
do nothing.

If `var2` is not passed, invert the spec of `var1`.

!!! warning

    The newly fixed variable will have its value set to:
    
        1. its value in the most recent solution
        2. its start value
        3. 0
        
    in that order. Use `set_value` or `@set` to change the value of the newly fixed variable.

See also [`@set`](@ref), [`@specs`](@ref)
"""
function flip(var1::VariableRef, var2::VariableRef)
    isfixed = (is_fixed(var1), is_fixed(var2))
    xor(isfixed...) || return
    fix(isfixed[1] ? var2 : var1)
    free(isfixed[1] ? var1 : var2)
    return
end
flip(var::VariableRef) = is_fixed(var) ? free(var) : fix(var)

"""
    connect(var_lhs::VariableRef, var_rhs::VariableRef)

Given two variables `var_lhs` and `var_rhs`, at least one of which is fixed, add a new constraint:

    var_lhs == var_rhs

If both variables are fixed, free `var_rhs`. Otherwise free the single fixed variable. If both
variables are free, don't add a new constraint and return `nothing`.

Return the newly created ConstraintRef.
See also [`@set`](@ref)
"""
function connect(var_lhs::VariableRef, var_rhs::VariableRef)
    is_fixed(var_lhs) || is_fixed(var_rhs) || return nothing
    m = owner_model(var_lhs)
    (m == owner_model(var_rhs)) || return nothing
    unfix(is_fixed(var_rhs) ? var_rhs : var_lhs)
    eq_sym = Symbol(:cn, gensym())
    eq_name = string(eq_sym)
    m[eq_sym] = @constraint(m, var_lhs == var_rhs, base_name=eq_name)
    return m[eq_sym]
end

"""
Representation of a material stream in a mass balance model.

    Stream(name::Symbol, fs, comps; basis=FRAC)

Create a new `Stream` with `name` and components `comps` in flowsheet `fs`. The keyword
argument `basis` can be `FRAC` or `FLOW`, with `FRAC` being the default.

Creating a `Stream` does not create any variables or equations in the model. The stream becomes available
for use as an argument to the constructors of blocks, which are subtyped of `AbstractBlock`. Blocks are
responsible for creating stream variables and equations. When a block creates stream variables 
it uses the basis of each stream to determine what variables to create. If
`basis=FRAC` a mass fraction variable is created for every component in the stream; if `basis=FLOW` a
component mass flow rate variable is created for every component.

!!! warning

    None of the built-in blocks create equations to sum the mass fractions or the component mass flow rates.

### Examples

Create two inlet streams and one outlet stream, then create a `Mixer` block that mixes the two inlet streams.

```julia
julia> m = Model(); fs = Flowsheet();

julia> in1 = Stream(:in1, fs, @components(h2, ch4))
Stream(in1, fs=index, basis=FRAC, comps=[h2, ch4])

julia> in2 = Stream(:in2, fs, @components(h2, ch4, c2h4), basis=FLOW)
Stream(in2, fs=index, basis=FLOW, comps=[h2, ch4, c2h4])

julia> out = Stream(:out, fs, @components(h2, ch4, c2h4))
Stream(out, fs=index, basis=FRAC, comps=[h2, ch4, c2h4])

julia> mix1 = Mixer(m, :mix1, fs, [in1, in2], out)
Mixer(in=[in1, in2], out=[out])

julia> print_vars(mix1)
              Name              |    Value   |    Lower   |    Upper   |    Start   |Spec
--------------------------------|------------|------------|------------|------------|----
mix1_in1_mass                   |           1|            |            |            |F
mix1_in1_h2_massfrac            |         0.5|            |            |            |F
mix1_in1_ch4_massfrac           |         0.5|            |            |            |F
mix1_in2_mass                   |           1|            |            |            |F
mix1_in2_h2_mass                |     0.33333|            |            |            |F
mix1_in2_ch4_mass               |     0.33333|            |            |            |F
mix1_in2_c2h4_mass              |     0.33333|            |            |            |F
mix1_out_mass                   |            |            |            |            |
mix1_out_h2_massfrac            |            |            |            |            |
mix1_out_ch4_massfrac           |            |            |            |            |
mix1_out_c2h4_massfrac          |            |            |            |            |
11 variables
```

See also [`@stream`](@ref), [`@streams`](@ref), [`Flowsheet`](@ref), [`@components`](@ref), [`StreamBasis`](@ref)
"""
mutable struct Stream
    name::Symbol
    fs
    basis::StreamBasis
    comps::OrderedSet{Symbol}
    to::Union{AbstractBlock, Nothing}
    from::Union{AbstractBlock, Nothing}

    function Stream(name::Symbol,
                    fs,
                    comps;
                    basis=FRAC)
        @assert(!isempty(comps), "Stream $name must contain at least one component.")
        self = new(name, fs, basis, comps, nothing, nothing)
        fs.streams[name] = self
        return self
    end
end

"""
    is_frac(strm::Stream)

Return `true` if the `strm` basis=`FRAC`, otherwise return false.

See also [`is_flow`](@ref), [`StreamBasis`](@ref)
"""
is_frac(strm::Stream) = strm.basis == FRAC

"""
    is_flow(strm::Stream)

Return `true` if the `strm` basis=`FLOW`, otherwise return false.

See also [`is_frac`](@ref), [`StreamBasis`](@ref)
"""
is_flow(strm::Stream) = strm.basis == FLOW

function Base.show(io::IO, strm::Stream)
    s = "Stream($(strm.name), fs=$(strm.fs.name), basis=$(strm.basis)"
    isnothing(strm.from) || (s *= ", from=$(strm.from.name)")
    isnothing(strm.to) || (s *= ", to=$(strm.to.name)")
    length(strm.comps) > 0 && (s *= ", components=[" * join(string.(strm.comps), ", ") * ']')
    s *= ')'
    print(io, s)
end

"""
    connect(strm::Stream)

Connect all the variables in the source and destination blocks of `strm`. If the stream
has no source, or no destination, do nothing.
"""
function connect(strm::Stream)
    (isnothing(strm.to) || isnothing(strm.from)) && return
    to_blk = strm.to
    eq = connect(strm.from.strm_vars[strm.name][:total_mass], strm.to.strm_vars[strm.name][:total_mass])
    isnothing(eq) && return
    push!(to_blk.eq_list, eq)
    for c in strm.comps
        eq = connect(strm.from.strm_vars[strm.name][:fx][c], strm.to.strm_vars[strm.name][:fx][c])
        isnothing(eq) || push!(to_blk.eq_list, eq)
    end
    return
end

"""
    connect(m::Model, strms::Vector{Stream})

Connect all the streams in the array `strms`.
"""
connect(strms::Vector{Stream}) = foreach(connect, strms)

"""
    copy_stream(strm::Stream)

Copy the values of the `strm` variables in `strm`'s `from` block into the values of the `strm` variables in
`strm`'s `to` block. If either the `to` or `from` field of `strm` is empty, do nothing.

See also [`copy_streams`](@ref)
"""
function copy_stream(strm::Stream)
    (isnothing(strm.from) || isnothing(strm.to)) && return
    from_blk_vars = strm.from.strm_vars[strm.name]
    to_blk_vars = strm.to.strm_vars[strm.name]
    set_value(to_blk_vars[:total_mass], get_value(from_blk_vars[:total_mass]))
    for (c, var) in from_blk_vars[:fx]
        set_value(to_blk_vars[:fx][c], get_value(var))
    end
    return
end

"""
    copy_streams(strms::Vector{Stream})

Invoke `copy_stream` for each element of `strms`.

See also [`copy_stream`](@ref)
"""
copy_streams(strms::Vector{Stream}) = foreach(copy_stream, strms)

"""
    Flowsheet(name::Symbol=:index, parent=nothing)

Create a flowsheet. If `name` is omitted it defaults to `:index`. If specified, `parent` must be another `Flowsheet`;
The new flowsheet will be a child of the parent flowsheet. If `parent` is omitted, the flowsheet has no parent.

A `Flowsheet` serves as a container for streams, blocks (subtypes of `AbstractBlock`),
and other flowsheets. The names of a flowsheet's ancestors are embedded in the variable and equation names of the
blocks and streams contained in the flowsheet. For example, if a flowsheet named `:unit1` is a child of the index
flowsheet (a flowsheet created with no parent), the variable and equation names in that flowsheet will start with
`unit1_`. If the `:unit1` flowsheet has a child named `reactors` the names will start with `unit1_reactors_`. 
Variables and equations in a flowsheet with no parent have no prefix.

### Examples

Create a `:index` flowsheet with no parent:
```julia
julia> fs = Flowsheet()
Flowsheet(
   name=index   )
```
Add a `Stream` to the `:index` flowsheet:
```julia
julia> in1 = Stream(:in1, fs, @components(h2,ch4))
Stream(in1, fs=index, basis=FRAC, components=[h2, ch4])

julia> fs
Flowsheet(
   name=index,
   streams=[in1]   )
```

Create a flowsheet named `:unit1` whose parent is `:index`:
```julia
julia> fs_1 = Flowsheet(:unit1, fs)
Flowsheet(
   name=unit1,
   parent=index   )

julia> fs
Flowsheet(
   name=index,
   children=[unit1],
   streams=[in1]   )
```

See also [`Stream`](@ref), [`AbstractBlock`](@ref)
"""
struct Flowsheet
    name::Symbol
    path::AbstractString
    prefix::AbstractString
    parent::Union{Flowsheet, Nothing}
    children::Vector{Flowsheet}
    blocks::OrderedDict{Symbol, AbstractBlock}
    streams::OrderedDict{Symbol, Stream}

    function Flowsheet(name::Symbol=:index,
                       parent=nothing)
        name_str = string(name)
        if isnothing(parent)
            path = name_str
            prefix = ""
        else
            path = parent.path * "_" * name_str
            prefix = parent.prefix * name_str * "_"
        end
        self = new(name, path, prefix, parent, Vector{Flowsheet}(), OrderedDict{Symbol, AbstractBlock}(), OrderedDict{Symbol, Stream}())
        !isnothing(parent) && push!(parent.children, self)
        return self
    end
end

function Base.show(io::IO, fs::Flowsheet)
    s = "Flowsheet(\n   name=$(fs.name)"
    isnothing(fs.parent) || (s *= ",\n   parent=$(fs.parent.name)")
    isempty(fs.children) || (s *= ",\n   children=[" * join(string.(getfield.(fs.children, :name)), ", ") * ']')
    isempty(fs.blocks)   || (s *= ",\n   blocks=[" * join(string.(getfield.(values(fs.blocks), :name)), ", ") * ']')
    isempty(fs.streams)  || (s *= ",\n   streams=[" * join(string.(getfield.(values(fs.streams), :name)), ", ") * ']')
    s *= "   )"
    print(io, s)
end

"""
    connect(fs::Flowsheet)

Connect all the streams in the flowsheet `fs`.
"""
function connect(fs::Flowsheet)
    connect(collect(values(fs.streams)))
    for fsht in fs.children
        connect(fsht)
    end
    return
end

"""
    make_stream_vars!(m::Model, strm::Stream, prefix::AbstractString, var_list::Vector{VariableRef})
        -> Dict{Symbol, Any}

Create stream variables for `strm` using `prefix` at the beginning of the variable names. Add the
new variables to `var_list`. Return the stream variables in a dictionary.

See also [`Stream`](@ref), [`set_stream_var_specs!`](@ref), [`make_var!`](@ref)
"""
function make_stream_vars!(m::Model,
                           strm::Stream,
                           prefix::AbstractString,
                           var_list::Vector{VariableRef})
    vars = Dict{Symbol, Any}()
    prefix_ = prefix * string(strm.name) * "_"
    vars[:total_mass] = make_var!(m, prefix_ * "mass", var_list)
    v = vars[:fx] = Dict{Symbol, VariableRef}()
    suffix = is_frac(strm) ? "_massfrac" : "_mass"
    for c in strm.comps
        v[c] = make_var!(m, prefix_ * string(c) * suffix, var_list)
    end
    return vars
end

"""
    make_stream_vars!(m::Model, strms::Vector{Stream}, prefix::AbstractString, var_list::Vector{VariableRef})
        -> Dict{Symbol, Dict}

Create stream variables for all the streams in the array `strms` and return them in a dictionary.
"""
function make_stream_vars!(m::Model,
                           strms::Vector{Stream},
                           prefix::AbstractString,
                           var_list::Vector{VariableRef})
    vars = Dict{Symbol, Dict}()
    for strm in strms
        vars[strm.name] = make_stream_vars!(m, strm, prefix, var_list)
    end
    return vars
end

"""
    set_stream_var_specs!(strms::Vector{Stream}, vars::Dict{Symbol, Dict})

Fix the variables stored in `vars`, for each stream in `strms`, equal to default values.

See also [`Stream`](@ref), [`make_stream_vars!`](@ref), [`fix`](@ref)
"""
function set_stream_var_specs!(strms::Vector{Stream}, vars::Dict{Symbol, Dict})
    for s in strms
        fix(vars[s.name][:total_mass], 1.0)
        f = 1.0/length(s.comps)
        for c in s.comps
            fix(vars[s.name][:fx][c], f)
        end
    end
    return nothing
end

"""
    make_eq!(m::Model, name::AbstractString, con::ConstraintRef, eq_list::Vector{ConstraintRef})

Register the equation `con` in model `m` using `name` as the base name. Add the equation to `eq_list`.

See also [`make_var!`](@ref)
"""
function make_eq!(m::Model,
                  name::AbstractString,
                  con::ConstraintRef,
                  eq_list::Vector{ConstraintRef})
    m[Symbol(name)] = con
    push!(eq_list, con)
    return nothing
end

"""
    print_eqs([io::IO], m::Model)

Print all the equations in model `m`.

See also [`print_vars`](@ref)
"""
function print_eqs(io::IO, m::Model)
    foreach(all_constraints(m, include_variable_in_set_constraints = true)) do eq println(io, eq) end
    println(io, "$(num_constraints(m, count_variable_in_set_constraints = true)) equations")
end
print_eqs(m::Model) = print_eqs(stdout, m)

"""
    print_eqs([io::IO], blk::AbstractBlock)

Print all the equations in block `blk`.
"""
function print_eqs(io::IO, blk::AbstractBlock)
    foreach(blk.eq_list) do eq println(io, eq) end
    println("$(length(blk.eq_list)) equations")
end
print_eqs(blk::AbstractBlock) = print_eqs(stdout, blk)

"""
    print_eqs([io::IO], m::Model, eq::Symbol)

Print the equation registered in model `m` as `m[:eq]`
"""
print_eqs(io::IO, m::Model, eq::Symbol) = println(io, constraint_by_name(m, string(eq)))
print_eqs(m::Model, eq::Symbol) = print_eqs(stdout, m, eq)

"""
    print_model([io::IO], model_or_block::Union{Model, AbstractBlock})

Print all the variables and equations in a model or block specified by `model_or_block`.

See also [`print_vars`](@ref), [`print_eqs`](@ref)
"""
function print_model(io::IO, model_or_block::Union{Model, AbstractBlock})
    print_vars(io, model_or_block)
    println(io)
    print_eqs(io, model_or_block)
end
print_model(model_or_block::Union{Model, AbstractBlock}) = print_model(stdout, model_or_block)

function Base.show(io::IO, blk::AbstractBlock)
    s = "$(nameof(typeof(blk)))(in=[$(blk.inlets[1].name)"
    for strm in blk.inlets[2:end]
        s *= ", $(strm.name)"
    end
    s *= "], out=[$(blk.outlets[1].name)"
    for strm in blk.outlets[2:end]
        s *= ", $(strm.name)"
    end
    s *= "])"
    print(io, s)
end

include("blocks.jl")
include("macros.jl")

end # module MassBalanceOpt
