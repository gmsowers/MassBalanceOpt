"""
    set_lower(var::VariableRef, lower::Real)

If `var` is not fixed, set its lower bound to `lower`.

See also [`set_upper`](@ref), [`delete_lower`](@ref)
"""
set_lower(var::VariableRef, lower::Real) = is_fixed(var) || set_lower_bound(var, lower)

"""
    set_upper(var::VariableRef, lower::Real)

If `var` is not fixed, set its upper bound to `upper`.

See also [`set_lower`](@ref), [`delete_upper`](@ref)
"""
set_upper(var::VariableRef, upper::Real) = is_fixed(var) || set_upper_bound(var, upper)

"""
    delete_lower(var::VariableRef)

If `var` is not fixed and has a lower bound, delete its lower bound.

See also [`set_lower`](@ref), [`delete_upper`](@ref)
"""
delete_lower(var::VariableRef) = !is_fixed(var) && has_lower_bound(var) && delete_lower_bound(var)

"""
    delete_upper(var::VariableRef)

If `var` is not fixed and has an upper bound, delete its upper bound.

See also [`set_upper`](@ref), [`delete_lower`](@ref)
"""
delete_upper(var::VariableRef) = !is_fixed(var) && has_upper_bound(var) && delete_upper_bound(var)

const m_fs_note = """assumes that the model is stored in a variable named `m` and the
active flowsheet is stored in a variable named `fs` in the current scope"""
const m_note = "assumes that the model is stored in a variable named `m` in the current scope"

function transform_variable_expr(ex)
    err = "Unable to transform expression `$(string(ex))`\n"
    if ex.head == :(=)
        err *= """Expressions using `=` must look like:
                              x = 1.0   -or-
                            1.0 = x     -or-
                              x = y
                where `1.0` is any number and `x` and `y` are symbols"""
        if ex.args[1] isa Symbol
            if ex.args[2] isa Real
                return esc(:(set_value(m[$(QuoteNode(ex.args[1]))], $(ex.args[2]))))
            elseif ex.args[2] isa Symbol
                return esc(:(connect(m[$(QuoteNode(ex.args[1]))], m[$(QuoteNode(ex.args[2]))])))
            else
                error(err)
            end
        elseif ex.args[1] isa Real
            if ex.args[2] isa Symbol
                return esc(:(set_value(m[$(QuoteNode(ex.args[2]))], $(ex.args[1]))))
            else
                error(err)
            end
        else
            error(err)
        end
    elseif ex.head == :call && ex.args[1] isa Symbol
        if ex.args[1] == :(~)
            err *= """Expressions using `~` must look like:
                                 ~x      -or-
                                 x ~ y
                    where `x` and `y` are symbols"""
            if ex.args[2] isa Symbol
                if length(ex.args) == 3 && ex.args[3] isa Symbol
                    return esc(:(flip(m[$(QuoteNode(ex.args[2]))], m[$(QuoteNode(ex.args[3]))])))
                elseif length(ex.args) == 2
                    return esc(:(flip(m[$(QuoteNode(ex.args[2]))])))
                else
                    error(err)
                end
            else
                error(err)
            end
        elseif ex.args[1] == :(+) || ex.args[1] == :(-)
            err *= """Expressions using `+` or `-` must look like:
                                    +x   -or-
                                    -x
                    where `x` is a symbol"""
            if ex.args[2] isa Symbol
                if length(ex.args) == 2
                    op = (ex.args[1] == :(+) ? fix : free)
                    return esc(:($op(m[$(QuoteNode(ex.args[2]))])))
                else
                    error(err)
                end
            else
                error(err)
            end
        elseif ex.args[1] == :(>) || ex.args[1] == :(>=)
            err *= """Binary expressions using `>` must look like:
                              x > 1.0   -or-
                            1.0 > x     -or-
                            Inf > x     -or-
                              x > -Inf
                    where `1.0` is any number and `x` is a symbol"""
            if ex.args[2] isa Symbol
                if ex.args[2] == :Inf && ex.args[3] isa Symbol
                    return esc(:(delete_upper(m[$(QuoteNode(ex.args[3]))])))
                elseif ex.args[3] isa Expr && ex.args[3].head == :call
                    subex = ex.args[3]
                    if subex.args[1] == :(-) && length(subex.args) == 2 && subex.args[2] == :Inf
                        return esc(:(delete_lower(m[$(QuoteNode(ex.args[2]))])))
                    else
                        error(err)
                    end
                elseif ex.args[3] isa Real
                    return esc(:(set_lower(m[$(QuoteNode(ex.args[2]))], $(ex.args[3]))))
                else
                    error(err)
                end
            elseif ex.args[2] isa Real && ex.args[3] isa Symbol
                    return esc(:(set_upper(m[$(QuoteNode(ex.args[3]))], $(ex.args[2]))))
            else
                error(err)
            end
        elseif ex.args[1] == :(<) || ex.args[1] == :(<=)
            err *= """Binary expressions using `<` must look like:
                              x < 1.0   -or-
                            1.0 < x     -or-
                           -Inf < x     -or-
                              x < Inf
                    where `1.0` is any number and `x` is a symbol"""
            if ex.args[2] isa Symbol
                if ex.args[3] == :Inf
                    return esc(:(delete_upper(m[$(QuoteNode(ex.args[2]))])))
                elseif ex.args[3] isa Real
                    return esc(:(set_upper(m[$(QuoteNode(ex.args[2]))], $(ex.args[3]))))
                else
                    error(err)
                end
            elseif ex.args[2] isa Expr && ex.args[2].head == :call
                subex = ex.args[2]
                if subex.args[1] == :(-) && length(subex.args) == 2 && subex.args[2] == :Inf
                    return esc(:(delete_lower(m[$(QuoteNode(ex.args[3]))])))
                else
                    error(err)
                end
            elseif ex.args[2] isa Real && ex.args[3] isa Symbol
                    return esc(:(set_lower(m[$(QuoteNode(ex.args[3]))], $(ex.args[2]))))
            else
                error(err)
            end
        else
            error(err)
        end
    elseif ex.head == :comparison
        err *= """Intervals must look like:
                    1.0 < x < 2.0   -or-
                    2.0 > x > 1.0   -or-
                   -Inf < x < 1.0   -or-
                    1.0 < x < Inf   -or-
                    Inf > x > 1.0   -or-
                    1.0 > x > -Inf  -or-
                   -Inf < x < Inf   -or-
                    Inf > x > -Inf
                where `1.0` and `2.0` are any numbers such that 1.0 < 2.0, and `x` is a symbol"""
        if length(ex.args) == 5 && ex.args[3] isa Symbol
            if (ex.args[2] == :(<) || ex.args[2] == :(<=)) && (ex.args[4] == :(<) || ex.args[4] == :(<=))
                if ex.args[1] isa Expr
                    if ex.args[1].head == :call && ex.args[1].args[1] == :(-) && ex.args[1].args[2] == :Inf
                        _low = -Inf
                    else
                        error(err)
                    end
                elseif ex.args[1] isa Real
                    _low = ex.args[1]
                else
                    error(err)
                end
                if ex.args[5] isa Symbol || ex.args[5] isa Real
                    _up = ex.args[5]
                else
                    error(err)
                end
            elseif (ex.args[2] == :(>) || ex.args[2] == :(>=)) && (ex.args[4] == :(>) || ex.args[4] == :(>=))
                if ex.args[5] isa Expr
                    if ex.args[5].head == :call && ex.args[5].args[1] == :(-) && ex.args[5].args[2] == :Inf
                        _low = -Inf
                    else
                        error(err)
                    end
                elseif ex.args[5] isa Real
                    _low = ex.args[5]
                else
                    error(err)
                end
                if ex.args[1] isa Symbol || ex.args[1] isa Real
                    _up = ex.args[1]
                else
                    error(err)
                end
            else
                error(err)
            end
            res = Expr(:block)
            if _low === -Inf
                push!(res.args, quote delete_lower(m[$(QuoteNode(ex.args[3]))]) end)
            else
                push!(res.args, quote set_lower(m[$(QuoteNode(ex.args[3]))], $(_low)) end)
            end
            if _up === :Inf
                push!(res.args, quote delete_upper(m[$(QuoteNode(ex.args[3]))]) end)
            else
                push!(res.args, quote set_upper(m[$(QuoteNode(ex.args[3]))], $(_up)) end)
            end
            if _low isa Real && _up isa Real && _up < _low
                error(err)
            end
            return esc(res)
        else
            error(err)
        end
    else
        error(err)
    end
    return nothing
end

function parse_reactions(ex::Expr)
    function parse_term!(term, _sign, coef)
        if term isa Symbol
            coef[term] = _sign
        elseif term isa Expr && term.args[1] == :(*)
            if term.args[2] isa Real
                if term.args[3] isa Symbol
                    coef[term.args[3]] = _sign * term.args[2]
                else
                    error("Component `$(term.args[3])` is not a symbol")
                end
            else
                error("Stoichiometric coefficient `$(term.args[2])` is not a number")
            end
        end
    end
    function parse_terms!(terms, side, coef)
        _sign = (side == :left ? -1 : 1)
        if terms isa Expr
            if terms.args[1] == :(+)
                foreach(terms.args[2:end]) do t parse_term!(t, _sign, coef) end
            elseif terms.args[1] == :(*)
                parse_term!(terms, _sign, coef)
            else
                err_tok = (terms.head == :call ? terms.args[1] : terms.head)
                if err_tok == :(=>)
                    error("Only one instance of the `=>` operator is allowed")
                else
                    error("Invalid operator `$err_tok`")
                end
            end
        else
            parse_term!(terms, _sign, coef)
        end
    end
    function parse_reaction(ex::Expr)
        coef = OrderedDict{Symbol, Real}()
        if ex.args[1] == :(=>)
            parse_terms!(ex.args[2], :left, coef)
            parse_terms!(ex.args[3], :right, coef)
        else
            err_tok = (ex.head isa Symbol ? ex.head : ex.args[1])
            error("Looking for `=>`, found `$err_tok` instead")
        end
        return coef
    end

    coefs = Vector{OrderedDict{Symbol, Real}}()
    if ex.head == :block
        for rx in ex.args
            if rx isa Expr
                push!(coefs, parse_reaction(rx))
            end
        end        
    elseif ex.head == :call
        push!(coefs, parse_reaction(ex))
    end
    return coefs
end

function transform_variable_exprs(ex)
    ex isa Expr || return
    if ex.head == :block
        blk = Expr(:block)
        for subex in ex.args
            if subex isa Expr
                push!(blk.args, transform_variable_expr(subex))
            end
        end
        return blk
    else
        return transform_variable_expr(ex)
    end
end

"""
    @stoic(expr)

Create an array of OrderedDicts, each of which encodes the stoichiometric coefficients of a single reaction.
The `expr` should be one or more balanced chemical reactions written in the form:

    aA + bB + ... => cC + dD

where `a`, `b`, etc. are real or integer stoichiometric coefficients, and `A`, `B`, etc. are species.
The returned array looks like:

    coef[1] = OrderedDict(:A => -a, :B => -b, :C => c, :D => d)   # coefficients for reaction 1
    coef[2] = OrderedDict(:A => -a, :B => -b, :C => c, :D => d)   # coefficients for reaction 2
      ...

Note that the coefficients in the OrderedDict are positive for products and negative for reactants.

### Examples
```julia
julia> coef = @stoic A + B => C
julia> coef
1-element Vector{OrderedCollections.OrderedDict{Symbol, Real}}:
 OrderedCollections.OrderedDict(:A => -1, :B => -1, :C => 1)

julia> @stoic begin
        A + B => C
      2A + 3B => 4C + D
    end
2-element Vector{OrderedCollections.OrderedDict{Symbol, Real}}:
OrderedCollections.OrderedDict(:A => -1, :B => -1, :C => 1)
OrderedCollections.OrderedDict(:A => -2, :B => -3, :C => 4, :D => 1)
```

See also [`StoicReactor`](@ref)
"""
macro stoic(ex)
    (ex isa Expr) || return
    return parse_reactions(ex)
end

"""
    @solve

Run the function `optimize!(m)` where `m` is a JuMP model in the current scope.
"""
macro solve()
    return esc(:(optimize!(m)))
end

"""
    @set expr

Do the operations described by `expr` (see examples below). `@set` $m_note.

### Examples
Set the value of the variable `var` to 100.0. If the variable is fixed,
set the fixed value to 100.0, otherwise set the start value to 100.0:
```julia
julia> @set var = 100.0
```

Connect variables `var1` and `var2`. This adds the equation `var1 == var2` to the model. At
least one of the variables must be fixed:
```julia
julia> @set var1 = var2
```

Set the upper bound on `var` to 100.0:
```julia
julia> @set var <= 100.0
```

Set the upper bound on `var` to 100.0:
```julia
julia> @set var <= 100.0
```

Delete the upper bound on `var`:
```julia
julia> @set var <= Inf
```

Delete the lower bound on `var`:
```julia
julia> @set -Inf <= var
```

Delete the lower and upper bounds on `var`:
```julia
julia> @set -Inf <= var <= Inf
```

Flip the specs on `var1` and `var2`. If `var1` is fixed and `var2` is free, free `var1` and fix `var2`.
If `var1` is free and `var2` is fixed, fix `var1` and free `var2`. If both variables have the same spec, do nothing:
```julia
julia> @set var1 ~ var2
```

Flip the spec on `var`. If `var` is fixed, free `var`. If `var` is free, fix `var`:
```julia
julia> @set ~var
```

Fix `var`:
```julia
julia> @set +var
```

Free `var`:
```julia
julia> @set -var
```

Set the lower and upper bounds on `var`:
```julia
julia> @set 1.0 < var < 2.0
julia> @set 2.0 > var > 1.0
```

Combine several operations in a `begin`/`end` block:
```julia
julia> @set begin
    x > 1.0
    y < 2.0
    1.0 < z < 2.0
    p ~ q
end
```

See also [`@values`](@ref), [`@specs`](@ref), [`@bounds`](@ref)
"""
macro set(ex)
    return transform_variable_exprs(ex)
end

"""
    @values expr

Alias for `@set`.

See also [`@set`](@ref), [`@specs`](@ref), [`@bounds`](@ref)
"""
macro values(ex)
    return transform_variable_exprs(ex)
end

"""
    @specs expr

Alias for `@set`.

See also [`@set`](@ref), [`@values`](@ref), [`@bounds`](@ref)

"""
macro specs(ex)
    return transform_variable_exprs(ex)
end

"""
    @bounds expr

Alias for `@set`.

See also [`@set`](@ref), [`@values`](@ref), [`@specs`](@ref)

"""
macro bounds(ex)
    return transform_variable_exprs(ex)
end

"""
    @components symbols...

Create and return a component set made up of one or more unquoted `symbols`.

### Examples
```julia
julia> comps1 = @components h2 ch4 c2h4
OrderedSet{Symbol} with 3 elements:
  :h2
  :ch4
  :c2h4

julia> comps1 = @components(h2, ch4, c2h4)
OrderedSet{Symbol} with 3 elements:
  :h2
  :ch4
  :c2h4
```

See also [`Stream`](@ref)
"""
macro components(syms...)
    return esc(:(OrderedSet($syms)))
end

function make_stream(name, comps, ex_basis=:(basis=FRAC))
    err = "Last argument of `@stream` must be `(basis=FRAC)` or `(basis=FLOW)`"
    if name isa Symbol && comps isa Symbol && ex_basis isa Expr && ex_basis.head == :(=)
        if ex_basis.args[2] == :FRAC
            return esc(:(Stream($(QuoteNode(name)), fs, $comps)))
        elseif ex_basis.args[2] == :FLOW
            return esc(:(Stream($(QuoteNode(name)), fs, $comps, basis=FLOW)))
        else
            error(err)
        end
    else
        error(err)
    end
end

"""
    @stream(name, comps, ex_basis=:(basis=FRAC))

Create a Stream named `name` with component set `comps` and `StreamBasis` `basis`.
`@stream` $m_fs_note.

### Examples
```julia
julia> feed = @stream(feed, comps_feed)

julia> prod = @stream(feed, comps_prod, basis=FLOW)
```

See also [`Stream`](@ref), [`@streams`](@ref), [`StreamBasis`](@ref)
"""
macro stream(name, comps, ex_basis=:(basis=FRAC))
    return make_stream(name, comps, ex_basis)
end

"""
    @streams(blk)

Create a tuple of Streams using the expressions in the `begin`/`end` block `blk`.
`@streams` $m_fs_note.

### Examples
```julia
julia> (ethfd, profd, cgcout) = @streams(begin
    ethfd, comps_ethfeed
    profd, comps_profeed
    cgcout, comps_cg, (basis=FLOW)
end)
```

See also [`Stream`](@ref), [`@stream`](@ref), [`StreamBasis`](@ref)
"""
macro streams(blk)
    (blk isa Expr && blk.head == :block) || error("Argument to `@streams` must be a begin/end block")
    new_blk = Expr(:tuple)
    for ex in blk.args
        if ex isa LineNumberNode continue end
        if ex isa Expr && ex.head == :tuple
            push!(new_blk.args, make_stream(ex.args...))
        else
            error("Invalid expression `$ex` in `@streams` begin/end block")
        end
    end
    return new_blk
end

"""
    @block(name, kind, args...)

Create a block named `name` with type `kind` in flowsheet `fs` using arguments `args...`. The value of
`kind` must be a subtype of `AbstractBlock`. `@block` $m_fs_note.

### Examples
```julia
julia> ctdc1 = @block(ctdc1, Separator, arxo, [offgas, dc1bt])

julia> arx = @block(arx, StoicReactor, dc2oh, arxo, arx_stoic, arx_mw, arx_conv)
```

See also [`Mixer`](@ref), [`Splitter`](@ref), [`Separator`](@ref), [`YieldReactor`](@ref),
    [`MultiYieldReactor`](@ref), [`StoicReactor`](@ref)
"""
macro block(name, kind, args...)
    if name isa Symbol && kind isa Symbol
        return esc(:($kind(m, $(QuoteNode(name)), fs, $(args...))))
    else
        error("First two arguments to `@block` must be symbols")
    end
end