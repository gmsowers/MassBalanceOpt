# MassBalanceOpt.jl

A package for formulating mass balance models of chemical processes and solving them with [JuMP](https://jump.dev/).

## Flowsheets, Streams, and Blocks
```@docs
AbstractBlock
Flowsheet
@components
Stream
StreamBasis
is_frac
is_flow
@stream
@streams
copy_stream
copy_streams
make_stream_vars!
set_stream_var_specs!
Mixer
Splitter
Separator
@stoic
StoicReactor
YieldReactor
MultiYieldReactor
@block
@Block_fields
@Block_init
@Block_finish
make_mass_flow_vars!
set_start_values(::AbstractBlock; ::Bool)
```

## Variables and Equations
```@docs
make_var!
fix(::VariableRef)
free
flip
connect
set_value(::VariableRef, ::Real)
set_lower
set_upper
delete_lower
delete_upper
@set
@values
@specs
@bounds
make_eq!
```

## Printing and Output
```@docs
print_vars
write_vars
print_fixed
print_free
print_bounds
print_active
print_eqs
print_model
```

## Solving Models
```@docs
@solve
eval_obj
```

## Index
```@index
```