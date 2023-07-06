# MassBalanceOpt.jl

This is a Julia package for formulating mass balance models of chemical processes and solving them with [JuMP](https://jump.dev/). If you don't have basic knowledge of JuMP this package probably won't be of much help.

This package is unregistered, which means it's not in the Julia package registry. If you happen to discover it and want to use it, do this:

1. Start Julia, press `]` to get to the Pkg REPL, activate your desired environment, and add the package:
```julia
(test) pkg> add "https://github.com/gmsowers/MassBalanceOpt"
```

2. Add the JuMP, Ipopt, and OrderedCollections packages:
```julia
(test) pkg> add JuMP Ipopt OrderedCollections
```
Expect precompilation of JuMP to take a long time.

3. Optionally, run the package tests:
```julia
(test) pkg> test MassBalanceOpt
```
You should see `Testing MassBalanceOpt tests passed`.

4. Press the backspace key to exit the Pkg REPL, and import the packages:
```julia
julia> using MassBalanceOpt, JuMP, Ipopt, OrderedCollections
```

5. You should now be able to create a Model and Flowsheet:
```julia
julia> m = Model; fs = Flowsheet()
Flowsheet(
   name=index   )
```

Package [documentation](https://gmsowers.github.io/MassBalanceOpt/).
