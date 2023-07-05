# MassBalanceOpt.jl

This is a package for formulating mass balance models of chemical processes and solving them with [JuMP](https://jump.dev/).
The user can build a model from a handful of unit operations ([`Mixer`](@ref), [`Splitter`](@ref), [`Separator`](@ref), 
[`StoicReactor`](@ref), [`YieldReactor`](@ref), [`MultiYieldReactor`](@ref)) connected by `Streams` and organized
into `Flowsheets`. The package also provides functions and macros for specifying and working with the model.
Familiarity with JuMP is a prerequisite for using this package.

The audience for the package is users of process simulators, particularly those that offer an equation-oriented solution mode, who want to solve optimization models of processes for which energy balances are either unneccesary or so straightforward that they can be formulated directly in JuMP. 

JuMP is a general purpose algebraic modeling language. The basic objects in a JuMP model are variables, variable bounds, constraints (referred to here as *equations*), objective functions, and solvers. MassBalanceOpt provides an additional layer of objects specific to process flowsheeting. Those objects are:

## Streams

A [`Stream`](@ref) represents a flow of material in a process. A stream contains a total mass flow rate and either a set of mass fractions (if the [`StreamBasis`](@ref) `= FRAC`, the default), or component mass flow rates (if the `StreamBasis = FLOW`). A stream's components are specified as an `OrderedSet` of symbols, usually created with the [`@components`](@ref) macro. Streams flow into and out of *blocks*, which are sets of equations that represent mass balance operations. By themselves, streams don't create or contain variables; the blocks do that. Streams can be created with the [`@stream`](@ref) or [`@streams`](@ref) macros, or by calling the `Stream` constructor directly.

## Blocks
A block is a subtype of [`AbstractBlock`](@ref). A block has inlet and outlet streams, and may have other input variables like split fractions or stoichiometric coefficients. A block contains the variables and equations that model a particular unit operation. Blocks can be created with the [`@block`](@ref) macro, or by calling the block's constructor. A block is created with some of its variables fixed to default values, so that the block is a self-contained model with the same number of equations and free variables. The inlet stream variables are initially fixed, so the blocks start out in a disconnected state. The function [`connect`](@ref) is used to connect a stream that flows out of one block and into another. This is a common practice in an equation-based process simulator.

## Flowsheets
A [`Flowsheet`](@ref) is a container for blocks and streams. Every model must have at least one flowsheet. Flowsheets can be created as children of an existing flowsheet, forming a tree structure. A flowsheet has a Julia symbol for a name. Flowsheet names are embedded in variable and equation names, so flowsheets provide namespaces that allow reuse of block and stream names in different flowsheets. A flowsheet created with the default constructor, `fs = Flowsheet()`, will have the name `:index` and has no parent. A child flowsheet can then be created with `fs_unit1 = Flowsheet(:unit1, fs)`. Variable names in the `:index` flowsheet have no prefix; variables in the `:unit1` flowsheet have the prefix `unit1_` prepended to their names. Thus you can create a stream named `feed` in the `:index` flowsheet and a stream named `feed` in the `:unit1` flowsheet without a variable name collision.

## Macros and utility functions
The package also supplies some macros and utility functions that make it easier to manipulate JuMP models. For example, instead of calling the JuMP function `set_lower_bound` to put a lower bound on a variable, you can use the macro [`@set`](@ref) (or if you want to be more explicit, [`@bounds`](@ref), which is simply an alias for `@set`), e.g.,
```julia
@set header_feed_mass < 100_000.0
```
sets the lower bound on variable `header_feed_mass`. You can do several assignments in one `@set` like this:
```julia
@set begin
    var1 = 1.0
    var2 > 0.0
    1.0 < var3 < 10.0
end
```
See [`@set`](@ref), [`@values`](@ref), [`@bounds`](@ref), [`@specs`](@ref).

You can print the variables in a block by calling
```julia
print_vars(block_name)
```
to get a table showing the solution value, bounds, start value, and an indication telling you if the variable is fixed or free. You can use the function
```julia
print_fixed(model_name)
```
to print only the fixed variables in the model. See [Printing and Output](@ref).

The package does not supply any functions or macros for creating or manipulating objective functions (except for a single function [`eval_obj`](@ref)).
