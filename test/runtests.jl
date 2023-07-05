using MassBalanceOpt
using JuMP, Ipopt
using OrderedCollections
using Test

testsets = [
    "set_start.jl",
    "macros.jl",
    "mixer.jl",
    "splitter.jl",
    "separator.jl",
    "yield_reactor.jl",
    "multi_yield_reactor.jl",
    "stoic_reactor.jl",
    "mixer+splitter.jl",
    "model.jl"
]

@testset verbose = true "MassBalanceOpt tests" begin
    for t in testsets
        include(t)
    end
end
nothing
