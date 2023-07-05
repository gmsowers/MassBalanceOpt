push!(LOAD_PATH,"../src/")
using Documenter, MassBalanceOpt, JuMP
makedocs(
    sitename="MassBalanceOpt.jl",
    format = Documenter.HTML(edit_link = "main"),
    doctest = false,
    # modules=MassBalanceOpt,
    pages=["Introduction" => "index.md",
           "Examples"     => "Examples.md",
           "API"          => "MassBalanceOpt.md"]
    )