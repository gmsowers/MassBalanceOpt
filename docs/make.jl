push!(LOAD_PATH,"../src/")
using Documenter, MassBalanceOpt, JuMP
makedocs(
    sitename="MassBalanceOpt.jl",
    # format = Documenter.HTML(prettyurls = false, disable_git = true, edit_link = nothing),
    doctest = false,
    # modules=MassBalanceOpt,
    pages=["Introduction" => "index.md",
           "Examples"     => "Examples.md",
           "API"          => "MassBalanceOpt.md"]
    )