push!(LOAD_PATH,"../src/")
using Documenter, MassBalanceOpt, JuMP

if isempty(ARGS)
    _format = Documenter.HTML(edit_link = "main")
elseif ARGS[1] == "pdf"
    _format = Documenter.LaTeX(platform = "native")
end

makedocs(
    sitename = "MassBalanceOpt.jl",
    format = _format,
    doctest = false,
    pages=["Introduction" => "index.md",
           "Examples"     => "Examples.md",
           "API"          => "MassBalanceOpt.md"]
    )