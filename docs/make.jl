using Documenter, GridapHybrid

makedocs(;
    modules=[GridapHybrid],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/santiagobadia/GridapHybrid.jl/blob/{commit}{path}#L{line}",
    sitename="GridapHybrid.jl",
    authors="Santiago Badia <santiago.badia@monash.edu>",
    assets=String[],
)

deploydocs(;
    repo="github.com/santiagobadia/GridapHybrid.jl",
)
