using HCBE
using Documenter

makedocs(;
    modules=[HCBE],
    authors="Inazuma110 <c011703534@edu.teu.ac.jp> and contributors",
    repo="https://github.com/Inazuma110/HCBE.jl/blob/{commit}{path}#L{line}",
    sitename="HCBE.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Inazuma110.github.io/HCBE.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Inazuma110/HCBE.jl",
)
