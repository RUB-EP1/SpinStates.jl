using Documenter
using SpinStates

makedocs(;
    sitename = "SpinStates.jl",
    modules = [SpinStates],
    authors = "Mikhail Mikhasenko",
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Spin orientation" => "spin_expectation.md",
            "Rotations" => "rotations.md",
            "Boosts" => "boosts.md",
            "Wigner rotation" => "wigner_rotation.md",
            "Dirac spinors" => "dirac_spinors.md",
        ],
        "API reference" => "api.md",
    ],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true"),
)

deploydocs(; repo = "github.com/RUB-EP1/SpinStates.jl.git", push_preview = true)
