using Documenter, MORFEInvariantManifold

makedocs(sitename="MORFEInvariantManifold")

makedocs(
    sitename = "Model Order Reduction in Julia",
    format = Documenter.HTML(),
    authors = "Andrea Opreni, Alessandra Vizzaccaro",
    modules = [MORFEInvariantManifold],
    pages = Any[
        "Home" => "index.md",
        "Overview" => "overview.md",
        "Tutorials" => "tutorials.md",
        "Supported Formats" => "supported_formats.md",
        "Integrators" => "integrators.md",
        "Library" => "library_all.md",
        "Authors" => "authors.md",
        "References" => "literature.md"
    ]
)