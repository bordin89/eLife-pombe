#=
Download GO taxon constraints
=#

dir = joinpath(ENV["POMBEAGEINGGENES"], "Data")

download(
    "https://raw.githubusercontent.com/geneontology/go-ontology/c13fd5123d798885420d65f3668c037953fe47db/src/taxon_constraints/only_in_taxon.tsv",
    joinpath(dir, "only_in_taxon.tsv")
    )

download(
    "https://raw.githubusercontent.com/geneontology/go-ontology/c13fd5123d798885420d65f3668c037953fe47db/src/taxon_constraints/never_in_taxon.tsv",
    joinpath(dir, "never_in_taxon.tsv")
    )
