#=
Generate files containing the FunFams that have proteins annotated with a particular GO
term.
=#

using PombeAgeingGenes, OBOParse, Distributed, DelimitedFiles

addprocs(3)

function funfams_with_goterm(ontoloy, goterm, dir)
    goterms = Set(descendants(ontology, goterm))
    ffs = String[]
    dir = joinpath(dir, "funfam_uniprot_goterms")

    ffs = @distributed vcat for file = readdir(dir)
        search_file(dir, file, goterms)
    end

    string.(ffs[.!(isnothing.(ffs))])
end

@everywhere function search_file(dir, file, goterms)
    open(joinpath(dir, file), "r") do io
        for line = eachline(io)
            l = split(line)

            # TODO only include IEA from UniProtKB-KW?
            # if l[3] == "IEA"
            #     if l[4] != "UniProtKB-KW"
            #         continue
            #     end
            # end

            if l[2] in goterms
                s = split(file, '.')
                return join([join(s[1:4], '.'); s[5:6]], '/')
            end
        end
    end
end

# Load files
dir = joinpath(ENV["POMBEAGEINGGENES"], "data", "funfam")
file = "funfam_uniprot_goterms"
cd(dir)
run(`tar zxf $(joinpath(dir, file*".tar.gz"))`)

ontology = load(GeneOntology.GO)

# Test
goterm = "GO:0030036"
funfams_with_goterm(ontology, goterm, dir)

# Iterate over all GO Slim terms
for goterm = GeneOntology.GO_SLIM_TERMS
    ffs = funfams_with_goterm(ontology, goterm, dir)
    writedlm(joinpath(dir, "funfams_with_goterms", replace(goterm, ":"=>"")*".csv"), ffs)
end

# Finally
rm(joinpath(dir, file), recursive=true)
