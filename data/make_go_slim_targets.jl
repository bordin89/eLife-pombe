#=
Make targets from GO Slim annotations for function prediction tasks.

GO annotations:
    ftp://ftp.pombase.org/pombe/annotations/Gene_ontology/gene_association.pombase.gz
    http://www.geneontology.org/page/guide-go-evidence-codes

GO Slim:
    https://www.pombase.org/browse-curation/fission-yeast-go-slim-terms
    http://www.geneontology.org/page/go-subset-guide
=#

using PombeAgeingGenes, PombeAgeingGenes.GeneOntology, DataFrames, CSVFiles

ontology = load(GO)
annotations = unique!(load(GOAnnotations))

s2d = slim2descendants(ontology)
d2s = descendants2slim(s2d)
slimanddescendants = Set(vcat(keys(s2d)..., values(s2d)...))

filter!(x->x.go ∈ slimanddescendants, annotations)

function maketargets(annotations)
    slims = map(x->x.go in GO_SLIM_TERMS ? Set([x.go]) : d2s[x.go], annotations)

    df = DataFrame(id=String[], go=String[])

    for i = 1:length(annotations)
        r = annotations[i]

        for t = slims[i]
            push!(df, (r.id, t))
        end
    end

    unique!(df)

    df[:go] = map(x->replace(x, "GO:"=>"GO"), df[:go])

    df[:value] = 1
    df = unstack(df, :id, :go, :value)

    for col = names(df)
        df[col] = coalesce.(df[col], 0)
    end

    df
end

targets = maketargets(annotations)

save("$(ENV["POMBEAGEINGGENES"])/data/goslim_targets.csv", targets)
