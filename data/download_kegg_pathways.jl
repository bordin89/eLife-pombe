using DelimitedFiles

dir = joinpath(ENV["POMBEAGEINGGENES"], "data", "kegg")
isdir(dir) || mkpath(dir)

organism = "spo"

fp = joinpath(dir, "kegg_pathways.tsv")
isfile(fp) || download("http://rest.kegg.jp/list/pathway/" * organism, fp)

for id = readdlm(fp, '\t')[:,1]
    download("http://togows.dbcls.jp/entry/pathway/$(id)/genes.json",
             joinpath(dir, id[6:end]*".json"))
end
