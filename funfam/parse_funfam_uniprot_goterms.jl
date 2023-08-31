#=
Parse the `funfam_uniprot_goterms` produced by querying the Proteins API with the members
of FunFams.
=#

using CSV

function goterms_associated_with_funfams(funfam_dir; include_uniprotkb_kw_iea=false)
    associations = asyncmap(readdir(funfam_dir), ntasks=200) do file
        process_file(joinpath(funfam_dir, file), include_uniprotkb_kw_iea=include_uniprotkb_kw_iea)
    end
    return vcat(collect.(associations)...)
end

function process_file(file; include_uniprotkb_kw_iea=false)
    ff = replace(replace(basename(file), ".relatives.go"=>""), ".FF."=>"/FF/")
    associations = Set{Tuple{String,String}}()
    for line in eachline(file)
        l = split(line)
        ec = l[3]
        goterm = l[2]
        if ec == "IEA"
            if include_uniprotkb_kw_iea
                if l[4] != "UniProtKB-KW"
                    continue
                end
            else
                continue
            end
        end
        push!(associations, (ff, goterm))
    end
    return associations
end

# NB original data stored at $POMBEAGEINGGENES/Data/funfam/funfam_uniprot_goterms.tar.gz
funfam_dir = joinpath(ENV["LFS"], "pombage", "funfam_uniprot_goterms")

# excluding IEA terms
associations = goterms_associated_with_funfams(funfam_dir)

CSV.write(
    joinpath(ENV["POMBEAGEINGGENES"], "Data", "funfam", "funfam_uniprot_goterms_parsed.csv"),
    associations,
    writeheader=false)

# including IEA terms
associations = goterms_associated_with_funfams(funfam_dir; include_uniprotkb_kw_iea=true)

CSV.write(
    joinpath(ENV["POMBEAGEINGGENES"], "Data", "funfam", "funfam_uniprot_goterms_parsed_uniprotkb_kw_iea.csv"),
    associations,
    writeheader=false)
