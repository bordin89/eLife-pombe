#=
Analysis of GO term information content for priority unstudied predictioncs vs annotations already
in Pombase.
=#

using PombeAgeingGenes
using CSV
using DataFrames
using DataStructures
using HypothesisTests
using StatsPlots

priority_unstudied_genes_df = CSV.File(joinpath(ENV["POMBEAGEINGGENES"], "data", "priority_unstudied_genes.tsv")) |> DataFrame
priority_unstudied_genes = collect(priority_unstudied_genes_df[:,1])

df = CSV.File(joinpath(@__DIR__, "predictions_combined.csv")) |> DataFrame
df = @in(df, :gene_id, priority_unstudied_genes)

counter(df.pombase_evidence_code)

df_new = @where(df, ismissing.(:pombase_evidence_code))
df_old = @where(df, .!ismissing.(:pombase_evidence_code))

mean(df_old.information_content)
mean(df_new.information_content)

median(df_old.information_content)
median(df_new.information_content)

violin(
    [
        df_old.information_content,
        df_new.information_content,
    ],
    ylabel = "Information content",
    xticks = (
        1:2,
        ["Annotations", "Predictions"],
    ),
    legend = false,
    c = :grey80,
)

savefig(joinpath(@__DIR__, "information_content_annotations_vs_predictions_priority_unstudied.pdf"))

pvalue(
    MannWhitneyUTest(
        df_new.information_content,
        df_old.information_content,
    ),
    tail = :right,
)
