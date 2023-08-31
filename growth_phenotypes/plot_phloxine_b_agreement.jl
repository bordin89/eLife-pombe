#=
Could show this plot in the thesis chapter to explain why I pooled conditions w/ and w/o
phloxine B, but the correlation (~0.4) is not great.
=#

using PombeAgeingGenes
using DataFrames
using DataFramesMeta
using Statistics

df = load(GrowthPhenotypesNoOutliers)

df.phlox = parse.(Bool, df.phlox)
df = @where(df, :condition .== "YES_32")

df_with = @where(df, :phlox .== true)
df_wout = @where(df, :phlox .== false)


df_with = by(df_with, :id, size = :size => median)
df_wout = by(df_wout, :id, size = :size => median)

df = join(df_with, df_wout, on = :id, makeunique = true)

cor(df.size, df.size_1)

using Plots

scatter(
    df.size,
    df.size_1,
)

plot!(
    [-2, 2],
    [-2, 2],
)
