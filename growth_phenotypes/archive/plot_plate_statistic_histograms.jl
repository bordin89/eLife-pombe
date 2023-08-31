#=
Plot histograms of plate statistics.
=#

cd(@__DIR__)

using PombeAgeingGenes, Statistics, Plots

plotlyjs()

df = load(GrowthPhenotypesNoOutliers)

platesummary = by(df, [:condition, :assay_plate, :scan_date, :repeat_number,
                       :library_version, :phlox]) do x
    y = x[:size]
    (mean = mean(y), var = var(y), cv = coefvar(y), snr = snr(y), missingcolonies = 1536 - size(x,1))
end

opts = (ylabel="Frequency (no. plates)", legend=false, bins=100)

histogram(platesummary[:mean], xlabel="Mean size"; opts...)
histogram(platesummary[:var], xlabel="Variance"; opts...)
histogram(platesummary[:cv], xlabel="CV"; opts...)
histogram(platesummary[:snr], xlabel="SNR"; opts...)
histogram(platesummary[:missingcolonies], xlabel="Missing colonies"; opts...)
