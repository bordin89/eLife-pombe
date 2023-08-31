
using Plots

plot(iterrange,
    mean.(results),
    fillrange=(maximum.(results), minimum.(results)),
    xticks=1:5:nconditions-1,
    lc=:black,
    legend=false,
    xlabel="# conditions",
    ylabel="AUCPR",
    fillalpha=.5,
    )
