import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.spatial as sp, scipy.cluster.hierarchy as hc

sns.set(font_scale=.5)

ageingterms = {
    "FYPO0004344": [
        "SPBC1861.02", "SPBC2G2.06c", "SPBC691.03c", "SPBC685.04c", "SPBC1D7.03",
        "SPBC428.08c", "SPAC4G9.11c", "SPAC23A1.06c", "SPAC1805.07c", "SPBC947.05c",
        "SPAC17A5.16", "SPCC1753.02c", "SPBC32F12.03c", "SPAC1687.15", "SPCC4G3.09c",
        "SPAC1834.04", "SPBC947.09", "SPBC19C7.01", "SPAC694.06c", "SPAC8F11.03",
        "SPAC1F3.09", "SPAC18B11.04", "SPCC663.06c", "SPBC543.07", "SPAC1F8.06",
        "SPBC106.10", "SPAC1A6.04c", "SPBC365.20c", "SPAC890.03", "SPBC8E4.02c",
        "SPAC16.01", "SPBC646.13", "SPCC4B3.12", "SPBC21C3.18", "SPAC20H4.03c",
        "SPBC17G9.09", "SPACUNK4.16c", "SPCC1020.08", "SPCC188.08c", "SPAC1783.02c",
        "SPAC13F5.04c", "SPAC23D3.03c", "SPAC8E11.05c", "SPBC1198.07c", "SPBC14C8.15",
        "SPBC18H10.18c", "SPBC1921.04c", "SPBC30D10.09c", "SPBC4B4.12c", "SPCC306.11",
        "SPCC320.03", "SPCC594.01", "SPCC594.02c", "SPCC794.03",
    ],
    "FYPO0001309": [
        "SPAC21E11.04", "SPCC16A11.08", "SPBC21C3.08c", "SPBC1D7.03", "SPCC1753.02c",
        "SPBC16D10.08c", "SPBC3H7.03c", "SPBC16E9.13", "SPBP4H10.11c", "SPAC17C9.02c",
        "SPAC821.07c", "SPBP35G2.11c", "SPAC806.07", "SPCC188.02", "SPCC16C4.11",
        "SPBC725.11c", "SPAC23C11.08", "SPBC3B8.02", "SPBC11B10.10c", "SPBC106.10",
        "SPAC26F1.10c", "SPBC1198.11c", "SPAC1B9.02c", "SPAC22E12.14c", "SPAC16E8.01",
        "SPCC162.12", "SPBP23A10.16", "SPBC30D10.10c", "SPAC1399.04c", "SPAC3A12.09c",
        "SPBC16E9.14c", "SPAC323.03c", "SPAC3H1.08c", "SPBP4H10.16c", "SPRRNA.47",
    ],
}

in_dir = "/Users/harry/Dropbox/UCL/pombage/data"
out_dir = "/Users/harry/Dropbox/UCL/pombage/Scripts/growth_phenotypes"

# conditions cor
# df = pd.read_csv(f"{in_dir}/Jan2019_BBSRC_results_no_outliers_wideform.csv")
df = pd.read_csv(f"{in_dir}/Oct2019_BBSRC_results_wideform.csv")
c = df.iloc[:,1:].corr()
linkage = hc.linkage(sp.distance.squareform(1 - c), method='average')

fig = sns.clustermap(
    c,
    row_linkage=linkage,
    col_linkage=linkage,
    cmap="RdBu_r", vmax=1., vmin=-1.,
    yticklabels=True, xticklabels=True,
    cbar_kws={"label": "PCC", "ticks": [-1., -.5, 0, .5, 1.]},
    figsize=(12,12),
    )

# fig = sns.clustermap(c,
#     metric="euclidean",
#     method="average",
#     cmap="RdBu_r", vmax=1., vmin=-1.,
#     yticklabels=True, xticklabels=True,
#     cbar_kws={"label": "PCC"},
#     figsize=(12,12),
#     )

fig.savefig(f"{out_dir}/clustermap_conditions.pdf")

# cor 2
c2 = c.corr() # np.corrcoef(c)
linkage = hc.linkage(sp.distance.squareform(1 - c2), method='average')

fig = sns.clustermap(
    c2,
    row_linkage=linkage,
    col_linkage=linkage,
    cmap="RdBu_r", vmax=1., vmin=-1.,
    yticklabels=True, xticklabels=True,
    cbar_kws={"label": "PCC", "ticks": [-1., -.5, 0, .5, 1.]},
    figsize=(12,12),
    )

# fig = sns.clustermap(c2,
#     metric="euclidean",
#     method="average",
#     cmap="RdBu_r", vmax=1., vmin=-1.,
#     yticklabels=True, xticklabels=True,
#     cbar_kws={"label": "PCC"},
#     figsize=(12,12),
#     )

fig.savefig(f"{out_dir}/clustermap_conditions_cor2.pdf")


# ids cor
# df = pd.read_csv(f"{in_dir}/Jan2019_BBSRC_results_no_outliers_wideform.csv")
df = pd.read_csv(f"{in_dir}/Oct2019_BBSRC_results_wideform.csv")
df = df.rename(columns={"id":"Genes"})
df = df.set_index("Genes")

# df = df.T
# df.columns = df.iloc[0]
# df = df.drop(df.index[0])

c = df.T.corr()
# c = tmp.astype(float).corr()

# cols = ["g" if (i in ageingterms["FYPO0004344"] or i in ageingterms["FYPO0001309"]) else "w" for i in c.index.values]

# c = np.corrcoef(df.astype(float).values.T)
linkage = hc.linkage(sp.distance.squareform(1 - c), method='average')

fig = sns.clustermap(
    c,
    row_linkage=linkage,
    col_linkage=linkage,
    cmap="RdBu_r", vmax=1., vmin=-1.,
    yticklabels=False, xticklabels=False,
    cbar_kws={"label": "PCC", "ticks": [-1., -.5, 0, .5, 1.]},
    figsize=(12,12),
    )

# fig = sns.clustermap(c,
#     metric="euclidean",
#     method="average",
#     cmap="RdBu_r", vmax=1., vmin=-1.,
#     yticklabels=False, xticklabels=False,
#     cbar_kws={"label": "PCC"},
#     # row_colors=cols, col_colors=cols,
#     figsize=(12,12),
#     )

fig.savefig(f"{out_dir}/clustermap_ids.png", dpi=150)

# cor2

c2 = c.corr()

# c2 = np.corrcoef(c)
linkage = hc.linkage(sp.distance.squareform(1 - c2), method='average')

fig = sns.clustermap(
    c2,
    row_linkage=linkage,
    col_linkage=linkage,
    cmap="RdBu_r", vmax=1., vmin=-1.,
    yticklabels=False, xticklabels=False,
    cbar_kws={"label": "PCC", "ticks": [-1., -.5, 0, .5, 1.]},
    figsize=(12,12),
)

fig.savefig(f"{out_dir}/clustermap_ids_cor2.png", dpi=150)

# ids x conditions
# df = pd.read_csv(f"{in_dir}/Jan2019_BBSRC_results_no_outliers_wideform.csv")
df = pd.read_csv(f"{in_dir}/Oct2019_BBSRC_results_wideform.csv")
df = df.rename(columns={"id":"Strains"})
x = df.set_index(df["Strains"]).drop("Strains", axis=1)
x.columns.name = "Conditions"

# cols = ["r" if (i in ageingterms["FYPO0004344"] or i in ageingterms["FYPO0001309"]) else "w" for i in x.index.values]

# Add column colouring for YES and EMM
c = sns.color_palette("colorblind")
col_colors = [c[1] if xi.startswith("YES") else c[2] for xi in x.columns]

fig = sns.clustermap(x,
    metric="euclidean",
    method="average",
    cmap="RdBu_r", vmin=-2., vmax=2.,
    yticklabels=False, xticklabels=False,
    col_colors=col_colors,
    # row_colors=cols,
    cbar_kws={"label": "Colony size", "ticks": [-2, -1, 0, 1, 2]},
    # figsize=(8,12),
)

fig.savefig(f"{out_dir}/clustermap_ids_conditions.png", dpi=150)


# ids x conditions trigitised
sns.set(font_scale=1)

df = pd.read_csv(f"{in_dir}/Oct2019_BBSRC_results_wideform.csv")
df = df.rename(columns={"id":"Strains"})
x = df.set_index(df["Strains"]).drop("Strains", axis=1)
x.columns.name = "Conditions"

# Trigitise
threshold = 0.32 # log2(1.25) == 0.32
x[x < -threshold] = -1
x[x > threshold] = 1
x[(-threshold <= x) & (x <= threshold)] = 0

# Add column colouring for YES and EMM
c = sns.color_palette("colorblind")
col_colors = [c[2] if xi.startswith("YES") else c[4] for xi in x.columns]

c = sns.xkcd_palette(["windows blue", "amber"])
fig = sns.clustermap(
    x,
    metric="euclidean",
    method="ward",
    # cmap="RdBu_r", vmin=-2, vmax=2,
    cmap = [c[0], "w", c[1]],
    yticklabels=False, xticklabels=False,
    col_colors=col_colors,
    # row_colors=cols,
    cbar_kws={"ticks": [-1, 0, 1], "label": "sensitive       resistant"},
    figsize=(8,12),
)

fig.savefig(f"{out_dir}/clustermap_ids_conditions_trigitised.png", dpi=300)


df = pd.read_csv(f"{in_dir}/Jan2019_BBSRC_results_no_outliers_wideform_trigitised.csv")
df = df.rename(columns={"id":"Genes"})
x = df.set_index(df["Genes"]).drop("Genes", axis=1)
x.columns.name = "Growth conditions"

cols = ["g" if (i in ageingterms["FYPO0004344"] or i in ageingterms["FYPO0001309"]) else "w" for i in x.index.values]

fig = sns.clustermap(x,
    # foo
    cmap="RdBu_r", vmin=-1.5, vmax=1.5,
    yticklabels=False, xticklabels=False,
    row_colors=cols,
    cbar_kws={"label": "sensitive, resistant", "ticks": [-1, 0, 1]},
    # figsize=(8,12),
    )

fig.savefig(f"{out_dir}/clustermap_ids_conditions_trigitised.png", dpi=600)
