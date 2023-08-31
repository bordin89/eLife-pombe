import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.spatial as sp, scipy.cluster.hierarchy as hc

sns.set(font_scale=1)

df = pd.read_csv(f"/Users/harry/Dropbox/UCL/pombage/data/network_embeddings/network_embeddings.csv")
df = df.rename(columns={"id":"Proteins"})
df = df.set_index("Proteins")

df = df.iloc[:, list(df.sum(axis=0) != 0)]

fig = sns.clustermap(
    df.iloc[:,1:],
    cmap="viridis",
    yticklabels=False, xticklabels=False,
    cbar_kws={"label": "Embedding value"},
    figsize=(12,12),
    )

fig.ax_heatmap.set_xlabel("Embedding dimensions")

fig.savefig(f"/Users/harry/Dropbox/UCL/pombage/Scripts/network_embeddings/network_embeddings.png", dpi=300)
