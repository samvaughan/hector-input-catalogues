import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 

df = pd.read_parquet("results/FinalInputCatalogues/HectorClusters_final_input.parquet")
clusters = np.array(['A0085', 'A0119', 'A0151', 'A2399', 'A3158', 'A3266', 'A3376', 'A3391/A3395', 'A3667', 'A3716'])

fig, axs = plt.subplots(ncols=5, nrows=2, figsize=(20, 8), constrained_layout=True)

for ax, cluster in zip(axs.ravel(), clusters):
    if cluster == 'A3391/A3395':
        ax.scatter(df.loc[df['NearClus'].isin(['A3391', 'A3395']), 'z'], df.loc[df.NearClus.isin(['A3391', 'A3395']), 'logMstarProxy_LS'], label='A3391/A33945')
    else:
        ax.scatter(df.loc[df['NearClus'] == cluster, 'z'], df.loc[df.NearClus == cluster, 'logMstarProxy_LS'], label=cluster)
    ax.set_title(cluster)

# Now add the Hector steps
for ax in axs.ravel():
    ax.plot([0.06, 0.095], [10.6, 10.6], c='k', linewidth=3.0)
    ax.plot([0.045, 0.06], [9.6, 9.6], c='k', linewidth=3.0)
    ax.plot([0.03, 0.045], [8.6, 8.6], c='k', linewidth=3.0)
    ax.plot([0.045, 0.045], [8.6, 9.6], c='k', linewidth=3.0)
    ax.plot([0.06, 0.06], [9.6, 10.6], c='k', linewidth=3.0)

fig.show()
