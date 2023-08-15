import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

clusters = [
    "A0119",
    "A0085",
    "A0151",
    "A3158",
    "A3266",
    "A3376",
    "A3391_A3395",
    "A3667",
    "A3716",
    "A2399",
]

df = pd.read_parquet(
    "results/MasterCatalogues/HectorClusters/HectorClusters_master_target_selected.parquet"
)

df_with_foreground = pd.read_parquet(
    "results/MasterCatalogues/HectorClusters/HectorClusters_target_selected_including_foreground.parquet"
)

# Add a column for the cluster each galaxy belongs to


# Add the cluster each one belongs to
region_information = pd.read_csv(
    "resources/RegionInformation/all_regions.csv", index_col="name"
)
region_information = region_information.loc[clusters]

for name, row in region_information.iterrows():

    region_mask = (
        (df.RA > row.min_RA)
        & (df.RA < row.max_RA)
        & (df.DEC > row.min_DEC)
        & (df.DEC < row.max_DEC)
    )
    df.loc[region_mask, "Cluster"] = name

    region_mask = (
        (df_with_foreground.RA > row.min_RA)
        & (df_with_foreground.RA < row.max_RA)
        & (df_with_foreground.DEC > row.min_DEC)
        & (df_with_foreground.DEC < row.max_DEC)
    )
    df_with_foreground.loc[region_mask, "Cluster"] = name

df = df.set_index("Cluster")
df_with_foreground = df_with_foreground.set_index("Cluster")


def make_cluster_plot(df, suptitle):

    fig, axs = plt.subplots(
        ncols=5, nrows=2, constrained_layout=True, figsize=(11.88, 5.48)
    )
    text = f"""# Average densities ({suptitle}):\n\n"""
    for name, ax in zip(clusters, axs.ravel()):

        target_counts, xedges, yedges, _ = stats.binned_statistic_2d(
            x=df.loc[name, "RA"],
            y=df.loc[name, "DEC"],
            values=df.loc[name, "RA"],
            statistic="count",
            bins=10,
        )

        bin_areas = np.ediff1d(xedges) * np.ediff1d(yedges)
        target_density = target_counts / bin_areas

        text += f"* {name}: {target_density.mean():.1f} galaxies per square degree ({target_density.mean() * np.pi:.1f} per 2dF FoV)\n"

        im = ax.pcolormesh(
            xedges, yedges, target_density.T * np.pi, vmin=0, vmax=200, cmap="plasma"
        )
        ax.set_title(name)
    fig.colorbar(im, label="N galaxies per 2dF FoV", ax=axs[:, -1:])
    fig.suptitle(suptitle)
    text += "\n"
    return fig, axs, text


fig_nf, axs_nf, nf_text = make_cluster_plot(df, "No Foreground galaxies")
fig_wf, axs_wf, wf_text = make_cluster_plot(
    df_with_foreground, "With Foreground galaxies"
)

print(nf_text)
print(wf_text)

with open("misc_analysis_results/no_foreground_galaxy_info.md", "w") as f:
    f.write(nf_text)

with open("misc_analysis_results/with_foreground_galaxy_info.md", "w") as f:
    f.write(wf_text)

# Save the plots and our info
fig_nf.savefig(
    "misc_analysis_results/cluster_target_density_no_foreground_galaxies.png"
)
fig_wf.savefig(
    "misc_analysis_results/cluster_target_density_with_foreground_galaxies.png"
)
