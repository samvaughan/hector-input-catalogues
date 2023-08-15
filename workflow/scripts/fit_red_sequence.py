import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from cmdstanpy import CmdStanModel
import astropy.units as u
import calc_kcor

smk = snakemake
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

run_on_subset = True
# smk = snakemake
# full_df = pd.read_parquet(smk.input.initial_cluster_catalogue)
full_df = pd.read_parquet(smk.input.cluster_cat_no_RS)
non_members = full_df.loc[full_df.mem_flag == 0]
df = full_df.loc[full_df.mem_flag == 1]

# Calculate the absolute R-band magnitude, including a K correction
k_cors = calc_kcor.calc_kcor("r", df.z, "g - r", df.mag_g - df.mag_r)
abs_mag_r = (
    df.mag_r.values * u.mag - cosmo.distmod(df.z.values) - k_cors.values * u.mag
).value
df["Absolute_mag_R"] = abs_mag_r

# Centre these values around the mean and make the y values
mean_abs_mag = np.mean(df.Absolute_mag_R)
x_all = df.Absolute_mag_R - mean_abs_mag
y_all = df.mag_g - df.mag_r

# Get rid of some galaxies which are very red
x = x_all.loc[y_all < 1.05]
y = y_all.loc[y_all < 1.05]
N = len(x)

data = dict(N=N, x=x, y=y)
if run_on_subset:
    N_subset = 1000
    indices = x.sample(n=N_subset).index
    x_subset = x.loc[indices]
    y_subset = y.loc[indices]
    data = dict(N=N_subset, x=x_subset, y=y_subset)

# Fit the stan model
sm = CmdStanModel(stan_file="fit_red_sequence.stan")
fit = sm.sample(data=data)

# Get the results as a pandas dataframe
samples = fit.draws_pd()

# Make some straight lines
xx = np.linspace(-4, 4)
yy_blue = (
    xx[:, None] * np.array(samples["slopes[1]"])[None, :]
    + np.array(samples["intercepts[1]"])[None, :]
)  # Blue Cloud
yy_red = (
    xx[:, None] * np.array(samples["slopes[2]"])[None, :]
    + np.array(samples["intercepts[2]"])[None, :]
)  # red sequence

# Make the probabilities
Red_Sequence_member_probability = 1 - np.exp(
    samples.loc[:, "log_Pr[1]":f"log_Pr[{data['N']}]"].mean()
)

# Make some plots
# Firstly, the red sequence membership probability on the color/magnitude diagram
plt.style.use("publication")
fig, ax = plt.subplots(figsize=(14, 6))
im = ax.scatter(
    data["x"] + mean_abs_mag,
    data["y"],
    c=Red_Sequence_member_probability.values,
    alpha=0.2,
    cmap="RdYlBu_r",
    rasterized=True,
    s=100,
)

# Make the colorbar, without any alpha
cmap = mpl.cm.RdYlBu_r
norm = mpl.colors.Normalize(vmin=0, vmax=1)
fig.colorbar(
    mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
    ax=ax,
    label=r"$p(\mathrm{red\,\,sequence\,\,member})$",
)

# Add on some straight line fits
# ax.plot(xx + mean_abs_mag, yy_blue.mean(axis=1), c='k', alpha=1, linestyle='dotted', linewidth=1.5)
ax.plot(
    xx + mean_abs_mag,
    yy_red.mean(axis=1),
    c="k",
    alpha=1,
    linestyle="solid",
    linewidth=3.5,
)
ax.plot(
    xx + mean_abs_mag,
    yy_red.mean(axis=1) - 2 * samples["scatter[2]"].mean(),
    c="0.2",
    linestyle="dashed",
    linewidth=1.5,
)
ax.plot(
    xx + mean_abs_mag,
    yy_red.mean(axis=1) + 2 * samples["scatter[2]"].mean(),
    c="0.2",
    linestyle="dashed",
    linewidth=1.5,
)

ax.set_xlabel("Absolute $r$-band magnitude (mag)")
ax.set_ylabel("$g - r$ colour (mag)")
ax.set_ylim(0.0, 1.2)
fig.savefig("red_sequence_clusters.pdf", bbox_inches="tight")


# Save the table
Red_Sequence_member_probability.index = data["x"].index
# Using the fact that the index is now correct
df["Red_Sequence_member_probability"] = Red_Sequence_member_probability

# Tweak the fact that some very red things have a smaller probability of being a red sequence member
df.loc[df.mag_g - df.mag_r > 1.05, "Red_Sequence_member_probability"] = 1

df["Red_Sequence_member_from_mixture_model"] = (
    df["Red_Sequence_member_probability"] > 0.4
)

# Add the column where we take everything above the red sequence minus 2 sigma of scatter to be a red sequence member
df["Red_Sequence_member_2_sigma_scatter"] = data["y"] > (
    data["x"] * samples["slopes[2]"].mean()
    + samples["intercepts[2]"].mean()
    - 2 * samples["scatter[2]"].mean()
)
df["Red_Sequence_member_2_sigma_scatter"] = df[
    "Red_Sequence_member_2_sigma_scatter"
].astype(bool)

df.to_parquet(smk.output.final_cluster_catalogue, index=False)
