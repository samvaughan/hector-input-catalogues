"""
We ended up making some changes to this target selection workflow
_after_ we'd already observed ~1000 galaxies. Some of the galaxies
we have observations for no longer make it into our final master catalogues, 
either because they were randomly not picked this time or because they
ended up being too close to the star cuts I've added.

This script adds the galaxies which we've already got observations
for back into the master catalogues. 
"""
import pandas as pd

smk = snakemake  # noqa

# Read in the catalogues
observed_galaxies_up_to_August2023 = pd.read_csv(
    smk.input.galaxies_observed_up_to_August_2023
)
master_catalogue = pd.read_parquet(smk.input.master_catalogue)

# Only select the appropriate galaxies to add in
observed_galaxies_in_region = observed_galaxies_up_to_August2023.loc[
    observed_galaxies_up_to_August2023["Field"] == "HectorClusters"
]
assert (
    len(observed_galaxies_in_region) > 0
), "Looks like we have no matching galaxies..."

# Merge the two catalogues and keep only things which were observed but not in the master catalogue
merged = pd.merge(
    master_catalogue["ID"],
    observed_galaxies_in_region,
    left_on="ID",
    right_on="ID",
    how="right",
    indicator=True,
    suffixes=("_x", None),
)
missing_galaxies = merged.loc[merged["_merge"] == "right_only"]
assert len(missing_galaxies) > 0, "Looks like we have no matching galaxies..."

# Now add the missing things back in
combined = pd.concat(
    (master_catalogue, missing_galaxies.drop(["Field", "Region", "_merge"], axis=1))
)
# Fix up some data types
combined["RS_member"] = combined["RS_member"].astype(bool)
combined["ClusterMember"] = combined["ClusterMember"].astype(bool)
combined["MW_analogue"] = combined["MW_analogue"].astype(bool)
combined["EO_wind_galaxy"] = combined["EO_wind_galaxy"].astype(bool)
# And save
combined.to_parquet(smk.output.master_catalogue_with_all_observations, index=False)
