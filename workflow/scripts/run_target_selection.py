"""
Run the target selection on our input catalogues.

Run on the entire, combined master catalogues for WAVES South, WAVES North and the clusters. We then break up into smaller regions in the next script.

NOTE: For the clusters, we use the redshifts of the cluster BCG, **not** the measured redshifts of individual member galaxies. Matt has already analysed the membership status of the objects, and we've already thrown away non-members in prepare_cluster_catalogues.py
"""

from hop import pipeline
import matplotlib.pyplot as plt
import shutil
from pathlib import Path
import pandas as pd
import utils
import numpy as np

# import scipy.stats as stats

smk = snakemake  # noqa

config_filename = smk.input.target_selection_config_file

HP = pipeline.HectorPipe(config_filename=config_filename)

# Make sure to have a random generator with known seed for repeatable results
random_generator = np.random.default_rng(12345)
fig, axs = HP.run_target_selection(
    random_generator=random_generator, save=False, plot=True
)
fig.savefig(smk.output.target_selection_plot, bbox_inches="tight")
plt.close("all")

# If we're dealing with the clusters, add back in the individual galaxy redshifts
if smk.wildcards["master_region"] == "HectorClusters":
    # Read the master input catalogue
    master_input_catalogue = pd.read_parquet(smk.input.input_catalogue)

    # Make sure our indices are all integer types
    master_input_catalogue.index = pd.to_numeric(master_input_catalogue.index)
    HP.target_selection.selection_function_sparsely_sampled.index = pd.to_numeric(
        HP.target_selection.selection_function_sparsely_sampled.index
    )
    HP.target_selection.selection_function.index = pd.to_numeric(
        HP.target_selection.selection_function.index
    )

    # Now just set the ['z'] columns to be equal- pandas takes care of the matching
    HP.target_selection.selection_function["z"] = master_input_catalogue["z"]
    HP.target_selection.selection_function_sparsely_sampled["z"] = (
        master_input_catalogue["z"]
    )


# If we're in the clusters, make sure the BCGs have been selected


# Go through and add the galaxy prorities to the catalogues
# Save the sparsely sampled catalogue and the un-sampled one too
print("Adding the galaxy priorities...")
for target_df, output_name in zip(
    [
        HP.target_selection.selection_function,
        HP.target_selection.selection_function_sparsely_sampled,
    ],
    [
        smk.output.master_target_selected_catalogue_not_subsampled,
        smk.output.master_target_selected_catalogue,
    ],
):
    # First select MW analogues and wind galaxies
    mw_analogue_mask = utils.select_MW_analogues(target_df)
    eo_wind_mask = utils.select_edge_on_wind_galaxies(target_df)
    low_mass_mask = utils.low_mass_galaxy_mask(target_df)

    target_df["MW_analogue"] = mw_analogue_mask
    target_df["EO_wind_galaxy"] = eo_wind_mask

    # Now make the priority and number of obs columns
    target_df["priority"] = 4
    target_df["N_observations_to_complete"] = 1

    # If we're in the clusters, downweight things within 1 r200 or outside 2R200, especially if they're red
    if smk.wildcards["master_region"] == "HectorClusters":
        one_to_two_r200_mask = utils.select_one_to_two_r200(target_df)
        within_r200_mask = utils.select_within_r200(target_df)
        outside_2_r200 = utils.select_outside_2r200(target_df)

        target_df.loc[within_r200_mask & target_df["RS_member"], "priority"] = 1
        target_df.loc[outside_2_r200 & target_df["RS_member"], "priority"] = 1

        target_df.loc[within_r200_mask & ~target_df["RS_member"], "priority"] = 8
        target_df.loc[one_to_two_r200_mask & ~target_df["RS_member"], "priority"] = 8

    # Now make sure that the edge-on wind galaxies or MW analogues are high priority
    target_df.loc[mw_analogue_mask, "priority"] = 8
    target_df.loc[eo_wind_mask, "priority"] = 8

    # And the low-mass repeats
    target_df.loc[low_mass_mask, "priority"] = 8
    target_df.loc[low_mass_mask, "N_observations_to_complete"] = 3

    target_df.to_parquet(output_name)

print("\tDone!")
# Clean up our output folder
base_output_folder = Path(smk.output.master_target_selected_catalogue).parent
folders_to_delete = [
    "Allocation",
    "Configuration",
    "DistortionCorrected",
    "Fibres",
    "FinalOutputs",
    "Logs",
    "Plots",
    "Tiles",
]
for folder in folders_to_delete:
    path = base_output_folder / folder
    try:
        shutil.rmtree(path)
    except FileNotFoundError:
        pass

# # Run some stats on the WAVES catalogues
# if smk.wildcards.master_region in ['WAVES_S', 'WAVES_N']:
