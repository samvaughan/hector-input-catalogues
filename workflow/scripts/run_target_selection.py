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

# import scipy.stats as stats

smk = snakemake

config_filename = smk.input.target_selection_config_file

HP = pipeline.HectorPipe(config_filename=config_filename)

fig, axs = HP.run_target_selection(save=False, plot=True)
fig.savefig(smk.output.target_selection_plot, bbox_inches="tight")
plt.close("all")

# If we're dealing with the clusters, add back in the individual galaxy redshifts
if smk.wildcards['master_region'] == 'HectorClusters':
    # Read the master input catalogue
    master_input_catalogue = pd.read_parquet(smk.input.input_catalogue)
    
    # Make sure our indices are all integer types
    master_input_catalogue.index = pd.to_numeric(master_input_catalogue.index)
    HP.target_selection.selection_function_sparsely_sampled.index = pd.to_numeric(HP.target_selection.selection_function_sparsely_sampled.index)
    HP.target_selection.selection_function.index = pd.to_numeric(HP.target_selection.selection_function.index)
    
    # Now just set the ['z'] columns to be equal- pandas takes care of the matching
    HP.target_selection.selection_function['z'] = master_input_catalogue['z']
    HP.target_selection.selection_function_sparsely_sampled['z'] = master_input_catalogue['z']


# Make sure to save the *sparsely sampled* catalogues!!
HP.target_selection.selection_function_sparsely_sampled.to_parquet(
    smk.output.master_target_selected_catalogue
)

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
