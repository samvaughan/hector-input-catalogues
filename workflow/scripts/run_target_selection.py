"""
Run the target selection on our input catalogues.

* Run on the entire, combined master catalogues for WAVES South, WAVES North and the clusters. We then break up into smaller regions in the next script.
"""

from hop import pipeline
import matplotlib.pyplot as plt
import shutil
from pathlib import Path
# import scipy.stats as stats

smk = snakemake

config_filename = smk.input.target_selection_config_file

HP = pipeline.HectorPipe(config_filename=config_filename)

fig, axs = HP.run_target_selection(save=False, plot=True)
fig.savefig(smk.output.target_selection_plot, bbox_inches="tight")
plt.close("all")

HP.target_selection.selection_function.to_parquet(
    smk.output.master_target_selected_catalogue
)

# Clean up our output folder
base_output_folder = Path(smk.output.master_target_selected_catalogue).parent
folders_to_delete = ["Allocation", "Configuration", "DistortionCorrected", "Fibres", "FinalOutputs", "Logs", "Plots", "Tiles"]
for folder in folders_to_delete:
    path = base_output_folder / folder
    shutil.rmtree(path)

# # Run some stats on the WAVES catalogues
# if smk.wildcards.master_region in ['WAVES_S', 'WAVES_N']:
