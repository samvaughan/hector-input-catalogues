"""
Take our master region catalogues and separate them up into each region.

Do this using RA and DEC limits set in a separate file
"""
import pandas as pd

smk = snakemake
region_information = pd.read_csv(smk.input.region_info_file, index_col='name')
master_df = pd.read_csv(smk.input.master_catalogue)

out_names = smk.output.region_target_catalogues

for (region, row), filename in zip(region_information.iterrows(), out_names):

    # if smk.params.component != 'Clusters':
    region_mask = (master_df.RA > row.min_RA) & (master_df.RA < row.max_RA) & (master_df.DEC > row.min_DEC) & (master_df.DEC < row.max_DEC)
    master_df.loc[region_mask].to_csv(filename)
    print(f"Selected {region_mask.sum()} galaxies for the {region} region")

    # else:
    #     raise NotImplementedError('Need to add the code for the clusters here')
