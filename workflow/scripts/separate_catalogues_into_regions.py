"""
Take our master region catalogues and separate them up into each region.

Do this using RA and DEC limits set in a separate file (resources/RegionInformation/all_regions.csv)
"""
import pandas as pd

smk = snakemake
region_information = pd.read_csv(smk.input.region_information, index_col="name")

# Load in all the catalogues we have
master_df = pd.DataFrame()
for catalogue in smk.input.master_catalogues:
    new_cat = pd.read_parquet(catalogue)
    master_df = pd.concat((master_df, new_cat))


output_filenames = smk.output.all_region_catalogues

for (region, row), filename in zip(region_information.iterrows(), output_filenames):

    print(f"Saving {region} to {filename}:")
    region_mask = (
        (master_df.RA > row.min_RA)
        & (master_df.RA < row.max_RA)
        & (master_df.DEC > row.min_DEC)
        & (master_df.DEC < row.max_DEC)
    )
    master_df.loc[region_mask].to_csv(filename, index=False)
    print(f"\tSelected {region_mask.sum()} galaxies for {region}")
