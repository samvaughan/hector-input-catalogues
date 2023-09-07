"""
Take our master cluster star catalogues and separate them up into each region.

Do this using RA and DEC limits set in a separate file (resources/RegionInformation/all_regions.csv). WE also need to rename some columns and add in a 'type' column.
"""
import pandas as pd

smk = snakemake  # noqa
region_information = pd.read_csv(smk.input.region_information, index_col="name")
cluster_region_information = region_information.loc[smk.params.cluster_region_names]

# Get our master star catalogues
master_guide_stars = pd.read_table(
    smk.input.master_guide_star_cat, delim_whitespace=True
)
master_standard_stars = pd.read_table(
    smk.input.master_standard_star_cat, delim_whitespace=True
)

output_standard_star_cats = smk.output.region_standard_star_catalogues
output_guide_star_cats = smk.output.region_guide_star_catalogues

rename_dict = dict(ROWID="ID", PSFMAG_R="r_mag", PRIORITY_FLAG="priority")
master_guide_stars.rename(columns=rename_dict, inplace=True)
master_standard_stars.rename(columns=rename_dict, inplace=True)

master_guide_stars["type"] = 2
master_standard_stars["type"] = 0

master_guide_stars.loc[:, ["pmRA", "pmDEC"]] = 0.0
master_standard_stars.loc[:, ["pmRA", "pmDEC"]] = 0.0

# Add the guide star priority- it's always 8
master_guide_stars["priority"] = 8

# Add letters to the start of each
master_guide_stars["ID"] = master_guide_stars.apply(lambda x: f"G{x.ID}", axis=1)
master_standard_stars["ID"] = master_standard_stars.apply(lambda x: f"S{x.ID}", axis=1)

columns_in_order = ["ID", "RA", "DEC", "r_mag", "type", "pmRA", "pmDEC", "priority"]

for (region, row), standard_filename, guide_filename in zip(
    region_information.iterrows(), output_standard_star_cats, output_guide_star_cats
):
    print(f"Saving stars for {region} to {guide_filename} and {standard_filename}")
    guide_mask = (
        (master_guide_stars.RA > row.min_RA)
        & (master_guide_stars.RA < row.max_RA)
        & (master_guide_stars.DEC > row.min_DEC)
        & (master_guide_stars.DEC < row.max_DEC)
    )
    master_guide_stars.loc[guide_mask, columns_in_order].to_csv(
        guide_filename, index=False
    )

    standard_mask = (
        (master_standard_stars.RA > row.min_RA)
        & (master_standard_stars.RA < row.max_RA)
        & (master_standard_stars.DEC > row.min_DEC)
        & (master_standard_stars.DEC < row.max_DEC)
    )

    master_standard_stars.loc[standard_mask, columns_in_order].to_csv(
        standard_filename, index=False
    )

    print(
        f"\t{region}: Selected {guide_mask.sum()} guide stars and {standard_mask.sum()} standard stars"
    )
