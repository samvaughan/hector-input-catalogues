# Running the Code

To run the code re-make the catalogues for all of the Hector regions, all you need to run is:

```bash
snakemake --cores 1 results/RegionCatalogues/{WAVES_N/WAVES_S}/{region_name}/{region_name}_Hector_target_galaxies.csv}
```

from this folder. `{region_name}` refers to G12, G15, etc, or a new region entirely (see [Adding a New Region](docs/adding_a_new_region.md)). Unlike for some of the other pipelines I've written for Hector, there isn't a config file you need to point to. The required input files you'll need for this to work and the steps that the pipeline will carry out are outlined below.

!!! danger
     The default choice of sub-sampling method during the target selection is _randomly selecting galaxies above a given stellar mass/colour_. I've tried very hard to ensure that I'm using a set random seed, such that if you ran this pipeline twice you'd get identical catalogues out. I can't guarantee that this will always be the case if you make changes to the code, however! **Please make a back up of the master catalogues for each region before running this code again**. If you overwrite the master catalogues they might not be recoverable, and this would make a huge mess going forward! This has already happened to me once, which is why I'm raising this point again (it's the reason we have a pipeline step to `add_previously_observed_galaxies_back_to_master_catalogue`- see below).

## Necessary input catalogues

To run this code, you'll need some master input catalogues. These have been created by other people and kindly contributed to the Hector survey. They;re not public, and so aren't stored on the github repo. As of February 2024, these are all stored on the Hector Cloud Storage site at `Hector/Stefania/HectorInputCatalogues.zip`, for which you'll need a Data Central account and to be a member of the Hector team.

In the folder `resources/InitialInputCatalogues` there should be three files: `COMBINED_ApMatch_MSTARV2.5_MattOwers_310123.parquet`, `WAVES_N_Hector.parquet` and `WAVES_S_Hector.parquet`. These contain the catalogues of galaxy targets for the Clusters, WAVES North and WAVES South. 

The cluster input catalogue was made by Matt Owers- talk to him if you have any questions! The WAVES catalogues were created by Sabine Bellstedt from the master WAVES catalogues, with some basic cuts to remove stars. 

You will also need the file `resources/Previous_Redshift_Cats/hemispec.v0.91recommended.parquet`. This is a compendium of literature redshifts put together by Ned Taylor for the 4HS survey. Ask him if you have any questions!

To de-select SAMI galaxies, I'm using Jesse van de Sande's table of kinematic measurements. But really any SAMI catalogue would do here. 

Finally, you'll need a catalogue of standard stars (also known as "F stars" due to their spectral type) for the clusters. Matt Owers has made these, and you should have the files `DeCaLs_DR9_Fstars.fits`, `DeCaLs_DR9_Guide_Stars.fits`, `HECTOR_CRS_fstars_RUN8_DR9.txt` and `HECTOR_CRS_guides_RUN8_DR9.txt` in the `resources/ClusterStarCatalogues/` folder. For the WAVES regions, the required star catalogues will be downloaded when the pipeline is run.

## Other necessary files

There are a few more files which are needed as inputs:

- Information about the redshifts of each Hector cluster, in `resources/ClusterInformation/HectorCluster_Information.txt`.
- The galaxies we've observed up to August 2023 during commissioning. These should be in `resources/Galaxies_Observed_up_to_August_2023/galaxies_observed_up_to_2023.csv`. When the master catalogues were updated in August 2023, some galaxies which had already been observed were randomly not selected for our master catalogues. This catalogue is used to add those back in.
- Observations from the Hector Redshift Survey are stored as `{region_name}_{date}_field{n}_autoz.fits`. These are read in and matched against the WAVES Catalogues during the pipeline.
- The folder `JSON_badclass_votes` contains the results of the people's votes of inspecting the cutout images for galaxies in the Hector cluster regions. If an image is consistently marked as 'bad' for some reason, we remove it from our catalogues.
- In `Misc`, we have two csv files which contain measurements of an F-star spectrum in various different colours. One file refers to the PANSTARRS filters and the other to Skymapper filters. We use this when selecting which stars in our catalogues are likely to be suitable for standards.
- The `RegionInformation` folder contains a single file called `all_regions.csv`. This gives the details of each Region in the Hector Survey (its centre, width, height, etc).

## Pipeline Steps

The pipeline undertakes the following steps:

- `make_redshift_catalogue_from_2dF_obs`: the rule takes our 2dF observations, combined them together and makes a redshift catalogue for use later. The script is `workflow/scripts/make_redshift_catalogue_from_observations.py`.
- `make_badclass_probabilities`: this rule calculates the probability that each of the galaxies in the cluster catalogue is "bad", for some reason (i.e. close to a star, the galaxy is actually an imaging artefact, etc). The script is `workflow/scripts/prepare_badclass_votes_table.py`.
- `prepare_WAVES_catalogue`: this rule takes some initial steps (such as turning flux columns into magnitudes) and makes a modified WAVES North and WAVES South catalogue. The script is `workflow/scripts/prepare_WAVES_catalogues.py`.
- `prepare_cluster_catalogue`: this rule performs similar initial preparatory steps for the cluster catalogue. The script is `workflow/scripts/prepare_cluster_catalogue.py`.
- `run_target_selection`: this rule takes our large master catalogues and applies the Hector selection function to them. We also remove galaxies which small half-light radii, as well as sub-sampling our catalogues to be either "flat in mass" (for the WAVES catalogues) or randomly sub-sampling the cluster red sequence (for the cluster catalogues) to come up with our master target-selected catalogues. The script is `workflow/scripts/run_target_selection.py`.
- `add_previously_observed_galaxies_back_to_master_catalogue`: this rule adds back in galaxies which were observed in 2023/2023 which then ended up being removed when I updated the catalogues and re-ran the random sampling. The script is `workflow/scripts/add_observed_galaxies_to_master_catalogues.py`.
- `separate_into_subregions`: this rule takes our master catalogues and separates them into the individual sub-regions (e.g. WAVES N --> G12, G15, etc). The script is `workflow/scripts/separate_catalogues_into_regions.py`.
- `make_star_catalogues`: this rule downloads star catalogues (standard stars and guide stars) for each region in the Hector survey. The script is `workflow/scripts/select_stars.py`.


There are also equivalent scripts which do the same thing for the cluster catalogue but including low-redshift foreground galaxies in the catalogues too.