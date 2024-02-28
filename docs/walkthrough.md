# Running the Code

To run the code re-make the catalogues for all of the Hector regions, all you need to run is:

```bash
snakemake --cores 1 results/RegionCatalogues/{WAVES_N/WAVES_S}/{region_name}/{region_name}_Hector_target_galaxies.csv}
```

from this folder. `{region_name}` refers to G12, G15, etc, or a new region entirely (see [Adding a New Region](docs/adding_a_new_region.md)). Unlike for some of the other pipelines I've written for Hector, there isn't a config file you need to point to. The required input files you'll need for this to work and the steps that the pipeline will carry out are outlined below.

**Important note**: the default choice of sub-sampling method during the target selection is _randomly selecting galaxies above a given stellar mass_. I've tried very hard to ensure that I'm using a set random seed, such that if you ran this pipeline twice you'd get identical catalogues out. I can't guarantee that this will always be the case if you make changes to the code, however! **Please make a back up of the master catalogues for each region before running this code again**. I'd *highly* recommend that you only use this code to make catalogues for _new_ regions, and that the catalogues for the clusters, G12, G15, G23, H01 and H03 remain fixed. Otherwise, galaxies which we've already observed might not be included in the re-made catalogues, and this would make a horrible mess going forward!

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