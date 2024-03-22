# Adding a new Region

As the Hector survey progresses, it will become necessary to create catalogues of targets, standard stars and guide stars for a new region of the sky. This page explains how to do that.

You'll need to make the following changes:

1. Make a backup of all of the current catalogues. You really don't want to overwrite anything!
2. Edit the file `resources/RegionInformation/all_regions.csv` to add the properties of the new region. You'll need to add a name, the minimum and maximum RA and Declination of the region, the area it covers, the "master region" that it comes from (e.g. `WAVES_S`, `WAVES_N` or `HectorClusters`) and the source of the photometry (where we get the guide and standard star catalogues from). Note that the pipeline at the moment can only download catalogues from the SkyMapper or PANSTARRS surveys at the moment.
3. Re-run the pipeline from the main folder: `snakemake --cores 1`. You should see the new region catalogues being made.
4. Check that the master catalogue you've just updated (e.g. the master WAVES N catalogue if you've added a WAVES N region) is the same as it was before, just with extra galaxies added which correspond to the new region of sky. 

!!! danger
    I've not tested adding a new region. There is a real chance of overwriting the master catalogues, which _might_ be unrecoverable due to the random sampling element. I don't _think_ that this is the case (I've added random seeds throughout the code to try and make everything reproducible) but I'm giving you a fair warning here to stress how important it is to **make a backup of the catalogues before you try adding a new region**.