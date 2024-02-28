# The Hector Input Catalogues

This site describes the code I've written to go from the initial catalogues containing millions of objects to the final catalogues of targets that will be part of the Hector survey. 

The Hector Survey is made up of two large patches of sky, named WAVES North and WAVES South, as well as 11 dense galaxy clusters[^1] at a range of Right Ascensions. To begin with, we've sub-divided the large WAVES North and South regions of sky into five smaller regions of roughly 60 square degrees: G12, G15, G23, H01 and H03. As the survey progresses, we'll add more of these.

[^1]: Really more like 10- the clusters A3391 and A3395 are merging, and are treated as a single field in the pipeline.

Along with the clusters, this makes 15 unique regions which this pipeline will create catalogues for. The ultimate outputs of this pipeline are the catalogues of galaxy targets, guide stars and standard stars for each of these regions.

The list of Hector regions and the sky coordinates of their centres (in degrees) is below, and can also be downloaded as a .csv file [here](files/HectorTargets.csv){:download="HectorRegions.csv"}

| Region Name | Right Ascension | Declination |
| ------ | ----------- | --------- |
| A0151  |   17.10920  |  -15.40920|
| A3158  |   55.77040  |  -53.65310|
| A3266  |   67.77460  |  -61.44360|
| A3376  |   90.15290  |  -40.03260|
| A3391  |   96.58590  |  -53.69330|
| A3395  |   96.88000  |  -54.43740|
| A3667  |  303.09170  |  -56.81520|
| A3716  |  312.86000  |  -52.70700|
| A2399  |  329.372605 |  -7.796920|
| A0119  |  14.067150  |  -1.255370|
| A0085  |  10.460211  |  -9.303184|
| G12    |  180        |  0.0      |
| G15    |  225        |  0.0      |
| G23    |  345        |  -32.5    |
| H01    |  15         |  30       |
| H03    |  45         |  30       |


