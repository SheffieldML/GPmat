
Matlab Datasets Toolbox
=======================


Release Information
-------------------

**Current release is 0.1371**.

As well as downloading the DATASETS software you need to obtain the toolboxes specified below. 

- [ndlutil](https://github.com/SheffieldML/ndlutil/) some general untilities.

Minor mods to lvmLoadData for Spellman data.

#### Version 0.137

Added loaders for speech grid corpus data.

#### Version 0.136

Added loader for Agarwal and Triggs silhouette data back in, slight modifications of robot wireless datasets.

#### Version 0.135

Added loader for Swiss Jura data and toy data for multiple output GPs.

#### Version 0.134

Added loader for the movielens data.

#### Version 0.133

Added some variants of the swiss roll "data".

#### Version 0.132

Added loader for Spellman data and the crabs data from Ripley's book. You will need the Spellman data from [here](http://genome-www.stanford.edu/cellcycle/data/rawdata/combined.txt) and the crabs data from [here](http://www.stats.ox.ac.uk/pub/PRNN/crabs.dat).

#### Version 0.131

Moved functionality from mappingLoadData, ivmLoadData and ncnmLoadData into mapLoadData.

#### Version 0.13

Added processing for CMU data base motion capture figure 35. To load, the data should be downloaded from the CMU motion capture data base, http://mocap.cs.cmu.edu, and placed in baseDir/mocap/cmu/35/, where baseDir is the directory where you have placed 'datasetsDirectory.m'. You need the amc files and the asf file.

Removed Example 86 motion 10 from the CMU motion capture data base.

#### Version 0.121

Modified robot traces example to reflect the fact that the data which is missing is genuinely missing rather than just reading low.

#### Release 0.12

There seems to have been a problem in earlier releases in terms of where the datasets were placed. This has now (hopefully) been fixed using the datasetsDirectory command which returns the directory in which that command is placed. All the data should be placed in that directory.

#### Release 0.11

A new data set has been added, motion capture in x,y,z, format of a run. This [data](http://accad.osu.edu/research/mocap/mocap_data.htm) is from Ohio state's Advanced Computing Center fro the Arts and Design web-site.

The mapLoadData.m command has been added to load in data for regression and classification problems.

#### Release 0.1

First release of the toolbox.

The toolbox currently contains two different files for loading datasets into matlab and a collection of matlab data sets.

Page updated on Wed May 19 09:36:44 2010


