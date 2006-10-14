This toolbox implements some generic functions used by other toolboxes. Originally it was spun out of the IVM 0.221 toolbox.


Version 0.152
-------------

There was a sign error in lnDiffCumGaussian. This has been fixed. The sign error was replicated in the NOISE and PRIOR toolboxes (a case of two wrongs making a 'right'). This means that earlier versions of this toolbox are incompatible with versions of PRIOR before 0.131 and versions of NOISE before 0.131.

Version 0.151
-------------

Sped up the sparseDiag command using spdiag, this makes the FGPLVM code run much faster.

Version 0.15
------------

Fixed inconsistencies in logdet and pdinv and corrected the amount of jitter that is added at first try. Placed the jitter addition in a new file jitChol which is called by both functions. 

Added kldivGaussian for computing Gaussian Kullback-Leibler divergences.

Added deg2rad for converting degrees to radians.

Version 0.142
-------------

Added sparseDiag and moved data loading methods into the new DATASETS toolbox. 

Version 0.14
------------

Moved in lvmLoadData and mappingLoadData as generic methods of loading in data sets. 

Version 0.131 Release Notes
---------------------------

Added files getline.m and tokenise.m which are string utilities for reading and processing files.

Version 0.13 Release Notes
--------------------------

pdinv now adds jitter which is proportional to the mean of the diagonal elements of the matrix.


