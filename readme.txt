The GP toolbox is an implementation of the GPs that uses the
Pseudo-Input method of Snelson and Ghahramani (NIPS 2005) for
sparsification as well as extensions given by Quinonero-Candela and
Rasmussen (JMLR 2005). Since version 0.132 it also provides an
implementation of Titsias's sparse variational approximation published
in AISTATS 2009.

Version 0.137
-------------

Minor updates to gpLoadResult for allowing different functions for loading in data.

Version 0.136
-------------

Changes to gpReadFromFID for compatibility with C++ release.

Version 0.135
-------------

Modifications by Carl Henrik Ek for compatability with the SGPLVM
toolbox.

Version 0.134
-------------

Updates to allow deconstruction of model files when writing to disk
(gpWriteResult, gpLoadResult, gpDeconstruct, gpReconstruct).

Version 0.133
-------------

Updates for running a GPLVM/GP using the data's inner product matrix
for Interspeech synthesis demos.

Version 0.132
-------------

Remove bug in gpExtractParam for sparse models where scale parameters
were in the wrong place. Moved in some examples from the oxford
toolbox. Added Titsias's variational sparse approximation.

Version 0.131
-------------

Changes to allow compatibility with SGPLVM and NCCA toolboxes.

Version 0.13
------------

Changes to allow more flexibility in optimisation of beta.

Version 0.12
------------

Various minor changes for enabling back constraints in hierarchical
GP-LVM models.

Version 0.11
------------

Changes include the use of the optimiDefaultConstraint('positive') to
obtain the function to constrain beta to be positive (which now
returns 'exp' rather than 'negLogLogit' which was previously the
default). Similarly default optimiser is now given by a command in
optimiDefaultOptimiser.

Version 0.1
-----------

The first version which is spun out of the FGPLVM toolbox. The
corresponding FGPLVM toolbox is 0.15.
