The GP toolbox is an implementation of the GPs that uses the Pseudo-Input method of Snelson and Ghahramani (NIPS 2005) for sparsification as well as extensions given by Quinonero-Candela and Rasmussen (JMLR 2005).


Version 0.13
------------

Changes to allow more flexibility in optimisation of beta.

Version 0.12
------------

Various minor changes for enabling back constraints in hierarchical GP-LVM models.

Version 0.11
------------

Changes include the use of the optimiDefaultConstraint('positive') to obtain the function to constrain beta to be positive (which now returns 'exp' rather than 'negLogLogit' which was previously the default). Similarly default optimiser is now given by a command in optimiDefaultOptimiser.

Version 0.1
-----------

The first version which is spun out of the FGPLVM toolbox. The corresponding FGPLVM toolbox is 0.15. 
