The FGPLVM toolbox is a faster implementation of the GP-LVM. It uses
sparse variants of Gaussian processes that include the Pseudo-Input
method of Snelson and Ghahramani (NIPS 2005) and the sparse
variational approach of Titsias (AISTATS 2009) for sparsification and
efficiency improvements.

Version 0.163
-------------

Changes for compatibility with new SGPLVM toolbox by Carl Henrik Ek.

Version 0.162
-------------

Added new files fgplvmWriteResults fgplvmLoadResults for saving
smaller model files.

Version 0.161
-------------

Updates for running a GPLVM when the inner produce matrix is used
(i.e. dimensionality much greater than data points).

Version 0.16
------------

Incorporate Michalis's variational sparse approximation in the toolbox.
Minor changes to fix reading of GPLVM files from latest C++ code.

Version 0.153
-------------

Changes to allow compatibility with SGPLVM and NCCA toolboxes.


Version 0.152
-------------

Bug fix from fgplvmReadFromFID where the values of model.m weren't
being computed correctly.

Version 0.151
-------------

In this version results for the CMU Mocap data set from Taylor et
al. of subject 35 running and walking are included, as well as some
minor changes to allow hierarchical GP-LVMs to be used.

Version 0.15
------------

This version splits the Gaussian process portion into a new GP
toolbox, the corresponding version is 0.1. Fixed bug in
gpDynamicsExpandParam, gpDynamicsExractParam and
gpDynamicsLogLikeGradient where 'fixInducing' option was not being
dealt with.

Fixed bug in fgplvmCreate.m where the back constraints were set up,
but the latent positions were not being set according to the back
constraints in the returned model.

Version 0.141
-------------

Changed GP-LVM default optimiser to scg rather than conjgrad. Added
fgplvmOptimiseSequence and dependent files. This is for optimising a
test sequence in the latent space, for the case where there are
dynamics on the model.

Version 0.14
------------

Carl Henrik Ek implemented multiple sequences in the gpDynamics model
used for dynamics in the GPLVM, this was refined and integrated by
Neil.

Fixed two bugs in gpPosteriorGradMeanVar which appeared if fitc was
used or the scales on the outputs were non-zero. This in turn affected
fgplvmOptimisePoint.

Default under back constraints switched to not optimise towards a PCA
initialisation.

Fixed bug in fgplvmReadFromFID where the old form of fgplvmCreate was
being called.

Version 0.132
-------------

Learning with missing data fully implemented across all models. Two
big speed improvements on the fitc approximation (thanks to Ed Snelson
for pointing out how slow it was!).

Version 0.131
-------------

Added learning with missing data for the FTC and reversible dynamics
model.

Version 0.13
------------

This version includes a much cleaner way of incorporating different
dynamics models. It is released in line two imminent reports on
learning large scale Gaussian processes and learning with back
constraints.

Version 0.11
------------

This version now includes the Snelson-Ghahramani approximation (called
FITC by Quinonero-Candela and Rasmussen) and the partially independent
training criterion (PITC). Additionally the approximations can be used
in standard Gaussian process regression.

Version 0.1
-----------

In the first release, only the projected latent variables
approximation of Seeger et al is implemented. The toolbox also
implements back-constraints as proposed by Lawrence and
Quinonero-Candela

The first release containing a couple of examples on the oil data
(demOil1.m and demOil2.m) and dynamics (demStick1.m and
demStick2.m). The toolbox can load in the C++ files with dynamics
associated and (through the mocapResultsCppBvh in the MOCAP toolbox)
can run motion capture files with dynamics.
