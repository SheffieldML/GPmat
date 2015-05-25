This toolbox provides various machine learning tools, either through wrapping other toolboxes (such as NETLAB) or providing the tool directly. It was designed originally as part of splitting the GPLVM and FGPLVM toolboxes.

Version 0.139
-------------

Fix to re-enable HGPLVM visualization.

Version 0.138
-------------

Minor tweak of model write result and model load result to allow specification of the data loading function.

Version 0.137
-------------

Release for latest version of VARGPLVM toolbox.

Version 0.136
-------------

Minor Mods.

Version 0.135
-------------

Minor mods.


Version 0.134
-------------

Added pmvu model.

Version 0.133
-------------

Added functionality for writing model files using modelDeconstruct commands to keep written files smaller.

Version 0.132
-------------

Add click visualise functionality for LVM visualization, Laplacian eigenmaps and wrapper for MVU. 

Version 0.1311
--------------

Minor change to lvmScatterPlot to fix bug caused when minimum values were positive.

Version 0.131
-------------

Minor changes to toolbox to fix reading in of files written by C++ code.

Version 0.13
------------

Added paramNameRegularExpressionLookup.m to regular expression match a parameter name in a model and return the associated indices. paramNameReverseLookup.m does the same thing but for the specific parameter name. Also added multimodel type, which allows for multi-task style learning of existing models. Added linear mapping type of model. 

Version 0.1291
--------------

Changes to modelOutputGrad.m, modelOut.m, kbrOutputGrad.m, kbrExpandParam.m, modelOptimise.m to allow compatibility with SGPLVM and NCCA toolboxes. Added a preliminary coding of LLE.


Version 0.129
-------------

Added dnet type model for GTM and density networks. Added various lvm helper files for doing nearest neighbour and plotting results for latent variable models. Added lmvu and mvu embedding wrapper. Added ppca model type. Added output gradients for model out functions (for magnification factor computation in dnet models). Added helpers for reading various models from FID mapmodel, matrix etc.).
Added rbfOutputGradX and visualisation for spring dampers type.

Version 0.128
-------------

Fixed Viterbi alignment algorithm, thanks to Raquel Urtasun for pointing out the problems with it.

Carl Henrik Ek added embeddings with maximum variance unfolding (landmark and normal) to the toolbox. Also several files modified by Carl to allow a single output dimension of a model to be manipulated.

Version 0.127
-------------

Minor modifications including adding file modelAddDynamics to replace fgplvmAddDynamics.

Version 0.126
-------------

Modified kbrParamInit to scale alpha weights and biases by number of data. Added 'dynamicsSliderChange' to lvmClassVisualise to allow visulisation of models with 'gpTime' style-dynamics.

Version 0.125
-------------

Added multimodel for learning multiple indepedent models with shared parameters.

Version 0.124
-------------

Added model gradient checker and added rbfperiodic function to provide a length scale for the gibbsperiodic kernel.

Version 0.123
-------------

Minor release in line with IVM toolbox 0.4.

Version 0.122
-------------

Added Hessian code for base model type and for MLP. Added Viterbi alignment code, viterbiAlign.

Version 0.121
-------------

Various minor bug fixes and changes which seem to have gone undocumented.

Version 0.12
------------

Extended model type to be a generic container module for optimising any model. Added model test for testing a created model. The code is still in a bit of flux though with some design decisions not made and some code untested.

Version 0.111
-------------

Fixed bug in kbr where bias parameter fields where still being referred to as b.Also acknowledged the fact that the isomap graph may not be fully connected in isomapEmbed, but don't yet deal with it properly. Finally added lleEmbed.m for wrapping the lle code.


Version 0.11
------------

Updated release for operation with FGPLVM toolbox 0.13. Structure of model creation changed and functions of the form modelOptions.m included for setting default options of the various machine learning models.

Version 0.1
-----------

The first release of the toolbox with various wrappers for NETLAB functions. Also latent variable model visualisation code was moved into this toolbox.
