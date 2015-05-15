MLTOOLS software
Version 0.126		Tuesday 13 Feb 2007 at 18:00
Copyright (c) 2007 Neil D. Lawrence

This toolbox provides various machine learning tools, either through wrapping other toolboxes (such as NETLAB) or providing the tool directly. It was designed originally as part of splitting the GPLVM and FGPLVM toolboxes.

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

MATLAB Files
------------

Matlab files associated with the toolbox are:

demMppca1.m: Demonstrate MPPCA on a artificial dataset.
imageModify.m: Helper code for visualisation of image data.
imageVisualise.m: Helper code for showing an image during 2-D visualisation.
isomapEmbed.m: Embed data set with Isomap.
kbrCreate.m: Create a KBR model.
kbrDisplay.m: Display parameters of the KBR model.
kbrExpandParam.m: Create model structure from KBR model's parameters.
kbrExtractParam.m: Extract parameters from the KBR model structure.
kbrOptimise.m: Optimise a KBR model.
kbrOptions.m: Create a default options structure for the KBR model.
kbrOut.m: Compute the output of a KBR model given the structure and input X.
kbrOutputGrad.m: Evaluate derivatives of KBR model outputs with respect to parameters.
kbrParamInit.m: KBR model parameter initialisation.
kpcaEmbed.m: Embed data set with kernel PCA.
linearCreate.m: Create a linear model.
linearDisplay.m: Display a linear model.
linearExpandParam.m: Update linear model with vector of parameters.
linearExtractParam.m: Extract weights from a linear model.
linearLogLikeGradients.m: Linear model gradients.
linearLogLikelihood.m: Linear model log likelihood.
linearOptimise.m: Optimise a linear model.
linearOptions.m: Options for learning a linear model.
linearOut.m: Obtain the output of the linear model.
linearOutputGrad.m: Evaluate derivatives of linear model outputs with respect to parameters.
linearOutputGradX.m: Evaluate derivatives of linear model outputs with respect to inputs.
linearParamInit.m: Initialise the parameters of an LINEAR model.
lleEmbed.m: Embed data set with LLE.
lvmClassVisualise.m: Callback function for visualising data in 2-D.
lvmScatterPlot.m: 2-D scatter plot of the latent points.
lvmScatterPlotColor.m: 2-D scatter plot of the latent points with color - for Swiss Roll data.
lvmTwoDPlot.m: Helper function for plotting the labels in 2-D.
mappingOptimise.m: Optimise the given model.
mlpCreate.m: Wrapper for NETLAB's mlp `net'.
mlpDisplay.m: Display the multi-layer perceptron model.
mlpExpandParam.m: Update mlp model with new vector of parameters.
mlpExtractParam.m: Extract weights and biases from an MLP.
mlpLogLikeGradients.m: Multi-layer perceptron gradients.
mlpLogLikeHessian.m: Multi-layer perceptron Hessian.
mlpLogLikelihood.m: Multi-layer perceptron log likelihood.
mlpOptimise.m: Optimise MLP for given inputs and outputs.
mlpOptions.m: Options for the multi-layered perceptron.
mlpOut.m: Output of an MLP model (wrapper for the NETLAB function mlpfwd).
mlpOutputGrad.m: Evaluate derivatives of mlp model outputs with respect to parameters.
mlpOutputGradX.m: Evaluate derivatives of mlp model outputs with respect to inputs.
mlpParamInit.m: Initialise the parameters of an MLP model.
mltoolsToolboxes.m: Load in the relevant toolboxes for the MLTOOLS.
modelCreate.m: Create a model of the specified type.
modelDisplay.m: Display a text output of a model.
modelExpandParam.m: Update a model structure with parameters.
modelExtractParam.m: Extract the parameters of a model.
modelGradient.m: Gradient of error function to minimise for given model.
modelGradientCheck.m: Check gradients of given model.
modelHessian.m: Hessian of error function to minimise for given model.
modelLogLikeGradients.m: Compute a model's gradients wrt log likelihood.
modelLogLikelihood.m: Compute a model log likelihood.
modelObjective.m: Objective function to minimise for given model.
modelOptimise.m: Optimise the given model.
modelOptions.m: Returns a default options structure for the given model.
modelOut.m: Give the output of a model for given X.
modelOutputGrad.m: Compute derivatives with respect to params of model outputs.
modelOutputGradX.m: Compute derivatives with respect to model inputs of model outputs.
modelParamInit.m: Initialise the parameters of the model.
modelPointLogLikelihood.m: Compute the log likelihood of a given point.
modelSamp.m: Give a sample from a model for given X.
modelTest.m: Run some tests on the specified model.
modelTieParam.m: Tie parameters of a model together.
mogCreate.m: Create a mixtures of Gaussians model.
mogEstep.m: Do an E-step on an MOG model.
mogLowerBound.m: Computes lower bound on log likelihood for an MOG model.
mogOptimise.m: Optimise an MOG model.
mogOptions.m: Sets the default options structure for MOG models.
mogUpdateCovariance.m: Update the covariances of an MOG model.
mogUpdateMean.m: Update the means of an MOG model.
mogUpdatePrior.m: Update the priors of an MOG model.
multimodelCreate.m: Create a MULTIMODEL model.
multimodelDisplay.m: Display parameters of the MULTIMODEL model.
multimodelExpandParam.m: Create model structure from MULTIMODEL model's parameters.
multimodelExtractParam.m: Extract parameters from the MULTIMODEL model structure.
multimodelLogLikeGradients.m: Gradient of MULTIMODEL model log likelihood with respect to parameters.
multimodelLogLikelihood.m: Log likelihood of MULTIMODEL model.
multimodelOptions.m: Create a default options structure for the MULTIMODEL model.
multimodelParamInit.m: MULTIMODEL model parameter initialisation.
ppcaEmbed.m: Embed data set with probabilistic PCA.
rbfCreate.m: Wrapper for NETLAB's rbf `net'.
rbfDisplay.m: Display an RBF network.
rbfExpandParam.m: Update rbf model with new vector of parameters.
rbfExtractParam.m: Wrapper for NETLAB's rbfpak.
rbfOptimise.m: Optimise RBF for given inputs and outputs.
rbfOptions.m: Default options for RBF network.
rbfOut.m: Output of an RBF model (wrapper for the NETLAB function rbffwd).
rbfOutputGrad.m: Evaluate derivatives of rbf model outputs with respect to parameters.
rbfperiodicCreate.m: Create a RBFPERIODIC model.
rbfperiodicDisplay.m: Display parameters of the RBFPERIODIC model.
rbfperiodicExpandParam.m: Create model structure from RBFPERIODIC model's parameters.
rbfperiodicExtractParam.m: Extract parameters from the RBFPERIODIC model structure.
rbfperiodicLogLikeGradients.m: Gradient of RBFPERIODIC model log likelihood with respect to parameters.
rbfperiodicLogLikelihood.m: Log likelihood of RBFPERIODIC model.
rbfperiodicOptions.m: Create a default options structure for the RBFPERIODIC model.
rbfperiodicOut.m: Compute the output of a RBFPERIODIC model given the structure and input X.
rbfperiodicOutputGrad.m: Evaluate derivatives of RBFPERIODIC model outputs with respect to parameters.
rbfperiodicOutputGradX.m: Evaluate derivatives of a RBFPERIODIC model's output with respect to inputs.
rbfperiodicParamInit.m: RBFPERIODIC model parameter initialisation.
smallrandEmbed.m: Embed data set with small random values.
spectrumModify.m: Helper code for visualisation of spectrum data.
spectrumVisualise.m: Helper code for showing an spectrum during 2-D visualisation.
vectorModify.m: Helper code for visualisation of vectorial data.
vectorVisualise.m:  Helper code for plotting a vector during 2-D visualisation.
viterbiAlign.m: Compute the Viterbi alignment.
