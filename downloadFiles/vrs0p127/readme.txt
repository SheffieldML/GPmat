MLTOOLS software
Version 0.127		Tuesday 22 May 2007 at 22:57
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

multimodelOptions.m: Create a default options structure for the MULTIMODEL model.
lleEmbed.m: Embed data set with LLE.
rbfperiodicOutputGrad.m: Evaluate derivatives of RBFPERIODIC model outputs with respect to parameters.
rbfOptimise.m: Optimise RBF for given inputs and outputs.
modelAddDynamics.m: Add a dynamics kernel to the model.
modelOutputGradX.m: Compute derivatives with respect to model inputs of model outputs.
rbfperiodicExtractParam.m: Extract parameters from the RBFPERIODIC model structure.
rbfperiodicOutputGradX.m: Evaluate derivatives of a RBFPERIODIC model's output with respect to inputs.
multimodelExtractParam.m: Extract parameters from the MULTIMODEL model structure.
rbfperiodicExpandParam.m: Create model structure from RBFPERIODIC model's parameters.
multimodelExpandParam.m: Create model structure from MULTIMODEL model's parameters.
rbfperiodicLogLikelihood.m: Log likelihood of RBFPERIODIC model.
mogEstep.m: Do an E-step on an MOG model.
kbrExpandParam.m: Create model structure from KBR model's parameters.
linearOutputGradX.m: Evaluate derivatives of linear model outputs with respect to inputs.
rbfExpandParam.m: Update rbf model with new vector of parameters.
kbrOptimise.m: Optimise a KBR model.
spectrumVisualise.m: Helper code for showing an spectrum during 2-D visualisation.
mlpExtractParam.m: Extract weights and biases from an MLP.
imageModify.m: Helper code for visualisation of image data.
linearParamInit.m: Initialise the parameters of an LINEAR model.
modelExtractParam.m: Extract the parameters of a model.
mlpOut.m: Output of an MLP model.
multimodelLogLikelihood.m: Log likelihood of MULTIMODEL model.
modelExpandParam.m: Update a model structure with parameters.
ppcaEmbed.m: Embed data set with probabilistic PCA.
kbrDisplay.m: Display parameters of the KBR model.
mappingOptimise.m: Optimise the given model.
mlpOutputGrad.m: Evaluate derivatives of mlp model outputs with respect to parameters.
mogOptions.m: Sets the default options structure for MOG models.
viterbiAlign.m: Compute the Viterbi alignment.
rbfOut.m: Output of an RBF model (wrapper for the NETLAB function rbffwd).
kbrOut.m: Compute the output of a KBR model given the structure and input X.
mlpExpandParam.m: Update mlp model with new vector of parameters.
modelGradient.m: Gradient of error function to minimise for given model.
rbfOutputGrad.m: Evaluate derivatives of rbf model outputs with respect to parameters.
kbrOutputGrad.m: Evaluate derivatives of KBR model outputs with respect to parameters.
modelOptimise.m: Optimise the given model.
multimodelDisplay.m: Display parameters of the MULTIMODEL model.
modelTieParam.m: Tie parameters of a model together.
linearOut.m: Obtain the output of the linear model.
mlpParamInit.m: Initialise the parameters of an MLP model.
rbfperiodicOut.m: Compute the output of a RBFPERIODIC model given the structure and input X.
linearOutputGrad.m: Evaluate derivatives of linear model outputs with respect to parameters.
linearExtractParam.m: Extract weights from a linear model.
imageVisualise.m: Helper code for showing an image during 2-D visualisation.
vectorVisualise.m:  Helper code for plotting a vector during 2-D visualisation.
modelParamInit.m: Initialise the parameters of the model.
mlpLogLikeGradients.m: Multi-layer perceptron gradients.
mlpOptions.m: Options for the multi-layered perceptron.
modelOut.m: Give the output of a model for given X.
mogOptimise.m: Optimise an MOG model.
mlpCreate.m: Multi-layer peceptron model.
linearCreate.m: Create a linear model.
mogLowerBound.m: Computes lower bound on log likelihood for an MOG model.
modelOutputGrad.m: Compute derivatives with respect to params of model outputs.
smallrandEmbed.m: Embed data set with small random values.
rbfperiodicDisplay.m: Display parameters of the RBFPERIODIC model.
modelCreate.m: Create a model of the specified type.
linearOptimise.m: Optimise a linear model.
multimodelCreate.m: Create a MULTIMODEL model.
kbrOptions.m: Create a default options structure for the KBR model.
rbfperiodicOptions.m: Create a default options structure for the RBFPERIODIC model.
mlpDisplay.m: Display the multi-layer perceptron model.
lvmClassVisualise.m: Callback function for visualising data in 2-D.
vectorModify.m: Helper code for visualisation of vectorial data.
rbfOptions.m: Default options for RBF network.
modelLogLikeGradients.m: Compute a model's gradients wrt log likelihood.
modelLogLikelihood.m: Compute a model log likelihood.
mogUpdatePrior.m: Update the priors of an MOG model.
linearLogLikelihood.m: Linear model log likelihood.
modelDisplay.m: Display a text output of a model.
mlpOptimise.m: Optimise MLP for given inputs and outputs.
mlpLogLikelihood.m: Multi-layer perceptron log likelihood.
mlpOutputGradX.m: Evaluate derivatives of mlp model outputs with respect to inputs.
mltoolsToolboxes.m: Load in the relevant toolboxes for the MLTOOLS.
linearOptions.m: Options for learning a linear model.
rbfperiodicCreate.m: Create a RBFPERIODIC model.
modelGradientCheck.m: Check gradients of given model.
lvmScatterPlotColor.m: 2-D scatter plot of the latent points with color - for Swiss Roll data.
multimodelParamInit.m: MULTIMODEL model parameter initialisation.
lvmTwoDPlot.m: Helper function for plotting the labels in 2-D.
linearExpandParam.m: Update linear model with vector of parameters.
isomapEmbed.m: Embed data set with Isomap.
mogCreate.m: Create a mixtures of Gaussians model.
rbfExtractParam.m: Wrapper for NETLAB's rbfpak.
mlpLogLikeHessian.m: Multi-layer perceptron Hessian.
rbfperiodicParamInit.m: RBFPERIODIC model parameter initialisation.
modelSamp.m: Give a sample from a model for given X.
rbfCreate.m: Wrapper for NETLAB's rbf `net'.
mogUpdateMean.m: Update the means of an MOG model.
modelOptions.m: Returns a default options structure for the given model.
modelTest.m: Run some tests on the specified model.
kbrParamInit.m: KBR model parameter initialisation.
lvmScatterPlot.m: 2-D scatter plot of the latent points.
multimodelLogLikeGradients.m: Gradient of MULTIMODEL model log likelihood with respect to parameters.
modelObjective.m: Objective function to minimise for given model.
spectrumModify.m: Helper code for visualisation of spectrum data.
kpcaEmbed.m: Embed data set with kernel PCA.
demMppca1.m: Demonstrate MPPCA on a artificial dataset.
modelPointLogLikelihood.m: Compute the log likelihood of a given point.
rbfDisplay.m: Display an RBF network.
kbrCreate.m: Create a KBR model.
rbfperiodicLogLikeGradients.m: Gradient of RBFPERIODIC model log likelihood with respect to parameters.
linearLogLikeGradients.m: Linear model gradients.
modelHessian.m: Hessian of error function to minimise for given model.
kbrExtractParam.m: Extract parameters from the KBR model structure.
linearDisplay.m: Display a linear model.
mogUpdateCovariance.m: Update the covariances of an MOG model.
