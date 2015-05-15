MLTOOLS software
Version 0.122		Friday 05 Jan 2007 at 23:05
Copyright (c) 2007 Neil D. Lawrence

This toolbox provides various machine learning tools, either through wrapping other toolboxes (such as NETLAB) or providing the tool directly. It was designed originally as part of splitting the GPLVM and FGPLVM toolboxes.


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

lvmScatterPlot.m: 2-D scatter plot of the latent points.
mlpLogLikelihood.m: Multi-layer perceptron log likelihood.
modelExpandParam.m: Update a model structure with parameters.
mlpOut.m: Output of an MLP model (wrapper for the NETLAB function mlpfwd).
ppcaEmbed.m: Embed data set with probabilistic PCA.
mogUpdateCovariance.m: Update the covariances of an MOG model.
kbrOptimise.m: Optimise a KBR model.
mogOptimise.m: Optimise an MOG model.
modelParamInit.m: Initialise the parameters of the model.
rbfDisplay.m: Display an RBF network.
linearOutputGradX.m: Evaluate derivatives of linear model outputs with respect to inputs.
viterbiAlign.m: Compute the Viterbi alignment.
mogEstep.m: Do an E-step on an MOG model.
rbfOutputGrad.m: Evaluate derivatives of rbf model outputs with respect to parameters.
kbrOutputGrad.m: Evaluate derivatives of KBR model outputs with respect to parameters.
modelExtractParam.m: Extract the parameters of a model.
linearOut.m: Obtain the output of the linear model.
mlpExpandParam.m: Update mlp model with new vector of parameters.
mlpOptimise.m: Optimise MLP for given inputs and outputs.
isomapEmbed.m: Embed data set with Isomap.
mlpLogLikeGradients.m: Multi-layer perceptron gradients.
mogCreate.m: Create a mixtures of Gaussians model.
linearLogLikeGradients.m: Linear model gradients.
linearOptions.m: Options for learning a linear model.
spectrumVisualise.m: Helper code for showing an spectrum during 2-D visualisation.
modelHessian.m: Hessian of error function to minimise for given model.
modelOutputGrad.m: Compute derivatives with respect to params of model outputs.
modelGradient.m: Gradient of error function to minimise for given model.
mlpOutputGrad.m: Evaluate derivatives of mlp model outputs with respect to parameters.
imageModify.m: Helper code for visualisation of image data.
demMppca1.m: Demonstrate MPPCA on a artificial dataset.
vectorVisualise.m:  Helper code for plotting a vector during 2-D visualisation.
modelDisplay.m: Display a text output of a model.
kbrOut.m: Compute the output of a KBR model given the structure and input X.
rbfExtractParam.m: Wrapper for NETLAB's rbfpak.
modelOptimise.m: Optimise the given model.
lvmClassVisualise.m: Callback function for visualising data in 2-D.
kbrOptions.m: Create a default options structure for the KBR model.
mogOptions.m: Sets the default options structure for MOG models.
mogUpdateMean.m: Update the means of an MOG model.
modelLogLikelihood.m: Compute a model log likelihood.
kpcaEmbed.m: Embed data set with kernel PCA.
rbfExpandParam.m: Update rbf model with new vector of parameters.
rbfOut.m: Output of an RBF model (wrapper for the NETLAB function rbffwd).
linearExpandParam.m: Update linear model with vector of parameters.
kbrExpandParam.m: Create model structure from KBR model's parameters.
modelOutputGradX.m: Compute derivatives with respect to model inputs of model outputs.
modelObjective.m: Objective function to minimise for given model.
mlpOptions.m: Options for the multi-layered perceptron.
mlpParamInit.m: Initialise the parameters of an MLP model.
rbfOptions.m: Default options for RBF network.
mlpLogLikeHessian.m: Multi-layer perceptron Hessian.
vectorModify.m: Helper code for visualisation of vectorial data.
mogLowerBound.m: Computes lower bound on log likelihood for an MOG model.
linearDisplay.m: Display a linear model.
mlpOutputGradX.m: Evaluate derivatives of mlp model outputs with respect to inputs.
mlpCreate.m: Wrapper for NETLAB's mlp `net'.
mappingOptimise.m: Optimise the given model.
modelCreate.m: Create a model of the specified type.
modelLogLikeGradients.m: Compute a model's gradients wrt log likelihood.
mogUpdatePrior.m: Update the priors of an MOG model.
linearCreate.m: Create a linear model.
lvmScatterPlotColor.m: 2-D scatter plot of the latent points with color - for Swiss Roll data.
modelOut.m: Give the output of a model for given X.
mlpExtractParam.m: Extract weights and biases from an MLP.
lleEmbed.m: Embed data set with LLE.
modelTieParam.m: Tie parameters of a model together.
modelPointLogLikelihood.m: Compute the log likelihood of a given point.
kbrParamInit.m: KBR model parameter initialisation.
linearOptimise.m: Optimise a linear model.
spectrumModify.m: Helper code for visualisation of spectrum data.
kbrCreate.m: Create a KBR model.
kbrDisplay.m: Display parameters of the KBR model.
lvmTwoDPlot.m: Helper function for plotting the labels in 2-D.
linearExtractParam.m: Extract weights from a linear model.
modelSamp.m: Give a sample from a model for given X.
linearOutputGrad.m: Evaluate derivatives of linear model outputs with respect to parameters.
linearParamInit.m: Initialise the parameters of an LINEAR model.
imageVisualise.m: Helper code for showing an image during 2-D visualisation.
modelTest.m: Run some tests on the specified model.
linearLogLikelihood.m: Linear model log likelihood.
modelOptions.m: Returns a default options structure for the given model.
rbfOptimise.m: Optimise RBF for given inputs and outputs.
kbrExtractParam.m: Extract parameters from the KBR model structure.
rbfCreate.m: Wrapper for NETLAB's rbf `net'.
mlpDisplay.m: Display the multi-layer perceptron model.
smallrandEmbed.m: Embed data set with small random values.
