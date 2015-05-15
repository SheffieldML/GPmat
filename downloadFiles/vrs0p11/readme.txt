MLTOOLS software
Version 0.11		Friday 20 Jan 2006 at 18:45
Copyright (c) 2006 Neil D. Lawrence

This toolbox provides various machine learning tools, either through wrapping other toolboxes (such as NETLAB) or providing the tool directly. It was designed originally as part of splitting the GPLVM and FGPLVM toolboxes.


Version 0.11
------------

Updated release for operation with FGPLVM toolbox 0.13. Structure of model creation changed and functions of the form modelOptions.m included for setting default options of the various machine learning models.

Version 0.1
-----------

The first release of the toolbox with various wrappers for NETLAB functions. Also latent variable model visualisation code was moved into this toolbox.

MATLAB Files
------------

Matlab files associated with the toolbox are:

linearCreate.m: Create a linear model.
isomapEmbed.m: Embed data set with Isomap.
kpcaEmbed.m: Embed data set with kernel PCA.
lvmClassVisualise.m: Callback function for visualising data in 2-D.
imageVisualise.m: Helper code for showing an image during 2-D visualisation.
kbrExpandParam.m: Update kernel based regression model with vector of parameters.
kbrExtractParam.m: Extract weights from a kernel based regression model.
kbrOptimise.m: Optimise a kernel based regression.
kbrOut.m: Obtain the output of the kernel based regression model.
kbrOutputGrad.m: Evaluate derivatives of kernel based regression model outputs with respect to parameters.
linearExpandParam.m: Update linear model with vector of parameters.
linearExtractParam.m: Extract weights from a linear model.
linearOptimise.m: Optimise a linear model.
linearOut.m: Obtain the output of the linear model.
linearOutputGrad.m: Evaluate derivatives of linear model outputs with respect to parameters.
mappingOptimise.m: Optimise the given model.
mlpCreate.m: Wrapper for NETLAB's mlp `net'.
mlpExpandParam.m: Update mlp model with new vector of parameters.
mlpOut.m: Output of an MLP model (wrapper for the NETLAB function mlpfwd).
mlpExtractParam.m: Wrapper for NETLAB's mlppak.
mlpOptimise.m: Optimise MLP for given inputs and outputs.
mlpOutputGrad.m: Evaluate derivatives of mlp model outputs with respect to parameters.
modelExpandParam.m: Update a model structure with parameters.
modelExtractParam.m: Extract the parameters of a model.
modelOptimise.m: Optimise the given model.
modelOut.m: Give the output of a model for given X.
modelOutputGrad.m: Compute derivatives with respect to params of model outputs.
rbfCreate.m: Wrapper for NETLAB's rbf `net'.
rbfExpandParam.m: Update rbf model with new vector of parameters.
rbfExtractParam.m: Wrapper for NETLAB's rbfpak.
rbfOptimise.m: Optimise RBF for given inputs and outputs.
rbfOut.m: Output of an RBF model (wrapper for the NETLAB function rbffwd).
rbfOutputGrad.m: Evaluate derivatives of rbf model outputs with respect to parameters.
lvmScatterPlot.m: 2-D scatter plot of the latent points.
ppcaEmbed.m: Embed data set with probabilistic PCA.
spectrumModify.m: Helper code for visualisation of spectrum data.
spectrumVisualise.m: Helper code for showing an spectrum during 2-D visualisation.
kbrCreate.m: Create a kernel based regression model.
vectorModify.m: Helper code for visualisation of vectorial data.
vectorVisualise.m:  Helper code for plotting a vector during 2-D visualisation.
lvmTwoDPlot.m: Helper function for plotting the labels in 2-D.
modelLogLikeGradients.m: Compute a model's gradients wrt log likelihood.
modelLogLikelihood.m: Compute a model log likelihood.
lvmScatterPlotColor.m: 2-D scatter plot of the latent points with color - for Swiss Roll data.
modelSamp.m: Give a sample from a model for given X.
modelDisplay.m: Display a text output of a model.
mlpDisplay.m: Display the multi-layer perceptron model.
mlpOptions.m: Options for the multi-layered perceptron.
kbrOptions.m: Kernel based regression options.
rbfDisplay.m: Display an RBF network.
linearOptions.m: Options for learning a linear model.
modelCreate.m: Create a model of the specified type.
rbfOptions.m: Default options for RBF network.
