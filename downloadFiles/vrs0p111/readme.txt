MLTOOLS software
Version 0.111		Saturday 18 Feb 2006 at 21:06
Copyright (c) 2006 Neil D. Lawrence

This toolbox provides various machine learning tools, either through wrapping other toolboxes (such as NETLAB) or providing the tool directly. It was designed originally as part of splitting the GPLVM and FGPLVM toolboxes.


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

imageVisualise.m: Helper code for showing an image during 2-D visualisation.
isomapEmbed.m: Embed data set with Isomap.
kbrCreate.m: Create a kernel based regression model.
kbrExpandParam.m: Update kernel based regression model with vector of parameters.
kbrExtractParam.m: Extract weights from a kernel based regression model.
kbrOptimise.m: Optimise a kernel based regression.
kbrOptions.m: Kernel based regression options.
kbrOut.m: Obtain the output of the kernel based regression model.
kbrOutputGrad.m: Evaluate derivatives of kernel based regression model outputs with respect to parameters.
kpcaEmbed.m: Embed data set with kernel PCA.
linearCreate.m: Create a linear model.
linearExpandParam.m: Update linear model with vector of parameters.
linearExtractParam.m: Extract weights from a linear model.
linearOptimise.m: Optimise a linear model.
linearOptions.m: Options for learning a linear model.
linearOut.m: Obtain the output of the linear model.
linearOutputGrad.m: Evaluate derivatives of linear model outputs with respect to parameters.
lleEmbed.m: Embed data set with LLE.
lvmClassVisualise.m: Callback function for visualising data in 2-D.
lvmScatterPlot.m: 2-D scatter plot of the latent points.
lvmScatterPlotColor.m: 2-D scatter plot of the latent points with color - for Swiss Roll data.
lvmTwoDPlot.m: Helper function for plotting the labels in 2-D.
mappingOptimise.m: Optimise the given model.
mlpCreate.m: Wrapper for NETLAB's mlp `net'.
mlpDisplay.m: Display the multi-layer perceptron model.
mlpExpandParam.m: Update mlp model with new vector of parameters.
mlpExtractParam.m: Wrapper for NETLAB's mlppak.
mlpOptimise.m: Optimise MLP for given inputs and outputs.
mlpOptions.m: Options for the multi-layered perceptron.
mlpOut.m: Output of an MLP model (wrapper for the NETLAB function mlpfwd).
mlpOutputGrad.m: Evaluate derivatives of mlp model outputs with respect to parameters.
modelCreate.m: Create a model of the specified type.
modelDisplay.m: Display a text output of a model.
modelExpandParam.m: Update a model structure with parameters.
modelExtractParam.m: Extract the parameters of a model.
modelLogLikeGradients.m: Compute a model's gradients wrt log likelihood.
modelLogLikelihood.m: Compute a model log likelihood.
modelOptimise.m: Optimise the given model.
modelOut.m: Give the output of a model for given X.
modelOutputGrad.m: Compute derivatives with respect to params of model outputs.
modelSamp.m: Give a sample from a model for given X.
ppcaEmbed.m: Embed data set with probabilistic PCA.
rbfCreate.m: Wrapper for NETLAB's rbf `net'.
rbfDisplay.m: Display an RBF network.
rbfExpandParam.m: Update rbf model with new vector of parameters.
rbfExtractParam.m: Wrapper for NETLAB's rbfpak.
rbfOptimise.m: Optimise RBF for given inputs and outputs.
rbfOptions.m: Default options for RBF network.
rbfOut.m: Output of an RBF model (wrapper for the NETLAB function rbffwd).
rbfOutputGrad.m: Evaluate derivatives of rbf model outputs with respect to parameters.
spectrumModify.m: Helper code for visualisation of spectrum data.
spectrumVisualise.m: Helper code for showing an spectrum during 2-D visualisation.
vectorModify.m: Helper code for visualisation of vectorial data.
vectorVisualise.m:  Helper code for plotting a vector during 2-D visualisation.
