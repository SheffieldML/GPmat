MLTOOLS software
Version 0.1		Tuesday 29 Nov 2005 at 15:41
Copyright (c) 2005 Neil D. Lawrence

This toolbox provides various machine learning tools, either through wrapping other toolboxes (such as NETLAB) or providing the tool directly. It was designed originally as part of splitting the GPLVM and FGPLVM toolboxes.

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
kbrOut.m: Obtain the output of the kernel based regression model.
kbrOutputGrad.m: Evaluate derivatives of kernel based regression model outputs with respect to parameters.
kpcaEmbed.m: Embed data set with kernel PCA.
linearCreate.m: Create a linear model.
linearExpandParam.m: Update linear model with vector of parameters.
linearExtractParam.m: Extract weights from a linear model.
linearOptimise.m: Optimise a linear model.
linearOut.m: Obtain the output of the linear model.
linearOutputGrad.m: Evaluate derivatives of linear model outputs with respect to parameters.
lvmClassVisualise.m: Callback function for visualising data in 2-D.
lvmScatterPlot.m: 2-D scatter plot of the latent points.
lvmTwoDPlot.m: Helper function for plotting the labels in 2-D.
mappingOptimise.m: Optimise the given model.
mlpCreate.m: Wrapper for NETLAB's mlp `net'.
mlpExpandParam.m: Update mlp model with new vector of parameters.
mlpExtractParam.m: Wrapper for NETLAB's mlppak.
mlpOptimise.m: Optimise MLP for given inputs and outputs.
mlpOut.m: Output of an MLP model (wrapper for the NETLAB function mlpfwd).
mlpOutputGrad.m: Evaluate derivatives of mlp model outputs with respect to parameters.
modelExpandParam.m: Update a model structure with parameters.
modelExtractParam.m: Extract the parameters of a model.
modelOptimise.m: Optimise the given model.
modelOut.m: Give the output of a model for given X.
modelOutputGrad.m: Compute derivatives with respect to params of model outputs.
ppcaEmbed.m: Embed data set with Isomap.
rbfCreate.m: Wrapper for NETLAB's rbf `net'.
rbfExpandParam.m: Update rbf model with new vector of parameters.
rbfExtractParam.m: Wrapper for NETLAB's rbfpak.
rbfOptimise.m: Optimise RBF for given inputs and outputs.
rbfOut.m: Output of an RBF model (wrapper for the NETLAB function rbffwd).
rbfOutputGrad.m: Evaluate derivatives of rbf model outputs with respect to parameters.
spectrumModify.m: Helper code for visualisation of spectrum data.
spectrumVisualise.m: Helper code for showing an spectrum during 2-D visualisation.
vectorModify.m: Helper code for visualisation of vectorial data.
vectorVisualise.m:  Helper code for plotting a vector during 2-D visualisation.
