MLTOOLS software
Version 0.1311		Thursday 14 May 2009 at 15:54

This toolbox provides various machine learning tools, either through wrapping other toolboxes (such as NETLAB) or providing the tool directly. It was designed originally as part of splitting the GPLVM and FGPLVM toolboxes.

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


MATLAB Files
------------

Matlab files associated with the toolbox are:

demMppca1.m: Demonstrate MPPCA on a artificial dataset.
kbrDisplay.m: Display parameters of the KBR model.
kbrOutputGrad.m: Evaluate derivatives of KBR model outputs with respect to parameters.
kbrParamInit.m: KBR model parameter initialisation.
kpcaEmbed.m: Embed data set with kernel PCA.
linearCreate.m: Create a linear model.
linearDisplay.m: Display a linear model.
linearExpandParam.m: Update linear model with vector of parameters.
linearExtractParam.m: Extract weights from a linear model.
linearLogLikeGradients.m: Linear model gradients.
linearParamInit.m: Initialise the parameters of an LINEAR model.
lleEmbed.m: Embed data set with LLE.
lvmClassVisualise.m: Callback function for visualising data in 2-D.
lvmScatterPlot.m: 2-D scatter plot of the latent points.
lvmScatterPlotColor.m: 2-D scatter plot of the latent points with color.
lvmTwoDPlot.m: Helper function for plotting the labels in 2-D.
mappingOptimise.m: Optimise the given model.
kbrExpandParam.m: Create model structure from KBR model's parameters.
mlpCreate.m: Multi-layer peceptron model.
mlpExtractParam.m: Extract weights and biases from an MLP.
mlpLogLikeGradients.m: Multi-layer perceptron gradients.
mlpLogLikeHessian.m: Multi-layer perceptron Hessian.
mlpLogLikelihood.m: Multi-layer perceptron log likelihood.
mlpOptimise.m: Optimise MLP for given inputs and outputs.
mlpOptions.m: Options for the multi-layered perceptron.
mlpOut.m: Output of an MLP model.
mlpOutputGrad.m: Evaluate derivatives of mlp model outputs with respect to parameters.
mlpOutputGradX.m: Evaluate derivatives of mlp model outputs with respect to inputs.
mlpParamInit.m: Initialise the parameters of an MLP model.
mltoolsToolboxes.m: Load in the relevant toolboxes for the MLTOOLS.
modelOutputGrad.m: Compute derivatives with respect to params of model outputs.
modelOutputGradX.m: Compute derivatives with respect to model inputs of model outputs.
modelParamInit.m: Initialise the parameters of the model.
modelPointLogLikelihood.m: Compute the log likelihood of a given point.
modelSamp.m: Give a sample from a model for given X.
modelTest.m: Run some tests on the specified model.
mogEstep.m: Do an E-step on an MOG model.
mogLowerBound.m: Computes lower bound on log likelihood for an MOG model.
mogOptimise.m: Optimise an MOG model.
imageModify.m: Helper code for visualisation of image data.
mogOptions.m: Sets the default options structure for MOG models.
imageVisualise.m: Helper code for showing an image during 2-D visualisation.
isomapEmbed.m: Embed data set with Isomap.
kbrCreate.m: Create a KBR model.
kbrExtractParam.m: Extract parameters from the KBR model structure.
lmvuEmbed.m: Embed data set with landmark MVU
mogUpdateCovariance.m: Update the covariances of an MOG model.
mogUpdateMean.m: Update the means of an MOG model.
mogUpdatePrior.m: Update the priors of an MOG model.
multimodelCreate.m: Create a MULTIMODEL model.
multimodelDisplay.m: Display parameters of the MULTIMODEL model.
multimodelExpandParam.m: Create model structure from MULTIMODEL model's parameters.
kbrOptimise.m: Optimise a KBR model.
kbrOptions.m: Create a default options structure for the KBR model.
kbrOut.m: Compute the output of a KBR model given the structure and input X.
linearLogLikelihood.m: Linear model log likelihood.
multimodelExtractParam.m: Extract parameters from the MULTIMODEL model structure.
multimodelOptions.m: Create a default options structure for the MULTIMODEL model.
multimodelLogLikeGradients.m: Gradient of MULTIMODEL model log likelihood with respect to parameters.
multimodelLogLikelihood.m: Log likelihood of MULTIMODEL model.
linearOptimise.m: Optimise a linear model.
linearOptions.m: Options for learning a linear model.
linearOut.m: Obtain the output of the linear model.
linearOutputGrad.m: Evaluate derivatives of linear model outputs with respect to parameters.
linearOutputGradX.m: Evaluate derivatives of linear model outputs with respect to inputs.
mlpDisplay.m: Display the multi-layer perceptron model.
mlpExpandParam.m: Update mlp model with new vector of parameters.
modelCreate.m: Create a model of the specified type.
modelDisplay.m: Display a text output of a model.
modelExpandParam.m: Update a model structure with parameters.
modelExtractParam.m: Extract the parameters of a model.
multimodelParamInit.m: MULTIMODEL model parameter initialisation.
ppcaEmbed.m: Embed data set with probabilistic PCA.
rbfCreate.m: Wrapper for NETLAB's rbf `net'.
rbfDisplay.m: Display an RBF network.
rbfExpandParam.m: Update rbf model with new vector of parameters.
rbfExtractParam.m: Wrapper for NETLAB's rbfpak.
rbfOptimise.m: Optimise RBF for given inputs and outputs.
rbfOptions.m: Default options for RBF network.
modelGradient.m: Gradient of error function to minimise for given model.
modelGradientCheck.m: Check gradients of given model.
modelHessian.m: Hessian of error function to minimise for given model.
modelLogLikeGradients.m: Compute a model's gradients wrt log likelihood.
modelLogLikelihood.m: Compute a model log likelihood.
modelObjective.m: Objective function to minimise for given model.
modelOptimise.m: Optimise the given model.
modelOptions.m: Returns a default options structure for the given model.
modelOut.m: Give the output of a model for given X.
modelTieParam.m: Tie parameters of a model together.
rbfOut.m: Output of an RBF model.
rbfOutputGrad.m: Evaluate derivatives of rbf model outputs with respect to parameters.
mogCreate.m: Create a mixtures of Gaussians model.
rbfperiodicCreate.m: Create a RBFPERIODIC model.
rbfperiodicDisplay.m: Display parameters of the RBFPERIODIC model.
rbfperiodicExpandParam.m: Create model structure from RBFPERIODIC model's parameters.
rbfperiodicExtractParam.m: Extract parameters from the RBFPERIODIC model structure.
rbfperiodicOutputGrad.m: Evaluate derivatives of RBFPERIODIC model outputs with respect to parameters.
rbfperiodicLogLikeGradients.m: Gradient of RBFPERIODIC model log likelihood with respect to parameters.
rbfperiodicLogLikelihood.m: Log likelihood of RBFPERIODIC model.
rbfperiodicOptions.m: Create a default options structure for the RBFPERIODIC model.
rbfperiodicOut.m: Compute the output of a RBFPERIODIC model given the structure and input X.
rbfperiodicOutputGradX.m: Evaluate derivatives of a RBFPERIODIC model's output with respect to inputs.
rbfperiodicParamInit.m: RBFPERIODIC model parameter initialisation.
smallrandEmbed.m: Embed data set with small random values.
spectrumModify.m: Helper code for visualisation of spectrum data.
vectorModify.m: Helper code for visualisation of vectorial data.
spectrumVisualise.m: Helper code for showing an spectrum during 2-D visualisation.
vectorVisualise.m:  Helper code for plotting a vector during 2-D visualisation.
viterbiAlign.m: Compute the Viterbi alignment.
modelAddDynamics.m: Add a dynamics kernel to the model.
mvuEmbed.m: Embed data set with MVU.
lvmClassVisualisePath.m: Latent variable model path drawing in latent space.
springDampersModify.m: Helper code for visualisation of springDamper data.
lleCreate.m: Locally linear embedding model.
springDampersVisualise.m: Helper code for showing an spring dampers during 2-D visualisation.
mogLogLikelihood.m: Mixture of Gaussian's log likelihood.
modelReadFromFile.m: Read model from a file FID produced by the C++ implementation.
modelReadFromFID.m: Load from a FID produced by C++ code.
doubleMatrixWriteToFID.m: Writes a double matrix to an FID.
mapmodelReadFromFID.m: Load from a FID produced by C++ code.
dnetCreate.m: Density network model.
dnetEstep.m: Do an E-step (update importance weights) on an Density Network model.
doubleMatrixReadFromFID.m: Read a full matrix from an FID.
modelWriteToFID.m: Write to a stream a given model.
matrixReadFromFID.m: Read a matrix from an FID.
dnetLogLikeGradients.m: Density network gradients.
dnetOptimise.m: Optimise an DNET model.
dnetOptions.m: Options for a density network.
dnetOut.m: Output of an DNET model.
dnetOutputGrad.m: Evaluate derivatives of dnet model outputs with respect to parameters.
dnetOutputGradX.m: Evaluate derivatives of DNET model outputs with respect to inputs.
dnetPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
dnetExpandParam.m: Update dnet model with new vector of parameters.
dnetExtractParam.m: Extract weights and biases from an DNET.
dnetGradient.m: Density Network gradient wrapper.
dnetLogLikelihood.m: Density network log likelihood.
dnetLowerBound.m: Computes lower bound on log likelihood for an DNET model.
dnetObjective.m: Wrapper function for Density Network objective.
dnetUpdateBeta.m: Do an M-step (update parameters) on an Density Network model.
dnetUpdateOutputWeights.m: Do an M-step (update parameters) on an Density Network model.
dnetTest.m: Test some settings for the density network.
mogPrintPlot.m: Print projection of MOG into two dimensions.
lvmLoadResult.m: Load a previously saved result.
lvmNearestNeighbour.m: Give the number of errors in latent space for 1 nearest neighbour.
lvmPrintPlot.m: Print latent space for learnt model.
lvmResultsDynamic.m: Load a results file and visualise them.
lvmVisualise.m: Visualise the manifold.
modelGetOutputWeights.m: Wrapper function to return output weight and bias matrices.
modelSetOutputWeights.m: Wrapper function to return set output weight and bias matrices.
ppcaCreate.m: Density network model.
ppcaOptions.m: Options for probabilistic PCA.
ppcaOut.m: Output of an PPCA model.
mogMeanCov.m: Project a mixture of Gaussians to a low dimensional space.
mogTwoDPlot.m: Helper function for plotting the labels in 2-D.
ppcaPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
rbfOutputGradX.m: Evaluate derivatives of a RBF model's output with respect to inputs.
mogProject.m: Project a mixture of Gaussians to a low dimensional space.
paramNameRegularExpressionLookup.m: Returns the indices of the parameter containing the given regular expression.
mogSample.m: Sample from a mixture of Gaussians model.
paramNameReverseLookup.m: Returns the index of the parameter with the given name.
findNeighbours.m: find the k nearest neighbours for each point in Y.
lleOptimise.m: Optimise an LLE model.
lvmScatterPlotNeighbours.m: 2-D scatter plot of the latent points with neighbourhood.
lfmResultsDynamic.m: Load a results file and visualise them.
lfmVisualise.m: Visualise the outputs in a latent force model
demOilLle2.m: Demonstrate LLE on the oil data.
lleOptions.m: Options for a density network.
demOilLle1.m: Demonstrate LLE on the oil data.
demOilLle3.m: Demonstrate LLE on the oil data.
demOilLle4.m: Demonstrate LLE on the oil data.
demSwissRollFullLle2.m: Demonstrate LLE on the oil data.
demSwissRollFullLle3.m: Demonstrate LLE on the oil data.
demSwissRollFullLle4.m: Demonstrate LLE on the oil data.
demSwissRollLle1.m: Demonstrate LLE on the oil data.
demSwissRollLle2.m: Demonstrate LLE on the oil data.
demSwissRollLle3.m: Demonstrate LLE on the oil data.
demSwissRollLle4.m: Demonstrate LLE on the oil data.
lfmClassVisualise.m: Callback function to visualize LFM in 2D
demSwissRollFullLle1.m: Demonstrate LLE on the oil data.
