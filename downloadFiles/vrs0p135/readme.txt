MLTOOLS software
Version 0.135		Thursday 03 Jun 2010 at 23:20

This toolbox provides various machine learning tools, either through wrapping other toolboxes (such as NETLAB) or providing the tool directly. It was designed originally as part of splitting the GPLVM and FGPLVM toolboxes.

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


MATLAB Files
------------

Matlab files associated with the toolbox are:

mlpExtractParam.m: Extract weights and biases from an MLP.
mvuEmbed.m: Embed data set with MVU.
linearParamInit.m: Initialise the parameters of an LINEAR model.
mogProject.m: Project a mixture of Gaussians to a low dimensional space.
doubleMatrixWriteToFID.m: Writes a double matrix to an FID.
linearExpandParam.m: Update linear model with vector of parameters.
rbfOut.m: Output of an RBF model.
mvuOptimise.m: Optimise an MVU model.
vectorModify.m: Helper code for visualisation of vectorial data.
modelHessian.m: Hessian of error function to minimise for given model.
leCreate.m: Laplacian eigenmap model.
dnetTest.m: Test some settings for the density network.
dnetObjective.m: Wrapper function for Density Network objective.
isomapReconstruct.m: Reconstruct an isomap form component parts.
rbfperiodicOutputGrad.m: Evaluate derivatives of RBFPERIODIC model outputs with respect to parameters.
mlpLogLikeGradients.m: Multi-layer perceptron gradients.
rbfDisplay.m: Display an RBF network.
rbfperiodicCreate.m: Create a RBFPERIODIC model.
linearCreate.m: Create a linear model.
lfmResultsDynamic.m: Load a results file and visualise them.
ppcaPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
modelTieParam.m: Tie parameters of a model together.
viterbiAlign.m: Compute the Viterbi alignment.
modelPointLogLikelihood.m: Compute the log likelihood of a given point.
ppcaPosteriorVar.m: Mean and variances of the posterior at points given by X.
lvmScatterPlot.m: 2-D scatter plot of the latent points.
ppcaEmbed.m: Embed data set with probabilistic PCA.
multimodelLogLikelihood.m: Log likelihood of MULTIMODEL model.
multimodelExtractParam.m: Extract parameters from the MULTIMODEL model structure.
ppcaOptions.m: Options for probabilistic PCA.
kbrExtractParam.m: Extract parameters from the KBR model structure.
mogOptions.m: Sets the default options structure for MOG models.
linearOptimise.m: Optimise a linear model.
dnetGradient.m: Density Network gradient wrapper.
kpcaEmbed.m: Embed data set with kernel PCA.
kbrCreate.m: Create a KBR model.
mlpLogLikelihood.m: Multi-layer perceptron log likelihood.
demSwissRollFullLle4.m: Demonstrate LLE on the oil data.
demSwissRollLle4.m: Demonstrate LLE on the oil data.
dnetLoadResult.m: Load a previously saved result.
mlpExpandParam.m: Update mlp model with new vector of parameters.
lvmClickVisualise.m: Visualise the manifold using clicks.
linearLogLikelihood.m: Linear model log likelihood.
mlpLogLikeHessian.m: Multi-layer perceptron Hessian.
mogOptimise.m: Optimise an MOG model.
rbfperiodicParamInit.m: RBFPERIODIC model parameter initialisation.
demSwissRollFullLle1.m: Demonstrate LLE on the oil data.
rbfperiodicExtractParam.m: Extract parameters from the RBFPERIODIC model structure.
findNeighbours.m: find the k nearest neighbours for each point in Y.
multimodelExpandParam.m: Create model structure from MULTIMODEL model's parameters.
mogCreate.m: Create a mixtures of Gaussians model.
rbfExpandParam.m: Update rbf model with new vector of parameters.
rbfperiodicOutputGradX.m: Evaluate derivatives of a RBFPERIODIC model's output with respect to inputs.
spectrumModify.m: Helper code for visualisation of spectrum data.
mogTwoDPlot.m: Helper function for plotting the labels in 2-D.
mogUpdateMean.m: Update the means of an MOG model.
multimodelParamInit.m: MULTIMODEL model parameter initialisation.
multimodelCreate.m: Create a MULTIMODEL model.
demOilLle1.m: Demonstrate LLE on the oil data.
dnetOutputGradX.m: Evaluate derivatives of DNET model outputs with respect to inputs.
demOilLle4.m: Demonstrate LLE on the oil data.
mogEstep.m: Do an E-step on an MOG model.
lleEmbed.m: Embed data set with LLE.
kbrDisplay.m: Display parameters of the KBR model.
linearOutputGrad.m: Evaluate derivatives of linear model outputs with respect to parameters.
modelGetOutputWeights.m: Wrapper function to return output weight and bias matrices.
mvuDeconstruct.m: break MVU in pieces for saving.
spectrumVisualise.m: Helper code for showing an spectrum during 2-D visualisation.
leOptions.m: Options for a Laplacian eigenmaps.
rbfOptimise.m: Optimise RBF for given inputs and outputs.
isomapCreate.m: isomap embedding model.
modelLogLikelihood.m: Compute a model log likelihood.
multimodelDisplay.m: Display parameters of the MULTIMODEL model.
modelExpandParam.m: Update a model structure with parameters.
lvmScatterPlotColor.m: 2-D scatter plot of the latent points with color.
mogUpdateCovariance.m: Update the covariances of an MOG model.
mogUpdatePrior.m: Update the priors of an MOG model.
spectralUpdateLaplacian.m: Update the Laplacian using graph connections.
modelSetOutputWeights.m: Wrapper function to return set output weight and bias matrices.
modelPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
modelGradientCheck.m: Check gradients of given model.
springDampersVisualise.m: Helper code for showing an spring dampers during 2-D visualisation.
paramNameRegularExpressionLookup.m: Returns the indices of the parameter containing the given regular expression.
demMppca1.m: Demonstrate MPPCA on a artificial dataset.
demSwissRollLle1.m: Demonstrate LLE on the oil data.
mogMeanCov.m: Project a mixture of Gaussians to a low dimensional space.
isomapEmbed.m: Embed data set with Isomap.
dnetReconstruct.m: Reconstruct an DNET form component parts.
lleOptimise.m: Optimise an LLE model.
mlpOut.m: Output of an MLP model.
lmvuEmbed.m: Embed data set with landmark MVU
demSwissRollLle3.m: Demonstrate LLE on the oil data.
rbfCreate.m: Wrapper for NETLAB's rbf `net'.
rbfOptions.m: Default options for RBF network.
vectorVisualise.m:  Helper code for plotting a vector during 2-D visualisation.
isomapOptions.m: Options for a isomap.
ppcaReconstruct.m: Reconstruct an PPCA form component parts.
modelPosteriorVar.m: variances of the posterior at points given by X.
mltoolsToolboxes.m: Load in the relevant toolboxes for the MLTOOLS.
kbrOptimise.m: Optimise a KBR model.
mvuOptions.m: Options for a MVU.
lvmClassVisualise.m: Callback function for visualising data.
modelCreate.m: Create a model of the specified type.
modelDisplay.m: Display a text output of a model.
lvmPrintPlot.m: Print latent space for learnt model.
lleOptions.m: Options for a locally linear embedding.
rbfOutputGradX.m: Evaluate derivatives of a RBF model's output with respect to inputs.
lfmVisualise.m: Visualise the outputs in a latent force model
kbrOptions.m: Create a default options structure for the KBR model.
lfmClassVisualise.m: Callback function to visualize LFM in 2D
kbrOut.m: Compute the output of a KBR model given the structure and input X.
kbrExpandParam.m: Create model structure from KBR model's parameters.
modelOut.m: Give the output of a model for given X.
leOptimise.m: Optimise an LE model.
demOilLle2.m: Demonstrate LLE on the oil data.
dnetWriteResult.m: Write a DNET result.
mapmodelReadFromFID.m: Load from a FID produced by C++ code.
linearExtractParam.m: Extract weights from a linear model.
demSwissRollFullLle2.m: Demonstrate LLE on the oil data.
mlpOptimise.m: Optimise MLP for given inputs and outputs.
multimodelOptions.m: Create a default options structure for the MULTIMODEL model.
modelSamp.m: Give a sample from a model for given X.
spectralUpdateX.m: Update the latent representation for spectral model.
rbfperiodicLogLikeGradients.m: Gradient of RBFPERIODIC model log likelihood with respect to parameters.
modelOutputGrad.m: Compute derivatives with respect to params of model outputs.
modelWriteToFID.m: Write to a stream a given model.
distanceWarp.m: Dynamic Time Warping Algorithm
rbfperiodicDisplay.m: Display parameters of the RBFPERIODIC model.
dnetEstep.m: Do an E-step (update importance weights) on an Density Network model.
linearOutputGradX.m: Evaluate derivatives of linear model outputs with respect to inputs.
rbfperiodicLogLikelihood.m: Log likelihood of RBFPERIODIC model.
leReconstruct.m: Reconstruct an LE form component parts.
dnetUpdateOutputWeights.m: Do an M-step (update parameters) on an Density Network model.
imageVisualise.m: Helper code for showing an image during 2-D visualisation.
modelReadFromFile.m: Read model from a file FID produced by the C++ implementation.
isomapDeconstruct.m: break isomap in pieces for saving.
lvmScatterPlotNeighbours.m: 2-D scatter plot of the latent points with neighbourhood.
smallrandEmbed.m: Embed data set with small random values.
demOilLle3.m: Demonstrate LLE on the oil data.
mvuReconstruct.m: Reconstruct an MVU form component parts.
mogLogLikelihood.m: Mixture of Gaussian's log likelihood.
lleReconstruct.m: Reconstruct an LLE form component parts.
ppcaCreate.m: Density network model.
modelOutputGradX.m: Compute derivatives with respect to model inputs of model outputs.
dnetOptimise.m: Optimise an DNET model.
modelExtractParam.m: Extract the parameters of a model.
rbfperiodicOptions.m: Create a default options structure for the RBFPERIODIC model.
mlpCreate.m: Multi-layer peceptron model.
mogLowerBound.m: Computes lower bound on log likelihood for an MOG model.
modelParamInit.m: Initialise the parameters of the model.
ppcaDeconstruct.m: break PPCA in pieces for saving.
dnetExpandParam.m: Update dnet model with new vector of parameters.
rbfperiodicExpandParam.m: Create model structure from RBFPERIODIC model's parameters.
lvmVisualise.m: Visualise the manifold.
rbfOutputGrad.m: Evaluate derivatives of rbf model outputs with respect to parameters.
modelReadFromFID.m: Load from a FID produced by C++ code.
dnetOptions.m: Options for a density network.
paramNameReverseLookup.m: Returns the index of the parameter with the given name.
modelLoadResult.m: Load a previously saved result.
leDeconstruct.m: break LE in pieces for saving.
isomapOptimise.m: Optimise an ISOMAP model.
dnetOut.m: Output of an DNET model.
modelLogLikeGradients.m: Compute a model's gradients wrt log likelihood.
lleDeconstruct.m: break LLE in pieces for saving.
modelOptions.m: Returns a default options structure for the given model.
matrixReadFromFID.m: Read a matrix from an FID.
multimodelLogLikeGradients.m: Gradient of MULTIMODEL model log likelihood with respect to parameters.
modelAddDynamics.m: Add a dynamics kernel to the model.
dnetPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
mlpOptions.m: Options for the multi-layered perceptron.
lvmResultsDynamic.m: Load a results file and visualise them.
linearDisplay.m: Display a linear model.
kbrOutputGrad.m: Evaluate derivatives of KBR model outputs with respect to parameters.
lvmLoadResult.m: Load a previously saved result.
lvmResultsClick.m: Load a results file and visualise them with clicks
mappingOptimise.m: Optimise the given model.
linearOptions.m: Options for learning a linear model.
dnetExtractParam.m: Extract weights and biases from an DNET.
dnetCreate.m: Density network model.
mlpDisplay.m: Display the multi-layer perceptron model.
mlpOutputGrad.m: Evaluate derivatives of mlp model outputs with respect to parameters.
modelTest.m: Run some tests on the specified model.
ppcaOut.m: Output of an PPCA model.
demSwissRollLle2.m: Demonstrate LLE on the oil data.
lvmClassClickVisualise.m: Callback function for visualising data in 2-D with clicks.
lvmScoreModel.m: Score model with a GP log likelihood.
dnetLowerBound.m: Computes lower bound on log likelihood for an DNET model.
lvmClassVisualisePath.m: Latent variable model path drawing in latent space.
kbrParamInit.m: KBR model parameter initialisation.
mogPrintPlot.m: Print projection of MOG into two dimensions.
linearLogLikeGradients.m: Linear model gradients.
modelObjective.m: Objective function to minimise for given model.
lvmNearestNeighbour.m: Give the number of errors in latent space for 1 nearest neighbour.
modelWriteResult.m: Write a model to file.
lleCreate.m: Locally linear embedding model.
dnetDeconstruct.m: break DNET in pieces for saving.
springDampersModify.m: Helper code for visualisation of springDamper data.
dnetLogLikelihood.m: Density network log likelihood.
lvmThreeDPlot.m: Helper function for plotting the labels in 3-D.
lvmTwoDPlot.m: Helper function for plotting the labels in 2-D.
linearOut.m: Obtain the output of the linear model.
mlpOutputGradX.m: Evaluate derivatives of mlp model outputs with respect to inputs.
modelGradient.m: Gradient of error function to minimise for given model.
dnetOutputGrad.m: Evaluate derivatives of dnet model outputs with respect to parameters.
lvmSetPlot.m: Sets up the plot for visualization of the latent space.
rbfperiodicOut.m: Compute the output of a RBFPERIODIC model given the structure and input X.
modelOptimise.m: Optimise the given model.
mogSample.m: Sample from a mixture of Gaussians model.
demSwissRollFullLle3.m: Demonstrate LLE on the oil data.
dnetUpdateBeta.m: Do an M-step (update parameters) on an Density Network model.
mlpParamInit.m: Initialise the parameters of an MLP model.
dnetLogLikeGradients.m: Density network gradients.
rbfExtractParam.m: Wrapper for NETLAB's rbfpak.
imageModify.m: Helper code for visualisation of image data.
doubleMatrixReadFromFID.m: Read a full matrix from an FID.
mvuCreate.m: Maximum variance unfolding embedding model.
