MLTOOLS software
Version 0.1291		Sunday 05 Oct 2008 at 23:06

This toolbox provides various machine learning tools, either through wrapping other toolboxes (such as NETLAB) or providing the tool directly. It was designed originally as part of splitting the GPLVM and FGPLVM toolboxes.

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

modelLogLikeGradients.m: Compute a model's gradients wrt log likelihood.
mogTwoDPlot.m: Helper function for plotting the labels in 2-D.
mlpOutputGrad.m: Evaluate derivatives of mlp model outputs with respect to parameters.
rbfOut.m: Output of an RBF model.
dnetLogLikeGradients.m: Density network gradients.
demSwissRollLle4.m: Demonstrate LLE on the oil data.
mogEstep.m: Do an E-step on an MOG model.
demMppca1.m: Demonstrate MPPCA on a artificial dataset.
ppcaEmbed.m: Embed data set with probabilistic PCA.
linearLogLikelihood.m: Linear model log likelihood.
demSwissRollFullLle4.m: Demonstrate LLE on the oil data.
imageModify.m: Helper code for visualisation of image data.
rbfCreate.m: Wrapper for NETLAB's rbf `net'.
rbfExtractParam.m: Wrapper for NETLAB's rbfpak.
demSwissRollLle2.m: Demonstrate LLE on the oil data.
mlpDisplay.m: Display the multi-layer perceptron model.
lvmVisualise.m: Visualise the manifold.
lleOptimise.m: Optimise an LLE model.
rbfOptions.m: Default options for RBF network.
mogLowerBound.m: Computes lower bound on log likelihood for an MOG model.
mogSample.m: Sample from a mixture of Gaussians model.
dnetOut.m: Output of an DNET model.
modelOutputGradX.m: Compute derivatives with respect to model inputs of model outputs.
mlpLogLikeHessian.m: Multi-layer perceptron Hessian.
linearOutputGrad.m: Evaluate derivatives of linear model outputs with respect to parameters.
kbrParamInit.m: KBR model parameter initialisation.
vectorVisualise.m:  Helper code for plotting a vector during 2-D visualisation.
doubleMatrixReadFromFID.m: Read a full matrix from an FID.
modelOut.m: Give the output of a model for given X.
kpcaEmbed.m: Embed data set with kernel PCA.
multimodelLogLikeGradients.m: Gradient of MULTIMODEL model log likelihood with respect to parameters.
mlpOutputGradX.m: Evaluate derivatives of mlp model outputs with respect to inputs.
demSwissRollFullLle3.m: Demonstrate LLE on the oil data.
modelObjective.m: Objective function to minimise for given model.
dnetGradient.m: Density Network gradient wrapper.
dnetPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
kbrCreate.m: Create a KBR model.
rbfperiodicOptions.m: Create a default options structure for the RBFPERIODIC model.
mogOptimise.m: Optimise an MOG model.
dnetUpdateBeta.m: Do an M-step (update parameters) on an Density Network model.
mogCreate.m: Create a mixtures of Gaussians model.
imageVisualise.m: Helper code for showing an image during 2-D visualisation.
mlpLogLikeGradients.m: Multi-layer perceptron gradients.
linearOut.m: Obtain the output of the linear model.
lvmClassVisualise.m: Callback function for visualising data in 2-D.
dnetOutputGrad.m: Evaluate derivatives of dnet model outputs with respect to parameters.
rbfDisplay.m: Display an RBF network.
lvmTwoDPlot.m: Helper function for plotting the labels in 2-D.
mlpExtractParam.m: Extract weights and biases from an MLP.
multimodelLogLikelihood.m: Log likelihood of MULTIMODEL model.
mogMeanCov.m: Project a mixture of Gaussians to a low dimensional space.
demSwissRollLle3.m: Demonstrate LLE on the oil data.
ppcaOptions.m: Options for probabilistic PCA.
modelExpandParam.m: Update a model structure with parameters.
modelReadFromFID.m: Load from a FID produced by C++ code.
lvmScatterPlot.m: 2-D scatter plot of the latent points.
springDampersVisualise.m: Helper code for showing an spring dampers during 2-D visualisation.
kbrExtractParam.m: Extract parameters from the KBR model structure.
rbfOutputGradX.m: Evaluate derivatives of a RBF model's output with respect to inputs.
rbfperiodicParamInit.m: RBFPERIODIC model parameter initialisation.
rbfperiodicCreate.m: Create a RBFPERIODIC model.
mogPrintPlot.m: Print projection of MOG into two dimensions.
mvuEmbed.m: Embed data set with MVU.
rbfperiodicLogLikelihood.m: Log likelihood of RBFPERIODIC model.
multimodelDisplay.m: Display parameters of the MULTIMODEL model.
rbfExpandParam.m: Update rbf model with new vector of parameters.
modelParamInit.m: Initialise the parameters of the model.
kbrDisplay.m: Display parameters of the KBR model.
ppcaCreate.m: Density network model.
linearLogLikeGradients.m: Linear model gradients.
lleCreate.m: Locally linear embedding model.
dnetEstep.m: Do an E-step (update importance weights) on an Density Network model.
rbfperiodicExpandParam.m: Create model structure from RBFPERIODIC model's parameters.
mogUpdateMean.m: Update the means of an MOG model.
modelTieParam.m: Tie parameters of a model together.
demSwissRollLle1.m: Demonstrate LLE on the oil data.
rbfperiodicDisplay.m: Display parameters of the RBFPERIODIC model.
mlpOptimise.m: Optimise MLP for given inputs and outputs.
lvmPrintPlot.m: Print latent space for learnt model.
modelDisplay.m: Display a text output of a model.
springDampersModify.m: Helper code for visualisation of springDamper data.
linearExpandParam.m: Update linear model with vector of parameters.
modelLogLikelihood.m: Compute a model log likelihood.
lvmScatterPlotColor.m: 2-D scatter plot of the latent points with color.
dnetUpdateOutputWeights.m: Do an M-step (update parameters) on an Density Network model.
doubleMatrixWriteToFID.m: Writes a double matrix to an FID.
rbfperiodicExtractParam.m: Extract parameters from the RBFPERIODIC model structure.
mlpLogLikelihood.m: Multi-layer perceptron log likelihood.
mogOptions.m: Sets the default options structure for MOG models.
findNeighbours.m: find the k nearest neighbours for each point in Y.
modelHessian.m: Hessian of error function to minimise for given model.
smallrandEmbed.m: Embed data set with small random values.
rbfperiodicOutputGrad.m: Evaluate derivatives of RBFPERIODIC model outputs with respect to parameters.
kbrOut.m: Compute the output of a KBR model given the structure and input X.
demOilLle2.m: Demonstrate LLE on the oil data.
mapmodelReadFromFID.m: Load from a FID produced by C++ code.
lvmClassVisualisePath.m: Latent variable model path drawing in latent space.
lvmScatterPlotNeighbours.m: 2-D scatter plot of the latent points with neighbourhood.
modelExtractParam.m: Extract the parameters of a model.
kbrOptions.m: Create a default options structure for the KBR model.
mogLogLikelihood.m: Mixture of Gaussian's log likelihood.
demOilLle1.m: Demonstrate LLE on the oil data.
modelPointLogLikelihood.m: Compute the log likelihood of a given point.
dnetLogLikelihood.m: Density network log likelihood.
kbrOutputGrad.m: Evaluate derivatives of KBR model outputs with respect to parameters.
modelGradientCheck.m: Check gradients of given model.
modelGradient.m: Gradient of error function to minimise for given model.
modelOptimise.m: Optimise the given model.
modelCreate.m: Create a model of the specified type.
vectorModify.m: Helper code for visualisation of vectorial data.
mappingOptimise.m: Optimise the given model.
kbrOptimise.m: Optimise a KBR model.
dnetExtractParam.m: Extract weights and biases from an DNET.
dnetLowerBound.m: Computes lower bound on log likelihood for an DNET model.
rbfperiodicLogLikeGradients.m: Gradient of RBFPERIODIC model log likelihood with respect to parameters.
mogProject.m: Project a mixture of Gaussians to a low dimensional space.
mogUpdatePrior.m: Update the priors of an MOG model.
modelAddDynamics.m: Add a dynamics kernel to the model.
mlpOut.m: Output of an MLP model.
lvmLoadResult.m: Load a previously saved result.
dnetOptimise.m: Optimise an DNET model.
mlpParamInit.m: Initialise the parameters of an MLP model.
linearOptimise.m: Optimise a linear model.
mogUpdateCovariance.m: Update the covariances of an MOG model.
multimodelExpandParam.m: Create model structure from MULTIMODEL model's parameters.
linearExtractParam.m: Extract weights from a linear model.
linearCreate.m: Create a linear model.
viterbiAlign.m: Compute the Viterbi alignment.
multimodelCreate.m: Create a MULTIMODEL model.
dnetExpandParam.m: Update dnet model with new vector of parameters.
modelSetOutputWeights.m: Wrapper function to return set output weight and bias matrices.
modelOutputGrad.m: Compute derivatives with respect to params of model outputs.
dnetOutputGradX.m: Evaluate derivatives of DNET model outputs with respect to inputs.
spectrumModify.m: Helper code for visualisation of spectrum data.
lleOptions.m: Options for a density network.
demSwissRollFullLle2.m: Demonstrate LLE on the oil data.
modelSamp.m: Give a sample from a model for given X.
ppcaPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
modelReadFromFile.m: Read model from a file FID produced by the C++ implementation.
kbrExpandParam.m: Create model structure from KBR model's parameters.
dnetObjective.m: Wrapper function for Density Network objective.
rbfperiodicOut.m: Compute the output of a RBFPERIODIC model given the structure and input X.
mltoolsToolboxes.m: Load in the relevant toolboxes for the MLTOOLS.
lvmResultsDynamic.m: Load a results file and visualise them.
demSwissRollFullLle1.m: Demonstrate LLE on the oil data.
mlpCreate.m: Multi-layer peceptron model.
rbfperiodicOutputGradX.m: Evaluate derivatives of a RBFPERIODIC model's output with respect to inputs.
linearOutputGradX.m: Evaluate derivatives of linear model outputs with respect to inputs.
isomapEmbed.m: Embed data set with Isomap.
modelWriteToFID.m: Write to a stream a given model.
linearParamInit.m: Initialise the parameters of an LINEAR model.
spectrumVisualise.m: Helper code for showing an spectrum during 2-D visualisation.
matrixReadFromFID.m: Read a matrix from an FID.
dnetTest.m: Test some settings for the density network.
multimodelParamInit.m: MULTIMODEL model parameter initialisation.
multimodelExtractParam.m: Extract parameters from the MULTIMODEL model structure.
modelOptions.m: Returns a default options structure for the given model.
mlpExpandParam.m: Update mlp model with new vector of parameters.
demOilLle3.m: Demonstrate LLE on the oil data.
lvmNearestNeighbour.m: Give the number of errors in latent space for 1 nearest neighbour.
modelTest.m: Run some tests on the specified model.
linearOptions.m: Options for learning a linear model.
demOilLle4.m: Demonstrate LLE on the oil data.
dnetOptions.m: Options for a density network.
lleEmbed.m: Embed data set with LLE.
rbfOutputGrad.m: Evaluate derivatives of rbf model outputs with respect to parameters.
modelGetOutputWeights.m: Wrapper function to return output weight and bias matrices.
multimodelOptions.m: Create a default options structure for the MULTIMODEL model.
lmvuEmbed.m: Embed data set with landmark MVU
dnetCreate.m: Density network model.
mlpOptions.m: Options for the multi-layered perceptron.
ppcaOut.m: Output of an PPCA model.
linearDisplay.m: Display a linear model.
rbfOptimise.m: Optimise RBF for given inputs and outputs.
