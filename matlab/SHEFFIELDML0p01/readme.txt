SHEFFIELDML software
Version 0.01		Wednesday 20 Feb 2013 at 09:48

This is the combined Sheffield ML toolbox. It contains MATLAB software for Gaussian processes.

Version 0.1
-----------

The first version which is spun out of the FGPLVM toolbox. The
corresponding FGPLVM toolbox is 0.15.


MATLAB Files
------------

Matlab files associated with the toolbox are:

sqexpKernExtractParam.m: Extract parameters from the SQEXP kernel structure.
noneKernExpandParam.m: Create kernel structure from NONE kernel's parameters.
disimKernGradX.m: Gradient of DISIM kernel with respect to a point x.
demCmu35gplvmReconstruct.m: Reconstruct right leg and body of CMU 35.
rbfperiodicKernExpandParam.m: Create kernel structure from RBFPERIODIC kernel's parameters.
dnetOutputGradX.m: Evaluate derivatives of DNET model outputs with respect to inputs.
fgplvmAddConstraint.m: Add latent constraints to FGPLVM model
modelLatentGradients.m: Gradients of the latent variables for dynamics models in the GPLVM.
readDoubleFromFID.m: Read a double from an FID.
matern32KernCompute.m: Compute the MATERN32 kernel given the parameters and X.
linardKernGradient.m: Gradient of LINARD kernel's parameters.
translateKernParamInit.m: TRANSLATE kernel parameter initialisation.
whitefixedXwhitefixedKernCompute.m: Compute a cross kernel between two WHITEFIXED kernels.
sdlfmaXsdrbfKernCompute.m: Cross kernel between a SDLFMA and a SDRBF kernels.
lleOptions.m: Options for a locally linear embedding.
polyardKernExtractParam.m: Extract parameters from the POLYARD kernel structure.
gpReversibleDynamicsLogLikeGradients.m: Gradients of the GP reversible dynamics wrt parameters.
whitehKernDiagGradX.m: Gradient of WHITEH kernel's diagonal with respect to X.
ggXgaussianKernGradient.m: Compute gradient between the GG and GAUSSIAN kernels.
spectralUpdateLaplacian.m: Update the Laplacian using graph connections.
lmcKernParamInit.m: LMC kernel parameter initialisation.
noneKernDiagGradX.m: Gradient of NONE kernel's diagonal with respect to X.
rbfwhiteKernExpandParam.m: Create kernel structure from RBF-WHITE kernel's
readBinaryDoubles.m: Read information from a binary file in as doubles.
gaussianPriorExtractParam.m: Extract params from Gaussian prior structure.
acclaimSemiLoadChannels.m: Load the channels from an AMC file for a subset
fgplvmExtractParam.m: Extract a parameter vector from a GP-LVM model.
lmcKernDiagGradient.m: Gradient of the LMC kernel's diagonal wrt parameters.
simwhiteXwhiteKernGradient.m: Compute gradient between the SIM-WHITE and WHITE kernels.
modelOut.m: Give the output of a model for given X.
rbfperiodicOutputGradX.m: Evaluate derivatives of a RBFPERIODIC model's output with respect to inputs.
noiseParamInit.m: Noise model's parameter initialisation.
cumGaussian.m: Cumulative distribution for Gaussian.
mlpOutputGrad.m: Evaluate derivatives of mlp model outputs with respect to parameters.
lfmGradientSigmaH3AA.m: Gradient of the function h_i(z) with respect \sigma.
plotMatrixOptions.m: Default options for plot matrix.
lmcKernExpandParam.m: Expands parameters into a LMC kernel structure.
leDeconstruct.m: break LE in pieces for saving.
gibbsKernGradX.m: Gradient of GIBBS kernel with respect to input locations.
robOneDynamicsCreate.m: Create the dynamics model. 
gpGradient.m: Gradient wrapper for a GP model.
ivm3dPlot.m: Make a 3-D or contour plot of the IVM.
spectrumModify.m: Helper code for visualisation of spectrum data.
traceProduct.m: Returns the trace of the product of two matrices.
fileKernGradient.m: Gradient of FILE kernel's parameters.
demRegressionOneIvm1.m: The data-set is sampled from a GP with known parameters.
fgplvmTestMissing.m: Make sure missing data likelihood match full ones.
velotransKernDisplay.m: Display parameters of the VELOTRANS kernel.
tensorKernExtractParam.m: Extract parameters from the TENSOR kernel structure.
ivmAddPoint.m: Add a point into the IVM representation.
probitNoiseLikelihood.m: Likelihood of the data under the PROBIT noise model.
sqexpKernDiagCompute.m: Compute diagonal of SQEXP kernel.
ngaussNoiseExpandParam.m: Create noise structure from NGAUSS noise's parameters.
orderedGradX.m: Gradient wrt x of log-likelihood for Ordered categorical model.
wienerKernDiagGradX.m: Gradient of WIENER kernel's diagonal with respect to X.
ratquadKernDiagGradient.m: Compute the gradient of the RATQUAD kernel's diagonal wrt parameters.
componentKernReadParamsFromFID.m: Read a component based kernel from a C++ file.
gpReversibleDynamicsExtractParam.m: Extract parameters from the GP reversible dynamics model.
cmpndKernExpandParamTransformSettings.m: Create kernel structure from CMPND kernel's parameter transformation settings.
nddisimKernGradient.m: Gradient of NDDISIM kernel's parameters.
skelVisualise.m: For drawing a skel representation of 3-D data.
translateKernExtractParam.m: Extract parameters from the TRANSLATE kernel structure.
ivmApproxGradX.m: Returns the gradient of the approxmate log-likelihood wrt x.
sdlfmKernComputeConstant.m: Compute constants for the SDLFM kernel
mgaussianNoiseExpandParam.m: Create noise structure from MGAUSSIAN noise's parameters.
rbfperiodicKernGradX.m: Gradient of RBFPERIODIC kernel with respect to a point x.
gaussianNoiseSites.m: Update the site parameters for the GAUSSIAN noise mode.
robOneDynamicsExtractParam.m: Extract parameters from the robot one dynamics model.
linard2KernDisplay.m: Display parameters of the LINARD2 kernel.
kernDisplay.m: Display the parameters of the kernel.
mlpardKernGradient.m: Gradient of MLPARD kernel's parameters.
rbfinfwhiteKernExtractParam.m: Extract parameters from the RBF-WHITE kernel
rbfardjitKernExpandParam.m: Create kernel structure from RBFARDJIT kernel's parameters.
stringSplit.m: Return separate parts of a string.
sdlfmXsdlfmvKernComputeBlock.m: Computes SDLFM kernel matrix for block i,j
biasKernParamInit.m: BIAS kernel parameter initialisation.
kernGetVariance.m: Get the signal associated with a the kernel.
demVowelsIsomap.m: Model the vowels data with a 2-D FGPLVM using RBF kernel.
whiteKernExtractParam.m: Extract parameters from the WHITE kernel structure.
sheatKernCompute.m: Compute a cross kernel between two SHEAT kernels.
disimKernGradient.m: Gradient of DISIM kernel's parameters.
sqexpKernExpandParam.m: Create kernel structure from SQEXP kernel's parameters.
mlpardKernExpandParam.m: Create kernel structure from MLPARD kernel's parameters.
uniformPriorGradient.m: Gradient wrt x of the uniform prior.
rbfard2KernExtractParam.m: Extract parameters from the RBFARD2 kernel structure.
laplacePriorLogProb.m: Log probability of Laplace prior.
ppcaEmbed.m: Embed data set with probabilistic PCA.
ivmGradX.m: Returns the gradient of the log-likelihood wrt x.
demRobotWirelessFgplvm4.m: Wireless Robot data from University of Washington with dynamics and back constraints.
noiseWriteToFID.m: Load from an FID written by the C++ implementation.
lnDiffErfs.m: Helper function for computing the log of difference
multiKernGradient.m: Gradient of MULTI kernel's parameters.
ppcaOut.m: Output of an PPCA model.
kernCorrelation.m: Compute the correlation matrix kernel given the parameters and X.
laplacePriorExtractParam.m: Extract params from Laplace prior structure.
gpReconstruct.m: Reconstruct an GP form component parts.
xyzhumanevaJoint2pos.m:
bvhWriteFile.m: Write a bvh file from a given structure and channels.
mlpParamInit.m: Initialise the parameters of an MLP model.
fgplvmTaylorAngleErrors.m: Helper function for computing angle errors for CMU 35 data.
skel2xyz.m: Compute XYZ values given skeleton structure and channels.
heatKernDisplay.m: Display parameters of the HEAT kernel.
rbfwhiteXwhiteKernGradient.m: Compute gradient between the RBF-WHITE and
gpDynamicsExtractParam.m: Extract parameters from the GP dynamics model.
treeGetWidths.m: give width of each level of tree.
simKernCompute.m: Compute the SIM kernel given the parameters and X.
fgplvmPrintPlot.m: Print latent space for learnt model.
lfmComputeH4AV.m: Helper function for computing part of the LFMAV kernel.
findDirectedNeighbours.m: find the k nearest neighbours for each point in Y preventing cycles in the graph.
disimSample.m: Sample from SIM kernel
fgplvmScatterPlotColor.m: 2-D scatter plot of the latent points with color - for Swiss Roll data.
gpDisplay.m: Display a Gaussian process model.
fgplvmPointGradient.m: Wrapper function for gradient of a single point.
rbfKernDiagGradX.m: Gradient of RBF kernel's diagonal with respect to X.
cmpndNoiseGradientParam.m: Gradient of CMPND noise's parameters.
dnetReconstruct.m: Reconstruct an DNET form component parts.
gaussianKernExpandParam.m: Create kernel structure from gaussian kernel's parameters.
demStickFgplvm1.m: Model the stick man using an RBF kernel.
printPlot.m: Print a plot to eps and png files.
wangPriorExtractParam.m: Extract params from Wang prior structure.
rbfperiodicLogLikeGradients.m: Gradient of RBFPERIODIC model log likelihood with respect to parameters.
diagKernGradient.m: Gradient of DIAG kernel's parameters.
matrixReadFromFID.m: Read a matrix from an FID.
rbfKernDiagCompute.m: Compute diagonal of RBF kernel.
srbfhKernCompute.m: Compute an SRBFH kernel.
ggXgaussianKernCompute.m: Compute a cross kernel between the GG and GAUSSIAN kernels.
simXsimKernDiagCompute.m: Diagonal of a cross kernel between two SIM kernels.
gaussianNoiseGradVals.m: Gradient of GAUSSIAN noise log Z with respect to input mean and variance.
whiteKernDiagGradient.m: Compute the gradient of the WHITE kernel's diagonal wrt parameters.
ngaussian.m: Compute a Gaussian with mean 0 and variance 1.
ngaussNoiseLikelihood.m: Likelihood of the data under the NGAUSS noise model.
xyzpoppeAnim.m: Animate point cloud of stick man from Poppe dataset.
simKernDisplay.m: Display parameters of the SIM kernel.
pcaEmbed.m: Embed data set with PCA.
lfmLogLikeGradients.m: Compute the gradients of the log likelihood of a LFM model.
demRegressionGp.m: Demonstrate Gaussian processes for regression.
imageVisualise.m: Helper code for showing an image during 2-D visualisation.
rbfDisplay.m: Display an RBF network.
invcmpndKernExpandParam.m: Create kernel structure from INVCMPND kernel's parameters.
demOilFgplvm5.m: Oil data with partially independent training conditional.
robThreeDynamicsCreate.m: Create the dynamics model. 
gpDynamicsSequenceLogLikeGradient.m: Log-likelihood gradient for of a sequence of the GP-LVM dynamics.
ratquadKernGradX.m: Gradient of RATQUAD kernel with respect to input locations.
disimKernDisplay.m: Display parameters of the DISIM kernel.
fgplvmPosteriorVar.m: Variances of the posterior at points given by X.
normuniPriorLogProb.m: Log probability of a normal uniform.
ivmLogLikelihood.m: Return the log-likelihood for the IVM.
lfmaXrbfKernGradient.m: Compute gradient between the LFMA and RBF kernels.
simwhiteKernDiagCompute.m: Compute the diagonal of the SIM-WHITE kernel.
scaleNoiseParamInit.m: Scale noise model's parameter initialisation.
kbrExpandParam.m: Create model structure from KBR model's parameters.
expKernParamInit.m: EXP kernel parameter initialisation.
probit3dPlot.m: Draw a 3D or contour plot for the probit.
demBrendanFgplvm1.m: Use the GP-LVM to model the Frey face data with FITC.
simwhiteKernDisplay.m: Display parameters of the SIM-WHITE kernel.
linard2KernParamInit.m: LINARD2 kernel parameter initialisation.
lvmVisualiseGeneral.m: Visualise the manifold.
multimodelDisplay.m: Display parameters of the MULTIMODEL model.
lfmGradientSigmaH3AV.m: Gradient of the function h_i(z) with respect \sigma.
sdlfmaXsdlfmaKernGradientBlock.m: Gradients of the parameters in block i,j
mlpOptimise.m: Optimise MLP for given inputs and outputs.
mlpOptions.m: Options for the multi-layered perceptron.
rbfOutputGrad.m: Evaluate derivatives of rbf model outputs with respect to parameters.
fileKernRead.m: Read kernel values from file or cache.
lmcKernDiagGradX.m: Gradient of LMC kernel's diagonal with respect to X.
simwhiteXrbfinfwhiteKernGradient.m: Compute a cross gradient between a
dnetDeconstruct.m: break DNET in pieces for saving.
robThreeDynamicsSetLatentValues.m: Set the latent values inside the model.
demFourWalksReconstruct.m: Reconstruct right leg of CMU 35.
demOilFgplvm7.m: Oil data with variational sparse approximation.
biasKernDiagGradient.m: Compute the gradient of the BIAS kernel's diagonal wrt parameters.
lfmjXlfmvKernCompute.m: Jolt and velocity LFM kernel  
demSilhouetteGp1.m: Model silhouette data with independent RBF GPs.
gaussianPriorLogProb.m: Log probability of Gaussian prior.
laplacePriorGradient.m: Gradient wrt x of the log Laplace prior.
whiteblockKernDiagGradient.m: WHITEBLOCK kernel's diagonal gradient wrt par.
ivmRun.m: Run the IVM on a given data set.
xyzpoppeDraw.m: Helper function for drawing data from Poppe.
lfmGradientUpsilon.m: Gradient of the function \upsilon(z) with respect to
heatKernGradient.m: Gradient of HEAT kernel's parameters.
matern32KernDisplay.m: Display parameters of the MATERN32 kernel.
fgplvmOptimise.m: Optimise the FGPLVM.
lfmVisualise.m: Visualise the outputs in a latent force model
gpReversibleDynamicsCreate.m: Create a reversible dynamics model. 
demSpgp1dGp2.m: Do a simple 1-D regression after Snelson & Ghahramani's example.
ngaussNoise3dPlot.m: Draws a 3D or contour plot for the NGAUSS noise model.
mlpCreate.m: Multi-layer peceptron model.
modelTest.m: Run some tests on the specified model.
fgplvmVisualise.m: Visualise the manifold.
lfmaXlfmKernGradient.m: Compute a cross gradient between a LFMA and a LFM.
identityTransform.m: 
cmdsRoadData.m: This script uses classical MDS to visualise some road distance data.
gaussianKernParamInit.m: Gaussian kernel parameter initialisation.
lfmjXlfmaKernCompute.m: Jolt and acceleration LFM kernel  
modelSamp.m: Give a sample from a model for given X.
lfmComputeH3.m: Helper function for computing part of the LFM kernel.
demCmu35gplvm4.m: Learn a GPLVM on CMU 35 data set.
wangPriorGradient.m: Gradient wrt x of the Wang prior.
multiKernExpandParam.m: Create kernel structure from MULTI kernel's parameters.
rbfard2KernDiagGradX.m: Gradient of RBFARD2 kernel's diagonal with respect to X.
scg2.m: Scaled conjugate gradient optimization like netlab's scg, with a slight modification for speeding it up.
ardKernGradient.m: Gradient of ARD kernel's parameters.
gpCovGradsTest.m: Test the gradients of the likelihood wrt the covariance.
whiteKernExpandParam.m: Create kernel structure from WHITE kernel's parameters.
xgamrnd.m: Draw a sample from the gamma distribution.
linardKernDiagCompute.m: Compute diagonal of LINARD kernel.
matern32KernParamInit.m: MATERN32 kernel parameter initialisation.
modelReadFromFID.m: Load from a FID produced by C++ code.
sdlfmKernCompute.m: Compute the SDLFM kernel given the parameters and X.
gpReadFromFID.m: Load from a FID produced by the C++ implementation.
rbfardjitKernGradX.m: Gradient of RBFARDJIT kernel with respect to input locations.
demFourWalks1.m: Model four seperate walsk using an RBF kernel and dynamics.
ratquadKernDiagCompute.m: Compute diagonal of RATQUAD kernel.
rbfXnoneKernCompute.m: Compute a cross kernel between RBF and NONE kernels.
lfmvXlfmKernGradient.m: Compute a cross gradient for a LFMVXLFM.
demSwissRollFullLle3.m: Demonstrate LLE on the oil data.
linard2KernGradX.m: Gradient of LINARD2 kernel with respect to input locations.
lfmvKernExtractParam.m: Extract parameters from the LFMV kernel structure.
kernPriorLogProb.m: Compute penalty terms associated with kernel priors.
gaussianNoisePointPlot.m: Plot the data-points for the GAUSSIAN noise model.
sdlfmvXsdlfmKernCompute.m: Cross kernel between a SDLFMV and a SDLFM kernels.
ncnmNoiseExpandParam.m: Expand null category noise model's structure from param vector.
mogPrintPlot.m: Print projection of MOG into two dimensions.
disimXdisimKernCompute.m: Compute a cross kernel between two DISIM kernels.
gaussianKernDiagCompute.m: Compute diagonal of gaussian kernel.
rbfwhiteKernParamInit.m: RBF-WHITE kernel parameter initialisation. The RBF-
demSilhouettePlot.m:
cmpndNoiseGradVals.m: Gradient of CMPND noise log Z with respect to input mean and variance.
demStickGp1.m: Demonstrate Gaussian processes for regression on stick man data.
rbfard2KernDiagCompute.m: Compute diagonal of RBFARD2 kernel.
ratquadKernParamInit.m: RATQUAD kernel parameter initialisation.
lfmComputeH4JP.m: Helper function for computing part of the LFMJP kernel.
mgaussianNoiseOut.m: Compute the output of the MGAUSSIAN noise given the input mean and variance.
modelCreate.m: Create a model of the specified type.
lnCumGaussSum.m: The log of the weighted sum of two cumulative Gaussians.
gpSubspaceExpandParam.m: 
mlpardKernDiagGradX.m: Gradient of MLPARD kernel's diagonal with respect to X.
writeIntToFID.m: Writes an integer to an FID.
sdlfmvKernDiagComputeBlock.m: Diagonal of a SDLFM kernel matrix for block i
xyzhumanevaHeadingAngle.m: 
dnetUpdateOutputWeights.m: Do an M-step (update parameters) on an Density Network model.
sdlfmXsdlfmKernGradient.m: Gradients of cross kernel between 2 SDLFM kernels.
leReconstruct.m: Reconstruct an LE form component parts.
rbfardKernDiagGradX.m: Gradient of RBFARD kernel's diagonal with respect to X.
gibbsKernDisplay.m: Display parameters of the GIBBS kernel.
mlpardKernCompute.m: Compute the MLPARD kernel given the parameters and X.
rbfinfwhiteXrbfinfwhiteKernGradient.m: Compute a cross gradient between two
kernFactors.m: Extract factors associated with transformed
ivmUpdateM.m: Update matrix M, L, v and mu.
lfmKernDiagGradX.m: Gradient of LFM kernel's diagonal with respect to X.
gpPosteriorGradMeanCovar.m: Gadient of the mean and variances of the posterior at points given by X.
lfmaaComputeUpsilonMatrix.m: Upsilon matrix acce. accel. with t1, t2 limits
priorWriteToFID.m: Write a prior to a C++ stream.
demEP1.m: Demonstrate Expectation propagation on a toy data set..
robTwoDynamicsExpandParam.m: Place the parameters vector into the model for first robot dynamics.
sdlfmvXsdrbfKernGradient.m: Gradients cross kernel between a SDLFM and SDRBF
printLatexOptions.m: Options for printing a plot to LaTeX.
lfmjXlfmKernCompute.m: Jolt and position LFM kernel  
stickVisualise.m: For drawing a stick representation of 3-D data.
multiKernDiagGradient.m: Compute the gradient of the MULTI kernel's diagonal wrt parameters.
simwhiteXrbfwhiteKernCompute.m: Compute a cross kernel between a SIM-WHITE
fgplvmPointSampleLogLikelihood.m:
modelTieParam.m: Tie parameters of a model together.
wienerKernCompute.m: Compute the WIENER kernel given the parameters and X.
gibbsperiodicKernGradient.m: Gradient of GIBBSPERIODIC kernel's parameters.
lvmResultsDynamic.m: Load a results file and visualise them.
ngaussNoiseGradientParam.m: Gradient of NGAUSS noise's parameters.
lfmapComputeUpsilonVector.m: Upsilon vector for acce. pos. with t1 limit
gaussianPriorExpandParam.m: Expand Gaussian prior structure from param vector.
gpDynamicsCreate.m: Create the dynamics model. 
nddisimKernExpandParam.m: Create kernel structure from NDDISIM kernel's parameters.
orderedNoiseExtractParam.m: Extract parameters from the ORDERED noise structure.
sdlfmKernExpandParam.m: Pass parameters from params to SDLFM kernel
tensorKernDiagGradX.m: Gradient of TENSOR kernel's diagonal with respect to X.
noiseUpdateNuG.m: Update nu and g for a given noise model.
ggwhiteKernExpandParam.m: Create kernel structure from GG white kernel's parameters.
ngaussNoiseNuG.m: Update nu and g parameters associated with noiseless Gaussian noise model.
acclaim2xyz.m: Compute XYZ values given skeleton structure and channels.
lfmComputeH4AP.m: Helper function for computing part of the LFMAP kernel.
sdlfmXsdlfmKernGradientBlockIGJ.m: 
rocholForeSub.m: Foreward substitute the representation of the rank one Cholesky.
disimXsimKernCompute.m: Compute a cross kernel between DISIM and SIM kernels.
sdrbfKernCompute.m: Compute the SDRBF kernel given the parameters and t1.
lvmScatterPlotNoVar.m: 2-D scatter plot of the latent points.
treeFindLeaves.m: Return indices of all leaf nodes in a tree structure.
leOptions.m: Options for a Laplacian eigenmaps.
rbfperiodicExtractParam.m: Extract parameters from the RBFPERIODIC model structure.
lmcKernExtractParam.m: Extract parameters from the LMC kernel struc.
ivmVirtual.m: Create virtual data points with the specified invariance.
whiteKernDisplay.m: Display parameters of the WHITE kernel.
simComputeTest.m: Test the file simComputeH.
rbfperiodicKernParamInit.m: RBFPERIODIC kernel parameter initialisation.
sdlfmvKernDisplay.m: Display parameters of the SDLFMV kernel.
ardKernCompute.m: Compute the ARD kernel given the parameters and X.
springDampersModify.m: Helper code for visualisation of springDamper data.
whitefixedKernDiagCompute.m: Compute diagonal of WHITEFIXED kernel.
parseWirelessData.m: Load wireless strength data.
fgplvmBackConstraintGrad.m: Gradient with respect to back constraints if present.
ncnmNoiseExtractParam.m: Extract parameters from null category noise model.
rbfard2KernDisplay.m: Display parameters of the RBFARD2 kernel.
xyzhumanevaAlign.m:
mogOptimise.m: Optimise an MOG model.
cmpndNoiseSites.m: Site updates for compound noise model.
ardKernExtractParam.m: Extract parameters from the ARD kernel structure.
indexardKernGradient.m: Gradient of INDEXARD kernel's parameters.
rbfperiodicLogLikelihood.m: Log likelihood of RBFPERIODIC model.
lfmGradientH41VP.m: Gradient of the function h_i(z) with respect to some of the
lfmaXlfmaKernCompute.m: Acceleration and acceleration LFM kernel  
sdlfmKernExtractParam.m: Extract parameters from the SDLFM kernel structure.
demRobotWirelessNavigate.m: Take some test data for the robot and navigate with it.
mocapParseText.m: Parse a motion capture text file.
rbfperiodicKernExtractParam.m: Extract parameters from the RBFPERIODIC kernel structure.
rbfKernDiagGradient.m: Compute the gradient of the RBF kernel's diagonal wrt parameters.
demCmu35gplvmFgplvm2.m: Learn a GPLVM on CMU 35 data set.
simKernDiagCompute.m: Compute diagonal of SIM kernel.
defaultOptions.m: The default options for optimisation.
ggKernParamInit.m: GG kernel parameter initialisation.
fgplvmPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
lfmComputeH3JA.m: Helper function for computing part of the LFMJA kernel.
cmpndKernGradX.m: Gradient of CMPND kernel with respect to a point x.
disimComputeHPrime.m: Helper function for comptuing part of the DISIM kernel.
stickModify.m: Helper code for visualisation of a stick man.
matern52KernGradX.m: Gradient of MATERN52 kernel with respect to input locations.
lvmVisualise.m: Visualise the manifold.
multimodelExtractParam.m: Extract parameters from the MULTIMODEL model structure.
normuniPriorGradient.m: Gradient wrt x of the log normal uniform prior.
multiKernExpandParamTransformSettings.m: Create kernel structure from MULTI kernel's parameter transformation settings.
gpOptimise.m: Optimise the inducing variable based kernel.
priorTest.m: Run some tests on the specified prior.
noiseOut.m: Give the output of the noise model given the mean and variance.
sdlfmXsdrbfKernGradientBlockIEJ.m: 
rbfperiodic2KernExtractParam.m: Extract parameters from the RBFPERIODIC2 kernel structure.
sdlfmaXsdrbfKernGradient.m: Gradients cross kernel between a SDLFM and SDRBF
smoothAngleChannels.m: Try and remove artificial discontinuities associated with angles.
orderedNoisePointPlot.m: Plot the data-points for the ORDERED noise model.
rbfardKernDisplay.m: Display parameters of the RBFARD kernel.
hessianCheck.m: Check Hessian of objective function.
invcmpndKernGradient.m: Gradient of INVERSE-PRESICION-CMPND kernel's parameters.
scaleNoiseExpandParam.m: Expand Scale noise structure from param vector.
lfmavComputeUpsilonMatrix.m: Upsilon matrix acce. vel. with t1, t2 limits
lfmComputeH3VV.m: Helper function for computing part of the LFMVXLFMV kernel.
whiteblockKernExtractParam.m: Extract parameters from WHITEBLOCK kernel str.
fgplvmPointObjectiveGradient.m: Wrapper function for objective and gradient of a single point in latent space and the output location..
gpCovGrads.m: Sparse objective function gradients wrt Covariance functions for inducing variables.
mocapResultsCppBvh.m: Load results from cpp file and visualise as a bvh format.
simwhiteKernCompute.m: Compute the SIM-WHITE kernel given the parameters, t1
ivmComputeLandM.m: Compute the L and M matrix.
gpExtractParam.m: Extract a parameter vector from a GP model.
polyKernExtractParam.m: Extract parameters from the POLY kernel structure.
robTwoDynamicsCreate.m: Create the dynamics model. 
orderedNoiseGradVals.m: Gradient of ORDERED noise log Z with respect to input mean and variance.
probitNoiseGradVals.m: Gradient of PROBIT noise log Z with respect to input mean and variance.
polyKernCompute.m: Compute the POLY kernel given the parameters and X.
rbfperiodicOptions.m: Create a default options structure for the RBFPERIODIC model.
whitehKernDisplay.m: Display parameters of the WHITEH kernel.
demCmu35gplvmReconstructTaylor.m: Reconstruct right leg of CMU 35.
xyzModify.m: Update visualisation of skeleton data.
matern32KernDiagGradX.m: Gradient of MATERN32 kernel's diagonal with respect to X.
rbfinfwhiteKernGradX.m: Gradient of RBF-WHITE kernel (with integration limits
lfmjvComputeUpsilonMatrix.m: Upsilon matrix jolt. vel. with t1, t2 limits
ivmOptimiseKernel.m: Optimise the kernel parameters.
sdlfmjMeanCompute.m: Jolt mean for the switching dynamical LFM model.
robOneDynamicsSetLatentValues.m: Set the latent values inside the model.
ivmEpUpdateM.m: Update matrix M, L, varSigma and mu for EP.
multiKernExtractParam.m: Extract parameters from the MULTI kernel structure.
kernGradX.m: Compute the gradient of the kernel wrt X.
linearOut.m: Obtain the output of the linear model.
ggwhiteXwhiteKernGradient.m: Compute gradient between the GGWHITE and WHITE kernels.
ivmComputeInfoChange.m: Compute the information change associated with each point.
sheatKernDiagGradient.m: Gradient of the parameters of diagonal of a SHEAT kernel.
lfmwhiteXrbfwhiteKernGradient.m: Compute a cross gradient between an LFM-WHITE
lfmaaGradientSigmaUpsilonMatrix.m: Gradient of upsilon matrix aa wrt sigma
ivmContour.m: Special contour plot showing decision boundary.
translateKernGradX.m: Gradient of TRANSLATE kernel with respect to a point x.
sqexpKernParamInit.m: SQEXP kernel parameter initialisation.
fgplvmAddDynamics.m: Add a dynamics kernel to the model.
rbfardjitKernCompute.m: Compute the RBFARD kernel given the parameters and X.
ivmDowndateM.m: Remove point from M, L, mu and varSigma.
lmcKernDisplay.m: Display parameters of the LMC kernel.
writeStringToFID.m: Writes a string to an FID.
lmcKernGradient.m: Gradient of LMC kernel's parameters.
demOilFgplvm9.m: Oil data with three dimensions and variational sparse approximation.
xyzankurAnimCompareMultiple.m: Animate many predictions and ground truth for
multiKernTest.m: Run some tests on the multiple output block kernel.
whiteKernDiagGradX.m: Gradient of WHITE kernel's diagonal with respect to X.
demOilLle1.m: Demonstrate LLE on the oil data.
linardKernDisplay.m: Display parameters of the LINARD kernel.
sdlfmXsdlfmvKernGradientBlockIGJ.m: 
whitehKernGradient.m: Gradient of WHITEH kernel's parameters.
rocholMultiply.m: Multiply by the rank one Cholesky.
gradLnDiffErfs.m: Compute the gradient of the log difference of two erfs.
biasKernDisplay.m: Display parameters of the BIASkernel.
linearOutputGrad.m: Evaluate derivatives of linear model outputs with respect to parameters.
matern32KernDiagGradient.m: Compute the gradient of the MATERN32 kernel's diagonal wrt parameters.
gammaPriorParamInit.m: Gamma prior model's parameter initialisation.
swissRollScatter.m: 3-D scatter plot with colors.
sdrbfKernExpandParam.m: Pass parameters from params to SDRBF kernel
fgplvmDynamicsPlot.m: 2-D scatter plot of the latent points.
demSwissRollLle2.m: Demonstrate LLE on the oil data.
demMppca1.m: Demonstrate MPPCA on a artificial dataset.
demOilFgplvm3.m: Oil data with deterministic training conditional.
generateCrescentData.m: Generate crescent data.
rbfardjitKernDiagCompute.m: Compute diagonal of RBFARDJIT kernel.
gibbsKernDiagCompute.m: Compute diagonal of GIBBS kernel.
deg2rad.m: Transform degrees to radians.
linearDisplay.m: Display a linear model.
demOilLle4.m: Demonstrate LLE on the oil data.
lfmComputeH3AA.m: Helper function for computing part of the LFMAA kernel.
wienerKernGradient.m: Gradient of WIENER kernel's parameters.
lfmaKernParamInit.m: LFMA kernel parameter initialisation. 
lfmGradientH42.m: Gradient of the function h_i(z) with respect to some of the
gaussianwhiteKernDisplay.m: Display parameters of the GAUSSIAN white kernel.
fgplvmResultsDynamic.m: Load a results file and visualise them.
negLogLogit.m: Function which returns the negative log of the logistic function.
ncnmLoadData.m: Load a dataset.
rbfperiodic2KernCompute.m: Compute the RBFPERIODIC2 kernel given the parameters and X.
rbfperiodicOutputGrad.m: Evaluate derivatives of RBFPERIODIC model outputs with respect to parameters.
ggwhiteKernDiagCompute.m: Compute diagonal of GG WHITE kernel.
bvhVisualise.m: For updating a bvh representation of 3-D data.
gaussianwhiteKernExpandParam.m: Create kernel structure from gaussian white 
whitefixedKernParamInit.m: WHITEFIXED kernel parameter initialisation.
lfmapGradientSigmaUpsilonVector.m: Gradient of upsilon vector ap wrt sigma
modelOutputGradX.m: Compute derivatives with respect to model inputs of model outputs.
rbfwhiteKernDiagCompute.m: Compute diagonal of RBF-WHITE kernel.
fgplvmDeconstruct.m: break FGPLVM in pieces for saving.
ggwhiteXgaussianwhiteKernGradX.m: Compute gradient between the GG white and
sdlfmXsdrbfKernCompute.m: Cross kernel between a SDLFM and a SDRBF kernels.
ngaussNoiseParamInit.m: NGAUSS noise parameter initialisation.
mogCreate.m: Create a mixtures of Gaussians model.
rbfardKernExpandParam.m: Create kernel structure from RBFARD kernel's parameters.
simwhiteKernGradX.m: Gradient of SIM-WHITE kernel with respect to a point t.
lfmwhiteXwhiteKernCompute.m: Compute a cross kernel between the LFM-WHITE
cmpndNoiseOut.m: Compute the output of the CMPND noise given the input mean and variance.
sdlfmaXsdlfmKernGradientBlockIEJ.m: 
expKernGradient.m: Gradient of EXP kernel's parameters.
nddisimKernDisplay.m: Display parameters of the NDDISIM kernel.
lfmvXrbfvKernCompute.m: Compute a cross kernel between the LFMV and RBFV kernels.
lfmGradientSigmaH4AP.m: Gradient of the function h_i(z) with respect \sigma.
rbfperiodicKernDisplay.m: Display parameters of the RBFPERIODIC kernel.
lfmvXlfmKernCompute.m: Velocity and position LFM kernel  
whiteblockKernDiagCompute.m: Compute diagonal of WHITEBLOCK kernel.
probitNoiseGradientParam.m: Gradient of PROBIT noise's parameters.
indexKernCompute.m: Compute the INDEX kernel given the parameters and X.
lfmGradientSigmaH3VP.m: Gradient of the function h_i(z) with respect \sigma.
kernSetWhite.m: Helper function to set the white noise in a kernel if it exists.
lvmLoadResult.m: Load a previously saved result.
dnetWriteResult.m: Write a DNET result.
gaussianNoiseParamInit.m: GAUSSIAN noise parameter initialisation.
lleDeconstruct.m: break LLE in pieces for saving.
diagKernDiagGradX.m: Gradient of DIAG kernel's diagonal with respect to X.
cmpndKernCompute.m: Compute the CMPND kernel given the parameters and X.
demSpgp1dGp1.m: Do a simple 1-D regression after Snelson & Ghahramani's example.
demSilhouettePlotTrue.m: Plot the true poses for the silhouette data.
probitNoisePointPlot.m: Plot the data-points for the PROBIT noise model.
modelExtractParam.m: Extract the parameters of a model.
mlpKernGradX.m: Gradient of MLP kernel with respect to input locations.
sdlfmvKernGradient.m: Gradient of SDLFM kernel's parameters.
lfmwhiteXlfmwhiteKernCompute.m: Compute a cross kernel between two LFM-WHITE
gpLogLikeGradients.m: Compute the gradients for the parameters and X.
rbfKernCompute.m: Compute the RBF kernel given the parameters and X.
modelGradientCheck.m: Check gradients of given model.
whitefixedKernCompute.m: Compute the WHITEFIXED kernel given the parameters and X.
rbfard2KernCompute.m: Compute the RBFARD kernel given the parameters and X.
lfmvXrbfKernGradient.m: Compute gradient between the LFMV and RBF kernels.
lfmwhiteComputeGradThetaH2.m: computes a portion of the LFM-WHITE kernel's gradient w.r.t. theta.
whitefixedKernExtractParam.m: Extract parameters from the WHITEFIXED kernel structure.
ggwhiteKernDisplay.m: Display parameters of the GG WHITE kernel.
mlpKernCompute.m: Compute the MLP kernel given the parameters and X.
pdinv.m: Invert a positive definite matrix.
xyzankurModify.m:  Helper function for modifying the point cloud from Agarwal and Triggs data.
demSwissRollFullLle5.m: Demonstrate LLE on the oil data.
multiKernCompute.m: Compute the MULTI kernel given the parameters and X.
vivmUspsResults.m: Summarise the USPS result files in LaTeX.
lfmComputeH3JV.m: Helper function for computing part of the LFMJV kernel.
simXsimKernGradient.m: Compute a cross gradient between two SIM kernels.
acclaimNumberOfFrames.m: Extract the number of frames.
sdlfmXsdlfmKernCompute.m: Compute a cross kernel between two SDLFM kernels.
invgammaPriorExtractParam.m: Extract params from inverse gamma prior structure.
gpTimeDynamicsCreate.m: Create the time dynamics model. 
gpReversibleDynamicsLatentGradients.m: Gradients of the X vector given the dynamics model.
heatKernCompute.m: Compute a kernel matrix for a HEAT kernel.
probitNoiseOut.m: Compute the output of the PROBIT noise given the input mean and variance.
lfmwhiteXlfmwhiteKernGradient.m: Compute a cross gradient between two
ngaussNoiseOut.m: Compute the output of the NGAUSS noise given the input mean and variance.
gibbsperiodicKernParamInit.m: GIBBSPERIODIC kernel parameter initialisation.
ngaussNoisePointPlot.m: Plot the data-points for the NGAUSS noise model.
lvmTwoDPlot.m: Helper function for plotting the labels in 2-D.
modelReadFromFile.m: Read model from a file FID produced by the C++ implementation.
mvuDeconstruct.m: break MVU in pieces for saving.
mvuOptions.m: Options for a MVU.
matern52KernDiagCompute.m: Compute diagonal of MATERN52 kernel.
xyzankurAnim.m: Animate point cloud of stick man from Agarwal & Triggs dataset.
xyzhumanevaVisualise2.m:
mlpLogLikeGradients.m: Multi-layer perceptron gradients.
logisticNormalPriorExpandParam.m: Expand logistic-normal prior structure from params.
cmpndNoiseNuG.m:  Update nu and g parameters associated with compound noise model.
mgaussianNoiseDisplay.m: Display parameters of the MGAUSSIAN noise.
lmcKernGradX.m: Gradient of LMC kernel with respect to input locations.
demVowelsFgplvm3.m: Model the vowels data with a 2-D FGPLVM using RBF kernel and back constraints, but without PCA initialisation.
ncnmNoiseLogLikelihood.m: Log-likelihood of data under null category noise model.
simXsimKernDiagGradient.m: Gradient for the diagonal between two SIM kernels.
preparePlot.m: Helper function for tidying up the plot before printing.
demRobotWirelessFgplvm3.m: Wireless Robot data from University of Washington with dynamics and no back constraints.
noneKernGradX.m: Gradient of NONE kernel with respect to a point x.
gaussianNoiseExpandParam.m: Create noise structure from GAUSSIAN noise's parameters.
xyzhumaneva2joint.m:
gaussianKernDiagGradient.m: Compute the gradient of the gaussian kernel's diagonal wrt parameters.
fgplvmOptimisePoint.m: Optimise the postion of a latent point.
kbrOut.m: Compute the output of a KBR model given the structure and input X.
skelReverseLookup.m: Return the number associated with the joint name.
gaussianKernGradient.m: Gradient of gaussian kernel's parameters.
lfmvXlfmvKernCompute.m: Velocity and velocity LFM kernel  
matern52KernDiagGradient.m: Compute the gradient of the MATERN52 kernel's diagonal wrt parameters.
rbfhKernParamInit.m: RBFH kernel parameter initialisation.
ivmEpLogLikelihood.m: Return the EP approximation to the log-likelihood.
lfmaXlfmvKernGradient.m: Compute a cross gradient between a LFMA and a LFMV.
demStickFgplvm3.m: Model the stick man using an RBF kernel and RBF kernel based back constraints.
modelWriteToFID.m: Write to a stream a given model.
mvuCreate.m: Maximum variance unfolding embedding model.
invgammaPriorGradient.m: Gradient wrt x of the log Gaussian prior.
gpDynamicsSamp.m: Sample from the dynamics for a given input.
sdlfmaXsdlfmvKernGradientBlock.m: Gradients of the parameters in block i,j
plotMatrix.m: Fill a given axis with a matrix plot.
lfmComputeH4JA.m: Helper function for computing part of the LFMJA kernel.
xyzVisualise.m: For drawing an xyz representation of 3-D data.
gpObjective.m: Wrapper function for GP objective.
gpDynamicsLogLikelihood.m: Give the log likelihood of GP dynamics.
ratquadKernExpandParam.m: Create kernel structure from RATQUAD kernel's parameters.
fgplvmObjective.m: Wrapper function for GP-LVM objective.
lfmaXrbfKernCompute.m: Compute cross kernel between the LFMA and RBF kernels.
datasetsDirectory.m: Returns directory where data is stored.
logisticNormalPriorLogProb.m: Log probability of logistic-normal prior.
lfmTest.m: Test the gradients of the LFM model.
multiKernDiagGradientBlock.m:
ggwhiteXggwhiteKernGradient.m: Compute a cross gradient between two GG WHITE kernels.
lvmSetPlotNoVar.m: A copy of lvmSetPlot where the variance in the input
lfmwhiteKernParamInit.m: LFM-WHITE kernel parameter initialisation.
nddisimKernParamInit.m: NDDISIM kernel parameter initialisation.
ivmOptions.m: Return default options for IVM model.
linearLogLikeGradients.m: Linear model gradients.
rbfXnoneKernGradient.m: Compute a cross gradient between RBF and DUMMY
rbfardKernDiagGradient.m: Compute the gradient of the RBFARD kernel's diagonal wrt parameters.
noiseDisplay.m: Display the parameters of the noise model.
fgplvmReconstruct.m: Reconstruct an FGPLVM from component parts.
multiKernCacheBlock.m:
heatKernDiagGradient.m: Gradient of the HEAT kernel's diagonal wrt parameters.
fgplvmSequenceObjective.m: Wrapper function for objective of a single sequence in latent space and the corresponding output sequence.
sdlfmXsdlfmKernGradientBlock.m: Gradients of the parameters in block i,j
mogLowerBound.m: Computes lower bound on log likelihood for an MOG model.
demCmu35gplvmFgplvm1.m: Learn a GPLVM on CMU 35 data set.
gpSubspaceCreate.m: 
disimXsimKernGradient.m: Compute gradient between the DISIM and SIM kernels.
gaussianwhiteKernGradX.m: Gradient of gaussian white kernel with respect 
tokenise.m: Split a string into separate tokens.
mlpKernExpandParam.m: Create kernel structure from MLP kernel's parameters.
sdlfmvXsdrbfKernComputeBlock.m: Cross kernel between SDLFM and SDRBF for i,j
polyKernParamInit.m: POLY kernel parameter initialisation.
lfmwhiteKernGradX.m: Gradient of LFM-WHITE kernel with respect to a point t.
polyardKernDisplay.m: Display parameters of the POLYARD kernel.
dnetTest.m: Test some settings for the density network.
sdlfmKernDiagComputeBlock.m: Diagonal of a SDLFM kernel matrix for block i
whiteKernGradX.m: Gradient of WHITE kernel with respect to input locations.
demThreeFiveIvm1.m: Try the IVM & NCNM on 3 vs 5.
priorExtractParam.m: Extract the prior model's parameters.
lfmvpGradientUpsilonMatrix.m: Gradient upsilon matrix vel. pos.
matern32KernGradient.m: Gradient of MATERN32 kernel's parameters.
orderedNoiseUpdateParams.m: Update parameters for ordered categorical noise model.
simwhiteXrbfinfwhiteKernCompute.m: Compute a cross kernel between a SIM-WHITE
ivmRemovePoint.m: Removes a given point from the IVM.
simwhiteXwhiteKernCompute.m: Compute a cross kernel between the SIM-WHITE
logisticNormalPriorParamInit.m: Logistic-normal prior model's parameter initialisation.
lfmGradientSigmaH4.m: Gradient of the function h_i(z) with respect \sigma.
cmpndKernDisplay.m: Display parameters of the CMPND kernel.
cmpndKernExtractParamTransformSettings.m: Extract parameter transform settings 
nddisimXndsimKernCompute.m: Compute a cross kernel between DISIM and SIM kernels with no decay in the SIM part.
gaussianKernGradX.m: Gradient of gaussian kernel with respect to input locations.
velotransKernParamInit.m: VELOTRANS kernel parameter initialisation.
noiseWriteParamsToFID.m: Write the noise parameters to a stream.
gpTimeDynamicsSetLatentValues.m: Set the latent values inside the model.
invcmpndKernCompute.m: Compute the INVERSE-PRECISION-CMPND kernel given the parameters and X.
noiseTest.m: Run some tests on the specified noise model.
optimiseParams.m: Optimise parameters.
sdlfmaXsdrbfKernGradientBlock.m: Gradients of the parameters in block i,j
simKernGradX.m: Gradient of SIM kernel with respect to each time point in t1.
wangPriorParamInit.m: Wang prior model's parameter initialisation.
heatKernExpandParam.m: Create kernel structure from HEAT kernel's parameters.
ivmPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
lfmvXrbfKernCompute.m: Compute a cross kernel between the LFMV and RBF kernels.
multimodelParamInit.m: MULTIMODEL model parameter initialisation.
modelLogLikelihood.m: Compute a model log likelihood.
gibbsKernExtractParam.m: Extract parameters from the GIBBS kernel structure.
whiteKernParamInit.m: WHITE kernel parameter initialisation.
mogLogLikelihood.m: Mixture of Gaussian's log likelihood.
lvmNearestNeighbour.m: Give the number of errors in latent space for 1 nearest neighbour.
sdlfmaXsdlfmaKernGradient.m: Gradients of cross kernel between 2 SDLFM kernels.
ngaussNoiseExtractParam.m: Extract parameters from the NGAUSS noise structure.
linearOptions.m: Options for learning a linear model.
noneKernDisplay.m: Display parameters of the NONE kernel.
readStringFromFID.m: Read an boolean from an FID.
cmpndNoiseExtractParam.m: Extract parameters from the CMPND noise structure.
multiKernFixBlocks.m:
kbrDisplay.m: Display parameters of the KBR model.
polyardKernExpandParam.m: Create kernel structure from POLYARD kernel's parameters.
modelExpandParam.m: Update a model structure with parameters.
xyzpoppeModify.m:
lfmGradientSigmaH3AP.m: Gradient of the function h_i(z) with respect \sigma.
uniformPriorExpandParam.m: Expand uniform prior structure from params.
ardKernGradX.m: Gradient of ARD kernel with respect to a point x.
ncnmNoiseGradientParam.m: Gradient of parameters for NCNM.
lvmClickVisualise.m: Visualise the manifold using clicks.
simComputeHStat.m: Helper function for computing part of the stationary version
ndsimKernGradient.m: Gradient of SIM kernel's parameters.
tensorKernDiagCompute.m: Compute diagonal of TENSOR kernel.
optimiDefaultOptimiser.m: Returns the default optimiser to be used.
polyardKernGradX.m: Gradient of POLYARD kernel with respect to input locations.
demBrendanFgplvm3.m: Use the GP-LVM to model the Frey face data with DTCVAR.
fgplvmDynamicsPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
demOilFgplvm4.m: Oil data with deterministic training conditional, and MLP back constraints.
fgplvmDynamicsRun.m: Runs auto regressive dynamics in a forward manner.
whiteKernCompute.m: Compute the white-noise (WHITE) kernel between
lfmvKernExpandParam.m: Create kernel structure from LFMV kernel's parameters.
lfmvpGradientUpsilonVector.m: Gradient upsilon vector vel. pos.
fgplvmPointLogLikelihood.m: Log-likelihood of a point for the GP-LVM.
imageModify.m: Helper code for visualisation of image data.
fgplvmSequenceLogLikeGradient.m: Log-likelihood gradient for of a sequence of the GP-LVM.
rbfinfwhiteKernDisplay.m: Display parameters of the RBF-WHITE kernel (with
simXsimKernCompute.m: Compute a cross kernel between two SIM kernels.
ivmReconstruct.m: Reconstruct an IVM form component parts.
cmpndKernExtractParam.m: Extract parameters from the CMPND kernel structure.
rbfhKernDiagCompute.m: Compute diagonal of RBFH kernel.
demUnlabelledOneIvm1.m: Test IVM code on a toy crescent data.
diagKernDiagCompute.m: Compute diagonal of DIAG kernel.
ivmPosteriorGradMeanVar.m: Gradient of mean and variances of the posterior wrt X.
sdlfmXsdlfmKernGradientICBlock.m: Partial derivatives initial conditions
demSilhouetteAverage.m: Shows the average of the poses.
lfmapGradientUpsilonMatrix.m: Gradient upsilon matrix accel. pos.
gradientCheck.m: Check gradients of objective function.
mocapLoadTextData.m: Load a motion capture data set from a text file.
doubleMatrixWriteToFID.m: Writes a double matrix to an FID.
heatXheatKernCompute.m: Compute a cross kernel between two HEAT kernels.
indexardKernExtractParam.m: Extract parameters from the INDEXARD kernel structure.
multiKernGradX.m: Gradient of MULTI kernel with respect to a point x.
lfmaXrbfvKernCompute.m: Compute cross kernel between the LFMA and RBFV kernels.
ppcaOptions.m: Options for probabilistic PCA.
invcmpndKernDiagGradient.m: Compute the gradient of the INVCMPND kernel's diagonal wrt parameters.
noneKernCompute.m: Compute the NONE kernel given the parameters and X.
gpDynamicsLatentGradients.m: Gradients of the X vector given the dynamics model.
modelLoadResult.m: Load a previously saved result.
getLocalGradAlphaOmega.m: Gradients of parameters in alpha and omega
demGpTwoSample.m: Test GP two sample code.
lfmvvComputeUpsilonMatrix.m: Upsilon matrix vel. vel. with t1, t2 limits
simwhiteXsimwhiteKernGradient.m: Compute a cross gradient between two
gpWriteResult.m: Write a GP result.
rbfhKernDisplay.m: Display parameters of the RBFH kernel.
gpDynamicsLogLikeGradients.m: Gradients of the GP dynamics wrt parameters.
ndsimKernDisplay.m: Display parameters of the NDSIM kernel.
gaussianwhiteKernGradient.m: Gradient of gaussian white kernel's parameters.
demSwissRollLle5.m: Demonstrate LLE on the oil data.
sdlfmXsdrbfKernGradientBlockIGJ.m: 
xyzhumanevaVisualiseModes.m:
demVowelsFgplvm2.m: Model the vowels data with a 2-D FGPLVM using RBF kernel.
lfmwhiteComputeH.m: Helper function for computing part of the LFM-WHITE
logisticNormalPriorSetBounds.m: Set logistic-normal prior bounds.
xyzhumanevaVisualise3d.m:
printLatexPlot.m: Print a plot to LaTeX.
fileKernGradX.m: Gradient of FILE kernel with respect to a point x.
sdlfmvXsdlfmKernGradientBlock.m: Gradients of the parameters in block i,j
ivmSelectPoint.m: Choose a point for inclusion or removal.
dexpKernGradX.m: Gradient of the double exponential kernel with respect to a
xyzhumanevaModify2.m:
invgammaPriorExpandParam.m: Expand inverse gamma prior structure from params.
lvmPrintPlot.m: Print latent space for learnt model.
lfmvvGradientUpsilonMatrix.m: Gradient upsilon matrix vel. vel.
normuniPriorExtractParam.m: Extract params from normal uniform prior structure.
ouKernDiagCompute.m: Compute diagonal of OU kernel (see ouKernCompute or
simwhiteKernDiagGradX.m: Gradient of SIM-WHITE kernel's diagonal w.r.t. t.
linKernGradient.m: Gradient of LIN kernel's parameters.
ndsimKernExpandParamTransformSettings.m: Create kernel structure from SIM kernel's parameters' transform settings.
ppcaPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
lfmComputeUpsilon.m: Helper function for comptuing part of the LFM kernel.
normuniPriorParamInit.m: Normal uniform prior model's parameter initialisation.
gaussianNoiseGradientParam.m: Gradient of GAUSSIAN noise's parameters.
ncnmNoiseDisplay.m: Display  parameters from null category noise model.
demSwissRollLle1.m: Demonstrate LLE on the oil data.
rocholFactorise.m: Rank one Cholesky factorise.
lfmaKernCompute.m: Compute the LFMA kernel given the parameters and X.
rbfardKernCompute.m: Compute the RBFARD kernel given the parameters and X.
xyzhumanevaModify.m:
lfmjpComputeUpsilonMatrix.m: Upsilon matrix jolt. pos. with t1, t2 limits
sdlfmvKernParamInit.m: SDLFMV kernel initialization
lnDiffCumGaussian.m: Log of the difference between two cumulative Gaussians.
orderedNoise3dPlot.m: Draws a 3D or contour plot for the ORDERED noise model.
gaussianNoiseLikelihood.m: Likelihood of the data under the GAUSSIAN noise model.
demOilFgplvm1.m: Oil data with fully independent training conditional.
rbfwhiteKernExtractParam.m: Extract parameters from the RBF-WHITE kernel
sdlfmvKernDiagCompute.m: Compute diagonal of a SDLFMV kernel.
noiseGradVals.m: Gradient of noise model wrt mu and varsigma.
ivmDowndateNuG.m: Downdate nu and g parameters associated with noise model.
disimComputeH.m: Helper function for comptuing part of the DISIM kernel.
matern32KernExtractParam.m: Extract parameters from the MATERN32 kernel structure.
whiteKernExpandParamTransformSettings.m: Create kernel structure from WHITE kernel's
gaussianwhiteKernDiagGradient.m: Compute the gradient of the gaussian white 
gpDynamicsPointLogLikelihood.m: Compute the log likelihood of a point under the GP dynamics model.
lfmKernGradX.m: Gradient of LFM kernel with respect to a point x.
noneKernExtractParam.m: Extract parameters from the NONE kernel structure.
lfmComputeH4JV.m: Helper function for computing part of the LFMJV kernel.
expKernCompute.m: Compute the EXP kernel given the parameters and X.
rbfperiodic2KernGradient.m: Gradient of RBFPERIODIC2 kernel's parameters.
whiteKernDiagCompute.m: Compute diagonal of WHITE kernel.
gaussianPriorParamInit.m: Gaussian prior model's parameter initialisation.
tensorKernGradient.m: Gradient of TENSOR kernel's parameters.
sdlfmaXsdlfmvKernGradient.m: Gradients of cross kernel between 2 SDLFM kernels.
sdlfmKernDiagCompute.m: Compute diagonal of a SDLFM kernel.
probitNoiseLogLikelihood.m: Log likelihood of the data under the PROBIT noise model.
probitNoiseExtractParam.m: Extract parameters from the PROBIT noise structure.
gpUpdateKernels.m: Update the kernels that are needed.
readBoolFromFID.m: Read a boolean from an FID.
mlpKernParamInit.m: MLP kernel parameter initialisation.
lleOptimise.m: Optimise an LLE model.
lfmaXlfmaKernGradient.m: Compute a cross gradient between a LFMA and a LFMA.
ggwhiteKernCompute.m: Compute the GG white kernel given the parameters and X.
gpTimeDynamicsExpandParam.m: Place the parameters vector into the model for GP time dynamics.
lfmapComputeUpsilonMatrix.m: Upsilon matrix acce. pos. with t1, t2 limits
noiseReadParamsFromFID.m: Read the noise parameters from C++ file FID.
lfmvpGradientSigmaUpsilonVector.m: Gradient of upsilon vector vp wrt sigma
sdlfmaXsdlfmKernGradientBlockILJ.m: 
demInterpolationGp.m: Demonstrate Gaussian processes for interpolation.
skelPlayData.m: Play skel motion capture data.
kernDiagGradX.m: Compute the gradient of the  kernel wrt X.
lfmGradientSigmaH3.m: Gradient of the function h_i(z) with respect \sigma.
gpReversibleDynamicsSamp.m: Sample from the dynamics for a given input.
kernExpandParam.m: Expand parameters to form a kernel structure.
bvhReadFile.m: Reads a bvh file into a tree structure.
skelConnectionMatrix.m: Compute the connection matrix for the structure.
fileKernParamInit.m: FILE kernel parameter initialisation.
vivmRunDataSetLearn.m: Try the virtual IVM on a data set and save the results.
modelLogLikeGradients.m: Compute a model's gradients wrt log likelihood.
lfmKernExtractParam.m: Extract parameters from the LFM kernel structure.
ivmSelectPoints.m: Selects the point for an IVM.
gpSubspaceExtractParam.m:
kernSetIndex.m: Set the indices on a compound kernel.
fgplvmLogLikelihood.m: Log-likelihood for a GP-LVM.
priorCreate.m: Create a prior structure given a type.
linardKernExpandParam.m: Create kernel structure from LINARD kernel's parameters.
scaleNoiseOut.m: A simple noise model that scales and centres the data.
fgplvmDisplay.m: Display an FGPLVM model.
modelSetLatentValues.m: Set the latent variables for dynamics models in the GPLVM.
priorSetBounds.m: Set the bounded prior model's bounds from bounds vector.
gpTimeDynamicsExtractParam.m: Extract parameters from the GP time dynamics model.
lfmKernParamInit.m: LFM kernel parameter initialisation. The latent force
gpTimeDynamicsLogLikelihood.m: Give the log likelihood of GP time dynamics.
rbfperiodic2KernGradX.m: Gradient of RBFPERIODIC2 kernel with respect to a point x.
rbfinfwhiteKernExpandParam.m: Create kernel structure from RBF-WHITE kernel's
tableRead.m: Read in data which has column titles in the first line and separated values in each other line.
rbfCreate.m: Wrapper for NETLAB's rbf `net'.
whitefixedKernExpandParam.m: Create kernel structure from WHITEFIXED kernel's parameters.
robThreeDynamicsExpandParam.m: Place the parameters vector into the model for first robot dynamics.
lfmaKernExpandParam.m: Create kernel structure from LFMA kernel's parameters.
lvmScatterPlotColor.m: 2-D scatter plot of the latent points with color.
linKernExpandParam.m: Create kernel structure from LIN kernel's parameters.
priorWriteParamsToFID.m: Write prior params from C++ written FID.
modelPointLogLikelihood.m: Compute the log likelihood of a given point.
sdlfmvXsdlfmvKernCompute.m: Compute a cross kernel between two SDLFMV kernels.
multimodelLogLikelihood.m: Log likelihood of MULTIMODEL model.
xyzhumanevaDraw.m:
simKernDiagGradX.m: Gradient of SIM kernel's diagonal with respect to the
cmpndKernParamInit.m: CMPND kernel parameter initialisation.
xyzhumanevaVisualise.m:
rbfKernDisplay.m: Display parameters of the RBF kernel.
matern52KernDiagGradX.m: Gradient of MATERN52 kernel's diagonal with respect to X.
sdlfmXsdlfmKernGradientBlockIEJ.m: 
matern32KernGradX.m: Gradient of MATERN32 kernel with respect to input locations.
sdlfmaMeanGradient.m: Gradients wrt parameters of the accel. mean SDLFM.
demCmu35TaylorNearestNeighbour.m: Recreate the Nearest Neighbour result from Taylor et al, NIPS 2006.
polyardKernParamInit.m: POLYARD kernel parameter initialisation.
xyzhumanevaAnim.m:
lfmaXlfmvKernCompute.m: Acceleration and velocity LFM kernel  
lvmScoreModel.m: Score model with a GP log likelihood.
scaleNoiseSites.m: Site updates for Scale noise model.
heatKernDiagCompute.m: Diagonal of a kernel matrix for a HEAT kernel.
isomapEmbed.m: Embed data set with Isomap.
xyzpoppeVisualise.m: Draw the Poppe figure return the graphics handle.
rbfhKernDiagGradX.m: Gradient of RBFH kernel's diagonal with respect to X.
polyardKernCompute.m: Compute the POLYARD kernel given the parameters and X.
indexKernDiagCompute.m: Compute diagonal of INDEX kernel.
demBrendanFgplvm4.m: Use the GP-LVM to model the Frey face data with DTCVAR and back constraints.
invCumGaussian.m: Computes inverse of the cumulative Gaussian.
ggwhiteKernGradient.m: Gradient of GG WHITE kernel's parameters.
kernExpandParamTransformSettings.m: Expand parameters' transform settings to form a kernel structure.
scatterPlot.m: 2-D scatter plot of labelled points.
lfmClassVisualise.m: Callback function to visualize LFM in 2D
sigmoidTransform.m: Constrains a parameter to be between 0 and 1.
matern52KernExtractParam.m: Extract parameters from the MATERN52 kernel structure.
invcmpndKernDiagCompute.m: Compute diagonal of INVCMPND kernel.
velotransKernDiagGradX.m: Gradient of VELOTRANS kernel's diagonal with respect to X.
ngaussNoiseGradVals.m: Gradient of NGAUSS noise log Z with respect to input mean and variance.
ivmGunnarData.m: Script for running experiments on Gunnar data.
ivmNegGradientNoise.m: Wrapper function for calling noise param gradients.
demThreeFiveResults.m: Plot results from the three vs five experiments.
gpTimeDynamicsDisplay.m: Display a GP time dynamics model.
rbfhKernExpandParam.m: Create kernel structure from RBFH kernel's parameters.
rbfinfwhiteXrbfinfwhiteKernCompute.m: Compute a cross kernel between two
whitehKernExpandParam.m: Create kernel structure from WHITEH kernel's parameters.
modelDisplay.m: Display a text output of a model.
gpPosteriorGradMeanVar.m: Gadient of the mean and variances of the posterior at points given by X.
ardKernDisplay.m: Display parameters of the ARD kernel.
sdlfmXsdlfmKernGradientIC.m: Computes partial derivative for init. const.
lfmComputeH3AP.m: Helper function for computing part of the LFMAP kernel.
whitefixedKernGradX.m: Gradient of WHITEFIXED kernel with respect to a point x.
fgplvmPointLogLikeGradient.m: Log-likelihood gradient for of a point of the GP-LVM.
lfmwhiteKernGradient.m: Gradient of LFM-WHITE kernel's parameters.
gpDynamicsExpandParam.m: Place the parameters vector into the model for GP dynamics.
ivmOptimiseIvm.m: Selects the points for an IVM model.
whitehKernDiagCompute.m: Compute diagonal of WHITEH kernel.
ggKernGradient.m: Gradient of GG kernel's parameters.
robThreeDynamicsLogLikelihood.m: Give the log likelihood of the robot three dynamics part.
demUspsIvm3.m: Try the ARD IVM on some digits data.
noiseReadFromFID.m: Load from an FID written by the C++ implementation.
gaussianNoiseLogLikelihood.m: Log likelihood of the data under the GAUSSIAN noise model.
robThreeDynamicsLatentGradients.m: Gradients of the X vector given the dynamics model.
ncnmNoiseParamInit.m: null category noise model's parameter initialisation.
disimKernCompute.m: Compute the DISIM kernel given the parameters and X.
rayleighSamp.m: Sample from a Rayleigh with a given sigma.
expKernDiagCompute.m: Compute diagonal of EXP kernel.
matern32KernDiagCompute.m: Compute diagonal of MATERN32 kernel.
ivmUspsResults.m: Summarise the USPS result files in LaTeX.
mlpLogLikelihood.m: Multi-layer perceptron log likelihood.
lleReconstruct.m: Reconstruct an LLE form component parts.
biasKernGradX.m: Gradient of BIAS kernel with respect to input locations.
rbfperiodicOut.m: Compute the output of a RBFPERIODIC model given the structure and input X.
invcmpndKernGradX.m: 
mvuEmbed.m: Embed data set with MVU.
matern32KernExpandParam.m: Create kernel structure from MATERN32 kernel's parameters.
demOptimiseGpTutorial.m: Shows that there is an optimum for the covariance function length scale.
lfmGradientSigmaH3VV.m: Gradient of the function h_i(z) with respect \sigma.
orderedNoiseGradientParam.m: Gradient of ORDERED noise's parameters.
gpTimeDynamicsLatentGradients.m: Gradients of the X vector given the time dynamics model.
scaleNoiseDisplay.m: Display the parameters of the scaled noise model.
rbfhKernDiagGradient.m: Gradient of the RBFH kernel's diagonal wrt parameters.
demUspsIvm1.m: Try the IVM on the USPS digits data with RBF kernel.
indexardKernParamInit.m: INDEXARD kernel parameter initialisation.
lvmSetPlot.m: Sets up the plot for visualization of the latent space.
sdlfmKernGradientMean.m: Gradients of the parameters of mean function cov.
kbrOptions.m: Create a default options structure for the KBR model.
translateKernDisplay.m: Display parameters of the TRANSLATE kernel.
invcmpndKernExtractParam.m: Extract parameters from the INVCMPND kernel structure.
sdlfmaXsdlfmaKernComputeBlock.m: Computes SDLFM kernel matrix for block i,j
sparseDiag.m: Create a diagonal matrix that is sparse from a vector.
isomapOptions.m: Options for a isomap.
kernReadParamsFromFID.m: Read the kernel parameters from C++ file FID.
lvmClassVisualise.m: Callback function for visualising data.
mapLoadData.m: Load a mapping model dataset (e.g. classification, regression).
ivmUpdateNuG.m: Update nu and g parameters associated with noise model.
fgplvmCovGradsTest.m: Test the gradients of the covariance.
tensorKernCompute.m: Compute the TENSOR kernel given the parameters and X.
indexardKernExpandParam.m: Create kernel structure from INDEXARD kernel's parameters.
demSilhouetteLinear1.m: Model silhouette data with independent linear models.
rbfperiodic2KernParamInit.m: RBFPERIODIC2 kernel parameter initialisation.
gradLogCumGaussian.m: Gradient of the log of the cumulative Gaussian.
xlsLoadData.m: Wrapper function for xlsread to get files from the datasets directory.
disimKernExpandParam.m: Create kernel structure from DISIM kernel's parameters.
mgaussianNoiseGradVals.m: Gradient of MGAUSSIAN noise log Z with respect to input mean and variance.
rbfhKernGradient.m: Gradient of RBFH kernel's parameters.
paramNameRegularExpressionLookup.m: Returns the indices of the parameter containing the given regular expression.
translateKernDiagGradX.m: Gradient of TRANSLATE kernel's diagonal with respect to X.
sdlfmvXsdlfmKernComputeBlock.m: Computes SDLFM kernel matrix for block i,j
fgplvmDynamicsFieldPlot.m: 2-D field plot of the dynamics.
bvhPlayFile.m: Play motion capture data from a bvh file.
disimXrbfKernGradient.m: Compute gradient between the DISIM and RBF kernels.
lfmapGradientSigmaUpsilonMatrix.m: Gradient of upsilon matrix ap wrt sigma
disimKernDiagCompute.m: Compute diagonal of DISIM kernel.
sdlfmvXsdlfmvKernGradientICBlock.m: Partial derivatives initial conditions
whiteXwhiteKernGradient.m: Compute a cross gradient between two WHITE kernels.
disimKernDiagGradX.m: Gradient of DISIM kernel's diagonal with respect to X.
rbfardKernGradX.m: Gradient of RBFARD kernel with respect to input locations.
modelOptions.m: Returns a default options structure for the given model.
rbfwhiteKernDisplay.m: Display parameters of the RBF-WHITE kernel.
kernTest.m: Run some tests on the specified kernel.
kernCreate.m: Initialise a kernel structure.
nddisimKernExtractParam.m: Extract parameters from the NDDISIM kernel structure.
gibbsperiodicKernDisplay.m: Display parameters of the GIBBSPERIODIC kernel.
linardKernCompute.m: Compute the LINARD kernel given the parameters and X.
ratquadKernGradient.m: Gradient of RATQUAD kernel's parameters.
ardKernDiagGradient.m: Compute the gradient of the ARD kernel's diagonal wrt parameters.
optimiDefaultConstraint.m: Returns function for parameter constraint.
diagKernDisplay.m: Display parameters of the DIAG kernel.
plot3Visualise.m:  Helper code for plotting a plot3 visualisation.
simKernParamInit.m: SIM kernel parameter initialisation.
rbfardjitKernParamInit.m: RBFARD2 kernel parameter initialisation.
ncnmNoiseGradVals.m: Compute gradient with respect to inputs to noise model.
translateKernGradient.m: Gradient of TRANSLATE kernel's parameters.
indexKernExpandParam.m: Create kernel structure from INDEX kernel's parameters.
lfmvKernDisplay.m: Display parameters of the LFMV kernel.
simComputeH.m: Helper function for comptuing part of the SIM kernel.
whiteXwhiteKernCompute.m: Compute a cross kernel between two WHITE kernels.
kernDiagGradient.m: Compute the gradient of the kernel's parameters for the diagonal.
gpMeanFunctionGradient.m: Compute the log likelihood gradient wrt the scales.
demOilFgplvm6.m: Oil data with partially independent training conditional, and MLP back constraints.
ouKernDisplay.m: Display parameters of the OU kernel (see ouKernCompute or
dnetLowerBound.m: Computes lower bound on log likelihood for an DNET model.
cell2num.m: Converts a cell array of strings to numbers.
readIntFromFID.m: Read an integer from an FID.
sdlfmaXsdrbfKernComputeBlock.m: Cross kernel between SDLFM and SDRBF for i,j
demSwissRollFullLle4.m: Demonstrate LLE on the oil data.
rbfard2KernDiagGradient.m: Compute the gradient of the RBFARD2 kernel's diagonal wrt parameters.
rbfKernGradient.m: Gradient of RBF kernel's parameters.
lfmwhiteXrbfwhiteKernCompute.m: Compute a cross kernel between an LFM-WHITE
lfmavGradientSigmaUpsilonMatrix.m: Gradient of upsilon matrix av wrt sigma
gpPosteriorVar.m: Variances of the posterior at points given by X.
rotationMatrix.m: Compute the rotation matrix for an angle in each direction.
sheatKernDiagCompute.m: Compute a diagonal for the SHEAT kernel matrix.
lfmComputeH4.m: Helper function for computing part of the LFM kernel.
uniformPriorLogProb.m: Log probability of uniform prior.
gaussianKernDiagGradX.m: Gradient of gaussian kernel's diagonal with respect to X.
sdlfmaXsdlfmvKernComputeBlock.m: Computes SDLFM kernel matrix for block i,j
robOneDynamicsExpandParam.m: Place the parameters vector into the model for first robot dynamics.
ncnmNoisePointPlot.m: Plot the data-points for null category noise model.
translateKernDiagCompute.m: Compute diagonal of TRANSLATE kernel.
multimodelExpandParam.m: Create model structure from MULTIMODEL model's parameters.
mogUpdateCovariance.m: Update the covariances of an MOG model.
computeKernel.m: Compute the kernel given the parameters and X.
diagKernCompute.m: Compute the DIAG kernel given the parameters and X.
lfmvpGradientSigmaUpsilonMatrix.m: Gradient of upsilon matrix vp wrt sigma
mapmodelReadFromFID.m: Load from a FID produced by C++ code.
cmpndNoise3dPlot.m: Draws a 3D or contour plot for the CMPND noise model.
fgplvmFieldPlot.m: 2-D field plot of the dynamics.
rbfperiodic2KernDiagCompute.m: Compute diagonal of RBFPERIODIC2 kernel.
simKernDiagGradient.m: Compute the gradient of the SIM kernel's diagonal wrt parameters.
sdrbfKernGradient.m: Gradient of SDRBF kernel's parameters.
xyzankur2joint.m: Converts data to xyz positions for each joint.
rbfardjitKernDisplay.m: Display parameters of the RBFARDJIT kernel.
sqexpKernGradient.m: Gradient of SQEXP kernel's parameters.
mgaussianNoiseExtractParam.m: Extract parameters from the MGAUSSIAN noise structure.
negNoiseLogLikelihood.m: Wrapper function for calling noise likelihoods.
demBrendanPpca1.m: Use PPCA to model the Frey face data with five latent dimensions.
gpPosteriorSample.m: Create a plot of samples from a posterior covariance.
cmpndNoiseLogLikelihood.m: Log likelihood of the data under the CMPND noise model.
printLatexText.m: Print a text string to file for latex input.
lfmvvGradientSigmaUpsilonMatrix.m: Gradient of upsilon matrix vv wrt sigma
sdlfmXsdrbfKernGradientBlock.m: Gradients of the parameters in block i,j
mlpKernDiagGradX.m: Gradient of MLP kernel's diagonal with respect to X.
rbfard2KernGradient.m: Gradient of RBFARD2 kernel's parameters.
whiteblockKernDisplay.m: Display parameters of the WHITEBLOCK kernel.
wienerKernDiagCompute.m: Compute diagonal of WIENER kernel.
simwhiteKernExtractParam.m: Extract parameters from the SIM-WHITE kernel
ouKernParamInit.m: Ornstein-Uhlenbeck (OU) kernel parameter initialisation.
whitefixedKernDiagGradX.m: Gradient of WHITEFIXED kernel's diagonal with respect to X.
safeSave.m: Safe save
simXsimComputeDiagHStat.m: Helper function for computing part of the stationary version
fileKernDiagCompute.m: Compute diagonal of FILE kernel.
demVowelsLle.m: Model the vowels data with a 2-D FGPLVM using RBF kernel.
diagKernExtractParam.m: Extract parameters from the DIAG kernel structure.
noneKernDiagCompute.m: Compute diagonal of NONE kernel.
modelObjective.m: Objective function to minimise for given model.
gpDeconstruct.m: break GP in pieces for saving.
gaussianNoiseDisplay.m: Display parameters of the GAUSSIAN noise.
laplacePriorParamInit.m: Laplace prior model's parameter initialisation.
fileKernExpandParam.m: Create kernel structure from FILE kernel's parameters.
sdrbfKernParamInit.m: SDRBF kernel initialization
ivmCovarianceGradient.m: The gradient of the likelihood approximation wrt the covariance.
invgammaPriorLogProb.m: Log probability of inverse gamma prior.
gpExpandParam.m: Expand a parameter vector into a GP model.
linearOptimise.m: Optimise a linear model.
probitNoise3dPlot.m: Draws a 3D or contour plot for the PROBIT noise model.
ndsimKernDiagCompute.m: Compute diagonal of NDSIM kernel.
priorExpandParam.m: Expand the prior model's parameters from params vector.
polyardKernGradient.m: Gradient of POLYARD kernel's parameters.
rbfperiodicExpandParam.m: Create model structure from RBFPERIODIC model's parameters.
robOneDynamicsLogLikeGradients.m: Gradients of the robot one dynamics wrt parameters.
lfmwhiteComputeGradThetaH1.m: computes a portion of the LFM-WHITE kernel's gradient w.r.t. theta.
kernPca.m: performs KPCA.
polyKernDiagGradX.m: Gradient of POLY kernel's diagonal with respect to X.
rbfinfwhiteKernParamInit.m: The RBF-WHITE-INF kernel is a convolutional
ivmReadFromFile.m: Load a file produced by the C++ implementation.
ouKernExtractParam.m: Extract parameters from the OU kernel structure (see
whitehKernGradX.m: Gradient of WHITEH kernel with respect to input locations.
ncnmNoiseSites.m: Site updates for null category model.
robThreeDynamicsExtractParam.m: Extract parameters from the robot three dynamics model.
gammaPriorExpandParam.m: Expand gamma prior structure from params.
gaussianwhiteKernExtractParam.m: Extract parameters from the gaussian white 
sdlfmvXsdlfmvKernComputeBlock.m: Computes SDLFM kernel matrix for block i,j
gibbsKernDiagGradient.m: Compute the gradient of the GIBBS kernel's diagonal wrt parameters.
simXsimComputeDiagH.m: Helper function for comptuing part of the SIM kernel.
ggwhiteKernParamInit.m: GG WHITE kernel parameter initialisation.
rbfOutputGradX.m: Evaluate derivatives of a RBF model's output with respect to inputs.
lvmScatterPlot.m: 2-D scatter plot of the latent points.
lfmComputeH.m: Helper function for computing part of the LFM kernel.
multimodelOptions.m: Create a default options structure for the MULTIMODEL model.
lvmResultsClick.m: Load a results file and visualise them with clicks
dexpKernDiagCompute.m: Compute diagonal of the double exponential kernel.
rbfperiodicKernDiagCompute.m: Compute diagonal of RBFPERIODIC kernel.
paramNameReverseLookup.m: Returns the index of the parameter with the given name.
expKernExpandParam.m: Create kernel structure from EXP kernel's parameters.
gpLogLikelihood.m: Compute the log likelihood of a GP.
linKernDiagCompute.m: Compute diagonal of LIN kernel.
noiseGradX.m: Returns the gradient of the log-likelihood wrt x.
whiteXwhiteKernGradX.m:
gpTimeDynamicsSequenceLogLikelihood.m: Return the log likelihood of a given latent sequence.
rbfperiodicKernCompute.m: Compute the RBFPERIODIC kernel given the parameters and X.
mlpExpandParam.m: Update mlp model with new vector of parameters.
demBrendanFgplvm5.m: Use the GP-LVM to model the Frey face data with DTCVAR and five latent dimensions..
rbfwhiteXrbfwhiteKernGradient.m: Compute a cross gradient between two
sdlfmaXsdlfmKernGradient.m: Gradients of cross kernel between 2 SDLFM kernels.
plot3Modify.m: Helper code for visualisation of 3-d data.
whitehKernCompute.m: Compute the WHITEH kernel given the parameters and X.
kbrExtractParam.m: Extract parameters from the KBR model structure.
indexKernExtractParam.m: Extract parameters from the INDEX kernel structure.
rbfwhiteXwhiteKernCompute.m: Compute a cross kernel between the RBF-WHITE
ndsimKernParamInit.m: SIM kernel parameter initialisation.
dnetUpdateBeta.m: Do an M-step (update parameters) on an Density Network model.
componentKernWriteParamsToFID.m: Write a component based kernel to a stream.
sdlfmXsdlfmvKernGradientICBlock.m: Partial derivatives initial conditions
sdrbfKernDisplay.m: Display parameters of the SDRBF kernel.
fgplvmReadFromFile.m: Load a file produced by the C++ implementation.
gaussianwhiteKernCompute.m: Compute the covariance of the output samples 
indexardKernCompute.m: Compute the INDEXARD kernel given the parameters and X.
biasKernExtractParam.m: Extract parameters from the BIAS kernel structure.
lfmExpandParam.m: Expand the given parameters into a LFM structure.
isomapDeconstruct.m: break isomap in pieces for saving.
whitefixedXwhitefixedKernGradient.m: Compute a cross gradient between two WHITEFIXED kernels.
rbfardjitKernExtractParam.m: Extract parameters from the RBFARD2 kernel structure.
rbfperiodic2KernDiagGradX.m: Gradient of RBFPERIODIC2 kernel's diagonal with respect to X.
cmpndNoiseDisplay.m: Display parameters of the CMPND noise.
sigmoidabTransform.m: Constrains a parameter to be between A and B
noiseGradientParam.m: Gradient wrt the noise model's parameters.
gpPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
mogOptions.m: Sets the default options structure for MOG models.
timeseriesdata.m: make a time series data set with the given window length.
gpDynamicsDisplay.m: Display a GP dynamics model.
ivmDowndateSites.m: Downdate site parameters.
mgaussianNoisePointPlot.m: Plot the data-points for the MGAUSSIAN noise model.
dnetOutputGrad.m: Evaluate derivatives of dnet model outputs with respect to parameters.
demSpgp1dGp5.m: Do a simple 1-D regression after Snelson & Ghahramani's example.
matern52KernParamInit.m: MATERN52 kernel parameter initialisation.
findAcyclicNeighbours2.m: find the k nearest neighbours for each point in Y preventing cycles in the graph.
treeFindRoots.m: Return indices of all root nodes in a tree structure.
mlpardKernGradX.m: Gradient of MLPARD kernel with respect to input locations.
ivmKernelObjective.m: Compute the negative of the IVM log likelihood approximation.
lfmKernGradient.m: Gradient of LFM kernel's parameters.
lfmvpComputeUpsilonDiagVector.m: Upsilon diag vector vel. pos. with t1, t2 limits
cumGamma.m: Cumulative distribution for gamma.
writeVersionToFID.m: Writes a version to an FID.
mocapConnections.m: Give a connection matrix for the motion capture data.
matern52KernDisplay.m: Display parameters of the MATERN52 kernel.
vectorVisualise.m:  Helper code for plotting a vector during 2-D visualisation.
disimKernParamInit.m: DISIM kernel parameter initialisation.
dnetPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
sdlfmaXsdlfmaKernCompute.m: Cross kernel between a SDLFMA and a SDLFMA kernels.
modelOutputGrad.m: Compute derivatives with respect to params of model outputs.
isomapOptimise.m: Optimise an ISOMAP model.
heatXrbfhKernCompute.m: Cross kernel between a HEAT and a RBF kernels.
gpTwoSample.m: Do Oliver Stegles simple two sample test on a data set.
dexpKernParamInit.m: The double exponential kernel is usually called
gpDataIndices.m: Return indices of present data.
whiteblockKernGradient.m: Gradient of WHITEBLOCK kernel's parameters.
ncnmNoiseOut.m: Ouput from null category noise model.
fgplvmWriteResult.m: Write a FGPLVM result.
fgplvmGradient.m: GP-LVM gradient wrapper.
treeSwapNode.m: Swap two nodes in the tree structure array.
priorReadFromFID.m: Read a prior from a C++ written FID.
lfmGradientSigmaH4VP.m: Gradient of the function h_i(z) with respect \sigma.
gpSubspaceOptimise.m:
matern52KernCompute.m: Compute the MATERN52 kernel given the parameters and X.
simXrbfKernGradient.m: Compute gradient between the SIM and RBF kernels.
indexardKernDisplay.m: Display parameters of the INDEXARD kernel.
optOptions.m: Give optimisation options for NETLAB.
dnetCreate.m: Density network model.
fileKernCompute.m: Compute the FILE kernel given the parameters and X.
priorGradient.m: Gradient of the prior with respect to its variables
linardKernGradX.m: Gradient of LINARD kernel with respect to input locations.
cmpndNoiseParamInit.m: CMPND noise parameter initialisation.
heatKernExtractParam.m: Extract parameters from the HEAT kernel structure.
simwhiteXsimwhiteKernCompute.m: Compute a cross kernel between two SIM-WHITE
probitNoiseParamInit.m: PROBIT noise parameter initialisation.
readVersionFromFID.m: Read version number from an FID.
modelParamInit.m: Initialise the parameters of the model.
jitChol.m: Do a Cholesky decomposition with jitter.
lfmwhiteKernExtractParam.m: Extract parameters from the LFM-WHITE kernel
ivmApproxLogLikelihood.m: Return the approximate log-likelihood for the IVM.
mlpOutputGradX.m: Evaluate derivatives of mlp model outputs with respect to inputs.
lmvuEmbed.m: Embed data set with landmark MVU
ppcaCreate.m: Density network model.
xyzankurAnimCompare.m: Animate a prediction and ground truth for stick man from Agarwal & Triggs dataset.
ggwhiteXwhiteKernCompute.m: Compute a cross kernel between a GG white and
sdlfmvXsdrbfKernCompute.m: Cross kernel between a SDLFMV and a SDRBF kernels.
kbrCreate.m: Create a KBR model.
gaussOverDiffCumGaussian.m: A Gaussian over difference of cumulative Gaussians.
robTwoDynamicsSetLatentValues.m: Set the latent values inside the model.
gibbsperiodicKernDiagCompute.m: Compute diagonal of GIBBSPERIODIC kernel.
demGpTwoSampleLifsh.m: Run GP two sample code on LifSh.
negNoiseGradientParam.m: Wrapper function for calling noise gradients.
rbfwhiteXrbfwhiteKernCompute.m: Compute a cross kernel between two RBF-WHITE
rbfardjitKernGradient.m: Gradient of RBFARD2 kernel's parameters.
cmpndKernExpandParam.m: Create kernel structure from CMPND kernel's parameters.
htkLoadMmf.m: File for loading synthesis data from HTK files.
expKernExtractParam.m: Extract parameters from the EXP kernel structure.
lfmComputeH4VV.m: Helper function for computing part of the LFMVXLFMV kernel.
translateKernCompute.m: Compute the TRANSLATE kernel given the parameters and X.
matern52KernExpandParam.m: Create kernel structure from MATERN52 kernel's parameters.
matern52KernGradient.m: Gradient of MATERN52 kernel's parameters.
sdlfmvXsdlfmvKernGradientBlock.m: Gradients of the parameters in block i,j
ardKernExpandParam.m: Create kernel structure from ARD kernel's parameters.
modelHessian.m: Hessian of error function to minimise for given model.
heatXheatKernGradient.m: Gradient wrt parameters between two HEAT kernels.
lnCumGaussian.m: log cumulative distribution for the normalised Gaussian.
velotransKernDiagCompute.m: Compute diagonal of VELOTRANS kernel.
robTwoDynamicsDisplay.m: Display the robot dynamics model. 
rbfOptimise.m: Optimise RBF for given inputs and outputs.
gradFuncWrapper.m: Wrapper function to enable use of Carl Rasmussen's minimze function.
noiseUpdateSites.m: Update site parameters for a given noise model.
biasKernDiagCompute.m: Compute diagonal of BIAS kernel.
ggwhiteXgaussianwhiteKernCompute.m: Compute a cross kernel between the GG white and GAUSSIAN white kernels.
whiteblockKernCompute.m: Compute the WHITEBLOCK kernel. 
lfmGradientH32.m: Gradient of the function h_i(z) with respect to some of the
fgplvmOptions.m: Return default options for FGPLVM model.
tensorKernSetIndex.m: Set the indices in the tensor kernel.
rocholhFactorise.m: Rank one Cholesky factorise.
velotransKernGradient.m: Gradient of VELOTRANS kernel's parameters.
sdlfmvKernCompute.m: Compute the SDLFMV kernel given the parameters and X.
orderedNoiseExpandParam.m: Create noise structure from ORDERED noise's parameters.
cgcarl.m: Wrapper for Carl Rasmussen's conjugate gradient implemntation.
velotransKernGradX.m: Gradient of VELOTRANS kernel with respect to a point x.
sdlfmXsdlfmvKernGradientBlock.m: Gradients of the parameters in block i,j
multiKernGradientBlock.m:
whitefixedKernDisplay.m: Display parameters of the WHITEFIXED kernel.
rbfinfwhiteKernDiagGradX.m: Gradient of RBF-WHITE kernel's (with integration
linKernGradX.m: Gradient of LIN kernel with respect to input locations.
noneKernGradient.m: Gradient of NONE kernel's parameters.
gplvmCmu35Animate.m: Animate the test data jointly with predictions.
lfmvvComputeUpsilonDiagVector.m: Upsilon vector vel. vel. with t1 = t2
rbfinfwhiteKernDiagCompute.m: Compute diagonal of RBF-WHITE kernel (with
orderedNoiseDisplay.m: Display parameters of the ORDERED noise.
wienerKernGradX.m: Gradient of WIENER kernel with respect to a point x.
lfmComputeTest.m: Test the file lfmComputeH.
dnetExpandParam.m: Update dnet model with new vector of parameters.
lfmvpComputeUpsilonVector.m: Upsilon vector for vel. pos. with t1 limit
probitNoiseDisplay.m: Display parameters of the PROBIT noise.
uniformPriorExtractParam.m: Extract params from uniform prior structure.
negLogLogitTransform.m: Constrains a parameter to be positive.
gpTest.m: Test the gradients of the gpCovGrads function and the gp models.
pskernelGradient.m: Gradient on likelihood approximation for point set IVM.
mlpOut.m: Output of an MLP model.
gibbsperiodicKernGradX.m: Gradient of GIBBSPERIODIC kernel with respect to a point x.
ppcaDeconstruct.m: break PPCA in pieces for saving.
gaussianNoise3dPlot.m: Draws a 3D or contour plot for the GAUSSIAN noise model.
sdlfmKernGradient.m: Gradient of SDLFM kernel's parameters.
linearExtractParam.m: Extract weights from a linear model.
lfmaKernDisplay.m: Display parameters of the LFMA kernel.
lvmScatterPlotNeighbours.m: 2-D scatter plot of the latent points with neighbourhood.
ardKernDiagGradX.m: Gradient of ARD kernel's diagonal with respect to X.
sdlfmKernParamInit.m: SDLFM kernel initialization
whitehKernDiagGradient.m: Compute the gradient of the WHITEH kernel's diagonal wrt parameters.
simKernGradient.m: Gradient of SIM kernel's parameters.
demUspsIvm2.m: Try the IVM on the USPS digits data with MLP kernel.
ivmOut.m: Evaluate the output of an IVM model.
acclaimReadSkel.m: Reads an ASF file into a skeleton structure.
demGpTwoSampleLif.m: Run GP two sample code on LIF.
mlpKernExtractParam.m: Extract parameters from the MLP kernel structure.
uniformPriorParamInit.m: Uniform prior model's parameter initialisation.
rbfOptions.m: Default options for RBF network.
linearParamInit.m: Initialise the parameters of an LINEAR model.
cmpndKernDiagGradient.m: Compute the gradient of the CMPND kernel's diagonal wrt parameters.
demStickFgplvm5.m: Model the stick man using an RBF kernel and regressive dynamics.
robOneDynamicsDisplay.m: Display the robot dynamics model. 
rbfOut.m: Output of an RBF model.
fgplvmViterbiSequence.m: Viterbi align a latent sequence.
gaussianwhiteKernDiagGradX.m: Gradient of gaussian white kernel's diagonal with respect to X.
sdlfmaXsdlfmaKernGradientBlockIGJ.m: 
lfmXlfmKernCompute.m: Compute a cross kernel between two LFM kernels.
polyKernGradX.m: Gradient of POLY kernel with respect to input locations.
robThreeDynamicsLogLikeGradients.m: Gradients of the robot three dynamics wrt parameters.
cmpndKernDiagCompute.m: Compute diagonal of CMPND kernel.
ggXggKernCompute.m: Compute a cross kernel between two GG kernels.
dexpKernExtractParam.m: Extract parameters from the double exponential's
gpReversibleDynamicsSetLatentValues.m: Set the latent values inside the model.
simXrbfKernCompute.m: Compute a cross kernel between the SIM and RBF kernels.
ivmDisplay.m: Display parameters of an IVM model.
lfmGradientSigmaH4VV.m: Gradient of the function h_i(z) with respect \sigma.
mogSample.m: Sample from a mixture of Gaussians model.
kbrOptimise.m: Optimise a KBR model.
gpWriteToFile.m: Write a file to be read by the C++ implementation.
ndsimKernExtractParam.m: Extract parameters from the SIM kernel structure.
findNeighbours.m: find the k nearest neighbours for each point in Y.
rbfperiodicCreate.m: Create a RBFPERIODIC model.
ivmReadFromFID.m: Load from a FID produced by the C++ implementation.
lfmwhiteKernDisplay.m: Display parameters of the LFM-WHITE kernel.
sdlfmvXsdlfmKernGradient.m: Gradients of cross kernel between 2 SDLFM kernels.
noiseLogLikelihood.m: Return the log-likelihood under the noise model.
velotransKernCompute.m: Compute the VELOTRANS kernel given the parameters and X.
sdlfmvXsdlfmKernGradientICBlock.m: Partial derivatives initial conditions
nddisimXndsimKernGradient.m: Compute a cross gradient between NDDISIM and NDSIM
springDampersVisualise.m: Helper code for showing an spring dampers during 2-D visualisation.
robTwoDynamicsExtractParam.m: Extract parameters from the robot two dynamics model.
demOrderedOneIvm1.m: Run a demonstration of the ordered categories noise model (linear data).
ggXggKernGradient.m: Compute a cross gradient between two GG kernels.
gpPosteriorMeanCovarTest.m: Test the gradients of the mean and covariance.
lfmKernCompute.m: Compute the LFM kernel given the parameters and X.
demGpCov2D.m: Simple demonstration of sampling from a covariance function.
sigmoid.m: The sigmoid function
mgaussianNoiseLogLikelihood.m: Log likelihood of the data under the MGAUSSIAN noise model.
priorReadParamsFromFID.m: Read prior params from C++ written FID.
xyzankurVisualise.m: Draw the Agarwal & Triggs figure return the graphics handle.
ggKernExpandParam.m: Create kernel structure from GG kernel's parameters.
demGpTwoSampleEB.m: Run GP two sample code on EB.
gpSample.m: Create a plot of samples from a GP.
gpTimeDynamicsSequenceLogLikeGradient.m: Log-likelihood gradient for of a sequence of the GP-LVM time dynamics.
gibbsperiodicKernExpandParam.m: Create kernel structure from GIBBSPERIODIC kernel's parameters.
skelModify.m: Update visualisation of skeleton data.
lfmUpdateKernels.m: Updates the kernel representations in the LFM structure.
invgammaPriorParamInit.m: Inverse gamma prior model's parameter initialisation.
ngaussNoiseSites.m: Site updates for noiseless Gaussian noise model.
gaussianNoiseExtractParam.m: Extract parameters from the GAUSSIAN noise structure.
gpBlockIndices.m: Return indices of given block.
mappingOptimise.m: Optimise the given model.
linKernCompute.m: Compute the LIN kernel given the parameters and X.
sdrbfKernExtractParam.m: Extract parameters from the SDRBF kernel structure.
mlpardKernDiagCompute.m: Compute diagonal of MLPARD kernel.
dnetLoadResult.m: Load a previously saved result.
dnetObjective.m: Wrapper function for Density Network objective.
biasKernExpandParam.m: Create kernel structure from BIAS kernel's parameters.
parseNobleKernelFile.m: Parse a kernel file from Bill Stafford Noble.
lfmGradientH42VP.m: Gradient of the function h_i(z) with respect to some of the
demGpCovFuncSample.m: Sample from some different covariance functions.
whiteXnoneKernGradient.m: Compute a cross gradient between WHITE and DUMMY kernels.
demCmu35gplvmFgplvm3.m: Learn a GPLVM on CMU 35 data set.
rbfperiodicParamInit.m: RBFPERIODIC model parameter initialisation.
xyzankurDraw.m: Helper function for drawing the point cloud from Agarwal and Triggs data.
wienerKernExtractParam.m: Extract parameters from the WIENER kernel structure.
isoctave.m: Returns true if the software running is Octave.
lfmComputeH3AV.m: Helper function for computing part of the LFMAV kernel.
sdlfmKernDisplay.m: Display parameters of the SDLFM kernel.
lfmjpComputeUpsilonVector.m: Upsilon vector jolt. pos. with t1, t2 limits
diagKernGradX.m: Gradient of DIAG kernel with respect to a point x.
invcmpndKernParamInit.m: INV.PRECISION-CMPND kernel parameter initialisation.
orderedNoiseLogLikelihood.m: Log likelihood of the data under the ORDERED noise model.
rbfhKernCompute.m: Compute the RBFH kernel given the parameters and X.
gibbsKernCompute.m: Compute the GIBBS kernel given the parameters and X.
dnetOut.m: Output of an DNET model.
rbfardjitKernDiagGradX.m: Gradient of RBFARDJIT kernel's diagonal with respect to X.
demUnlabelledOneIvm2.m: Test IVM code on a toy crescent data.
velotransKernExtractParam.m: Extract parameters from the VELOTRANS kernel structure.
dnetGradient.m: Density Network gradient wrapper.
spectrumVisualise.m: Helper code for showing an spectrum during 2-D visualisation.
demOilLle2.m: Demonstrate LLE on the oil data.
dnetLogLikeGradients.m: Density network gradients.
fgplvmClassVisualise.m: Callback function for visualising data in 2-D.
sdlfmvKernExtractParam.m: Extract parameters from the SDLFMV kernel structure.
gibbsKernSetLengthScaleFunc.m: Set the length scale function of the GIBBS kernel.
mogTwoDPlot.m: Helper function for plotting the labels in 2-D.
sdlfmaXsdlfmvKernCompute.m: Cross kernel between a SDLFMA and a SDLFMV kernels.
indexKernDisplay.m: Display parameters of the INDEX kernel.
lfmComputeH3JP.m: Helper function for computing part of the LFMJP kernel.
ppcaReconstruct.m: Reconstruct an PPCA form component parts.
orderedNoiseLikelihood.m: Likelihood of the data under the ORDERED noise model.
whitefixedKernDiagGradient.m: Compute the gradient of the WHITEFIXED kernel's diagonal wrt parameters.
expKernDisplay.m: Display parameters of the EXP kernel.
sdlfmXsdlfmKernGradientBlockILJ.m: 
stringSigFigs.m: Convert number to a string with a number of significant digits.
demSpgp1dGp4.m: Do a simple 1-D regression after Snelson & Ghahramani's example.
tensorKernDisplay.m: Display parameters of the TENSOR kernel.
acclaimGradient.m: computes the gradient of x,y,z locations wrt angles.
lfmaaComputeUpsilonDiagVector.m: Diag. of Upsilon matrix acce. accel. 
rbfard2KernParamInit.m: RBFARD2 kernel parameter initialisation.
xyzpoppe2joint.m:
polyKernGradient.m: Gradient of POLY kernel's parameters.
ivmInit.m: Initialise the IVM model.
robOneDynamicsLogLikelihood.m: Give the log likelihood of the robot one dynamics part.
ivmOptimise.m: Optimise the IVM.
gpPosteriorMeanCovar.m: Mean and covariances of the posterior at points given by X.
ggKernDisplay.m: Display parameters of the GG kernel.
lfmResultsDynamicWalking.m: Load a results file and visualise them.
gaussSamp.m: Sample from a Gaussian with a given covariance.
demOrderedTwoIvm1.m: Run a demonstration of the ordered categories noise model (circular data).
kernParamInit.m: Kernel parameter initialisation.
linard2KernDiagGradX.m: Gradient of LINARD2 kernel's diagonal with respect to X.
gaussianKernCompute.m: Compute the Gaussian kernel given the parameters and X.
ivmApproxLogLikeKernGrad.m: Gradient of the approximate likelihood wrt kernel parameters.
robOneDynamicsLatentGradients.m: Gradients of the X vector given the dynamics model.
kpcaEmbed.m: Embed data set with kernel PCA.
ndsimKernCompute.m: Compute the NDSIM kernel with no decay given the parameters and X.
rbfwhiteKernGradient.m: Gradient of RBF-WHITE kernel's parameters.
lfmOptions.m: Creates a set of default options for a LFM model.
lfmwhiteComputePsi.m: Helper function for comptuing part of the LFM-WHITE
simSample.m: Sample from SIM kernel
linardKernExtractParam.m: Extract parameters from the LINARD kernel structure.
whiteKernGradient.m: Gradient of white-noise (WHITE) kernel's
lfmaaGradientUpsilonMatrix.m: Gradient upsilon matrix accel. accel.
gaussianwhiteKernDiagCompute.m: Compute diagonal of gaussian white kernel.
uniformPriorSetBounds.m: Set uniform prior bounds.
sqexpKernCompute.m: Compute the SQEXP kernel given the parameters and X.
tensorKernExpandParam.m: Create kernel structure from TENSOR kernel's parameters.
simwhiteKernGradient.m: Gradient of SIM-WHITE kernel's parameters.
invcmpndKernDiagGradX.m: Gradient of INVCMPND kernel's diagonal with respect to X.
ivmCreate.m: Create a IVM model with the IVM sparse approximaiton.
kernReadFromFID.m: Load from an FID written by the C++ implementation.
biasKernDiagGradX.m: Gradient of BIAS kernel's diagonal with respect to X.
gpSequenceLogLikeGradient.m: Log-likelihood gradient for of a sequence of the GP-LVM.
fgplvmKernDynamicsSample.m: Sample a field from a given kernel.
linard2KernGradient.m: Gradient of LINARD2 kernel's parameters.
wangPriorExpandParam.m: Expand wang prior structure from params.
polyardKernDiagGradX.m: Gradient of POLYARD kernel's diagonal with respect to X.
writeDoubleToFID.m: Writes a double to an FID.
kernCompute.m: Compute the kernel given the parameters and X.
gpReversibleDynamicsOptions.m: Return default options for GP reversible dynamics model.
ivmSelectVisualise.m: Visualise the selected point.
rocholExtract.m: Extract the lower triangular matrix from the Cholesky structure.
ggKernDiagCompute.m: Compute diagonal of GG kernel.
rbfKernExpandParam.m: Create kernel structure from RBF kernel's parameters.
modelPosteriorVar.m: variances of the posterior at points given by X.
ppcaPosteriorVar.m: Mean and variances of the posterior at points given by X.
cmpndKernGradient.m: Gradient of CMPND kernel's parameters.
mlpKernGradient.m: Gradient of MLP kernel's parameters.
isomapCreate.m: isomap embedding model.
gaussianwhiteKernParamInit.m: Gaussian white kernel parameter initialisation.
fgplvmPointObjective.m: Wrapper function for objective of a single point in latent space and the output location..
fgplvmCmu35Animate.m: Animate the test data jointly with predictions.
tensorKernParamInit.m: TENSOR kernel parameter initialisation.
lfmComputeH3VP.m: Helper function for computing part of the LFMVXLFM kernel.
lvmThreeDPlot.m: Helper function for plotting the labels in 3-D.
ivm.m: Initialise an IVM model.
invcmpndKernDisplay.m: Display parameters of the INVCMPND kernel.
gpOut.m: Evaluate the output of an Gaussian process model.
leCreate.m: Laplacian eigenmap model.
mogUpdateMean.m: Update the means of an MOG model.
rbfardjitKernDiagGradient.m: Compute the gradient of the RBFARD2 kernel's diagonal wrt parameters.
ivmDeconstruct.m: break IVM in pieces for saving.
mvuOptimise.m: Optimise an MVU model.
gpTimeDynamicsOut.m: Evaluate the output of GPTIMEDYNAMICS.
simKernExpandParam.m: Create kernel structure from SIM kernel's parameters.
noneKernParamInit.m: NONE kernel parameter initialisation.  
sdlfmXsdlfmKernComputeBlock.m: Computes SDLFM kernel matrix for block i,j
gpComputeObservationLogLikelihood.m:  
kbrParamInit.m: KBR model parameter initialisation.
sqexpKernDiagGradX.m: Gradient of SQEXP kernel's diagonal with respect to X.
linKernDiagGradX.m: Gradient of LIN kernel's diagonal with respect to X.
linearExpandParam.m: Update linear model with vector of parameters.
lfmvKernDiagCompute.m: Compute diagonal of LFMV kernel.
centeringMatrix.m: returns the centering matrix for the given dimensionality.
fgplvmNearestNeighbour.m: Give the number of errors in latent space for 1 nearest neighbour.
demStickFgplvm4.m: Model the stick man using an RBF kernel and 3-D latent space.
indexKernGradient.m: Gradient of INDEX kernel's parameters.
lleCreate.m: Locally linear embedding model.
ouKernDiagGradX.m: Gradient of OU kernel's diagonal with respect to t (see
bvhPlayData.m: Play bvh motion capture data.
demOptimiseGp.m: Shows that there is an optimum for the covariance function length scale.
linKernDisplay.m: Display parameters of the LIN kernel.
translateKernExpandParam.m: Create kernel structure from TRANSLATE kernel's parameters.
whiteblockKernDiagGradX.m: Gradient of WHITEBLOCK kernel's diagonal wrt X.
lmcKernDiagCompute.m: Compute the diagonal of the LMC kernel.
invSigmoid.m: The inverse of the sigmoid function.
lfmVisualiseWalking.m: Visualise the outputs in a latent force model
linardKernDiagGradX.m: Gradient of LINARD kernel's diagonal with respect to X.
lfmwhiteXwhiteKernGradient.m: Compute gradient between the LFM-WHITE and
rbfinfwhiteXwhiteKernCompute.m: Compute a cross kernel between the RBF-WHITE
rbfwhiteKernDiagGradX.m: Gradient of RBF-WHITE kernel's diagonal w.r.t. t.
indexardKernDiagCompute.m: Compute diagonal of INDEXARD kernel.
whitehKernExtractParam.m: Extract parameters from the WHITEH kernel structure.
ivmOptimiseNoise.m: Optimise the noise parameters.
expKernDiagGradX.m: Gradient of EXP kernel's diagonal with respect to X.
linearOutputGradX.m: Evaluate derivatives of linear model outputs with respect to inputs.
sdlfmvXsdlfmvKernGradient.m: Gradients of cross kernel between 2 SDLFM kernels.
demBrendanFgplvm2.m: Use the GP-LVM to model the Frey face data with FITC and back constraints.
gpLoadResult.m: Load a previously saved result.
sdlfmXsdrbfKernGradient.m: Cross gradient between a SDLFM and a SDRBF kernels.
lfmaKernExtractParam.m: Extract parameters from the LFMA kernel structure.
demVowelsFgplvm1.m: Model the vowels data with a 2-D FGPLVM using RBF kernel and back constraints.
mlpKernDisplay.m: Display parameters of the MLP kernel.
fgplvmExpandParam.m: Expand a parameter vector into a GP-LVM model.
lfmaXlfmKernCompute.m: Acceleration and position LFM kernel  
demStickFgplvm2.m: Model the stick man using an RBF kernel and dynamics.
rbfardKernDiagCompute.m: Compute diagonal of RBFARD kernel.
rbfwhiteKernCompute.m: Compute the RBF-WHITE kernel given the parameters, t1
sdlfmvKernExpandParam.m: Pass parameters from params to SDLFMV kernel
ncnmNoiseLikelihood.m: Likelihood of data under null category noise model.
lfmComputeH4AA.m: Helper function for computing part of the LFMAA kernel.
rbfKernGradX.m: Gradient of RBF kernel with respect to input locations.
numsf2str.m: Convert number to a string with a number of significant digits.
tensorKernSlash.m: Tensor kernel created by removing ith component.
whiteblockKernParamInit.m: WHITE BLOCK kernel parameter initialisation.
whitehKernParamInit.m: WHITEH kernel parameter initialisation.
expTransform.m: Constrains a parameter to be positive through exponentiation.
multiKernDiagCompute.m: Compute diagonal of MULTI kernel.
whiteXnoneKernCompute.m: Compute a cross kernel between WHITE and NONE kernels.
mlpDisplay.m: Display the multi-layer perceptron model.
demWalkSitJogDynamicsLearn.m: Learn the stick man dynamics.
treeFindChildren.m: Given a tree that lists only parents, add children.
rbfExtractParam.m: Wrapper for NETLAB's rbfpak.
dexpKernCompute.m: Compute the double exponential kernel,
gaussianKernDisplay.m: Display parameters of the GAUSSIAN kernel.
wienerKernDisplay.m: Display parameters of the WIENER kernel.
lfmvKernParamInit.m: LFMV kernel parameter initialisation. 
lfmResultsDynamic.m: Load a results file and visualise them.
cmpndNoisePointPlot.m: Plot the data-points for the CMPND noise model.
sqexpKernDisplay.m: Display parameters of the SQEXP kernel.
mlpardKernDisplay.m: Display parameters of the MLPARD kernel.
sdlfmvXsdrbfKernGradientBlock.m: Gradients of the parameters in block i,j
logisticNormalPriorExtractParam.m: Extract params from logistic-normal prior structure.
disimXrbfKernCompute.m: Compute a cross kernel between the DISIM and RBF kernels.
disimKernDiagGradient.m: Compute the gradient of the DISIM kernel's diagonal wrt parameters.
multiKernDiagGradX.m: Gradient of MULTI kernel's diagonal with respect to X.
lmcKernCompute.m: Compute the LMC kernel given the parameters and X.
linearLogLikelihood.m: Linear model log likelihood.
lfmjXrbfKernCompute.m: Compute cross kernel between the LFMJ and RBF kernels.
fgplvmLogLikeGradients.m: Compute the gradients for the FGPLVM.
rbfperiodicKernGradient.m: Gradient of RBFPERIODIC kernel's parameters.
ncnmNoise3dPlot.m: Draw a 3D or contour plot for the NCNM noise model.
linard2KernExpandParam.m: Create kernel structure from LINARD2 kernel's parameters.
ratquadKernDisplay.m: Display parameters of the RATQUAD kernel.
fgplvmObjectiveGradient.m: Wrapper function for FGPLVM objective and gradient.
modelWriteResult.m: Write a model to file.
lfmavGradientUpsilonMatrix.m: Gradient upsilon matrix accel. vel.
ngaussNoiseDisplay.m: Display parameters of the NGAUSS noise.
demSpgp1dPlot.m: Plot results from 1-D sparse GP.
polyKernDisplay.m: Display parameters of the POLY kernel.
viterbiAlign.m: Compute the Viterbi alignment.
rbfardKernParamInit.m: RBFARD kernel parameter initialisation.
ivmMeshVals.m: Give the output of the IVM for contour plot display.
noiseLikelihood.m: Return the likelihood for each point under the noise model.
sqexpKernGradX.m: Gradient of SQEXP kernel with respect to a point x.
mogMeanCov.m: Project a mixture of Gaussians to a low dimensional space.
linard2KernDiagCompute.m: Compute diagonal of LINARD2 kernel.
rbfinfwhiteKernGradient.m: Gradient of the parameters of the RBF-WHITE kernel
rbfKernParamInit.m: RBF kernel parameter initialisation.
rbfard2KernExpandParam.m: Create kernel structure from RBFARD2 kernel's parameters.
xyzhumanevaRemovePart.m:
velotransKernExpandParam.m: Create kernel structure from VELOTRANS kernel's parameters.
lfmXrbfKernCompute.m: Compute a cross kernel between the LFM and RBF kernels.
ivmNegLogLikelihood.m: Wrapper function for calling IVM likelihood.
kbrOutputGrad.m: Evaluate derivatives of KBR model outputs with respect to parameters.
ggwhiteXgaussianwhiteKernGradient.m: Compute gradient between the GG white
demRobotTraces1.m: Wireless Robot data from University of Washington, with tailored dynamics.
ncnmTwoDPlot.m: Make a 2-D plot of the null category noise model.
rbfExpandParam.m: Update rbf model with new vector of parameters.
gammaPriorLogProb.m: Log probability of Gamma prior.
noise3dPlot.m: Draw a 3D or contour plot for the relevant noise model.
ratquadKernCompute.m: Compute the RATQUAD kernel given the parameters and X.
gpSubspaceOut.m:
dnetExtractParam.m: Extract weights and biases from an DNET.
ivmEpUpdatePoint.m: Do an EP update of a point.
robThreeDynamicsDisplay.m: Display the robot dynamics model. 
rocCurve.m: Draw ROC curve and return labels.
linKernParamInit.m: LIN kernel parameter initialisation.
simwhiteKernExpandParam.m: Create kernel structure from SIM-WHITE kernel's
kernelCenter.m: Attempts to Center Kernel Matrix
lfmXrbfvKernCompute.m: Compute a cross kernel between the LFM and RBFV kernels.
mogProject.m: Project a mixture of Gaussians to a low dimensional space.
fgplvmTest.m: Test the gradients of the gpCovGrads function and the fgplvm models.
demRobotWirelessFgplvm1.m: Wireless Robot data from University of Washington, without dynamics and without back constraints.
gibbsKernExpandParam.m: Create kernel structure from GIBBS kernel's parameters.
acclaimPlayFile.m: Play motion capture data from a asf and amc file.
robTwoDynamicsLogLikeGradients.m: Gradients of the robot two dynamics wrt parameters.
polyardKernDiagCompute.m: Compute diagonal of POLYARD kernel.
getline.m: Get a line from a file.
fgplvmSequenceLogLikelihood.m: Log-likelihood of a sequence for the GP-LVM.
cmpndKernSetIndex.m: Set the indices in the compound kernel.
logdet.m: The log of the determinant when argument is positive definite.
mlpKernDiagCompute.m: Compute diagonal of MLP kernel.
rocholTransMultiply.m: Multiply by the transposed version of the rank one Cholesky.
rbfard2KernGradX.m: Gradient of RBFARD2 kernel with respect to input locations.
kernDiagCompute.m: Compute the kernel given the parameters and X.
ndsimKernExpandParam.m: Create kernel structure from NDSIM kernel's parameters.
nddisimKernExpandParamTransformSettings.m: Create kernel structure from DISIM kernel's parameters.
dynamicsTest.m: Run some tests on the specified dynamics model.
lfmGradientH31.m: Gradient of the function h_i(z) with respect to some of the
lfmKernDisplay.m: Display parameters of the LFM kernel.
linard2KernCompute.m: Compute the LINARD2 kernel given the parameters and X.
isomapReconstruct.m: Reconstruct an isomap form component parts.
simKernExtractParam.m: Extract parameters from the SIM kernel structure.
demTwoClusters1.m:
ivmKernelGradient.m: Gradient of likelihood approximation wrt kernel parameters.
srbfhKernGradient.m: Gradient of the parameters of a SRBFH kernel.
ouKernExpandParam.m: Create kernel structure from OU kernel's parameters
expKernGradX.m: Gradient of EXP kernel with respect to a point x.
lfmSample.m: Sample from LFM kernel
modelGetOutputWeights.m: Wrapper function to return output weight and bias matrices.
diagKernParamInit.m: DIAG kernel parameter initialisation.
lfmGradientSigmaH.m: Gradient of the function h_i(z) with respect \sigma.
gibbsperiodicKernCompute.m: Compute the GIBBSPERIODIC kernel given the parameters and X.
pskernelObjective.m: Likelihood approximation for point set IVM.
gaussianPriorGradient.m: Gradient wrt x of the log Gaussian prior.
gibbsperiodicKernDiagGradient.m: Compute the gradient of the GIBBSPERIODIC kernel's diagonal wrt parameters.
lfmCreate.m: Create a LFM model.
ggwhiteXggwhiteKernCompute.m: Compute a cross kernel between two GG white kernels.
ncnmContour.m: Special contour plot showing null category region.
kernExtractParamTransformSettings.m: Extract parameter transform settings from kernel structure.
ggKernDiagGradient.m: Compute gradient of the diagonal of GG kernel.
gibbsKernGradient.m: Gradient of GIBBS kernel's parameters.
whitefixedKernGradient.m: Gradient of WHITEFIXED kernel's parameters.
demSpgp1dGp3.m: Do a simple 1-D regression after Snelson & Ghahramani's example.
lfmvXlfmvKernGradient.m: Compute a cross gradient between a LFMV and a LFMV.
disimXdisimKernGradient.m: Compute a cross gradient between two DISIM kernels.
findAcyclicNeighbours.m: find the k nearest neighbours for each point in Y preventing cycles in the graph.
lfmwhiteKernCompute.m: Compute the LFM-WHITE kernel given the parameters, t1
fileKernExtractParam.m: Extract parameters from the FILE kernel structure.
sdlfmXsdrbfKernComputeBlock.m: Cross kernel between SDLFM and SDRBF for i,j
zeroAxes.m: A function to move the axes crossing point to the origin.
writeBoolToFID.m: Writes a boolean to an FID.
ggXgaussianKernGradX.m: Compute gradient between the GG and GAUSSIAN
lfmaKernDiagCompute.m: Compute diagonal of LFMAXLFMA kernel.
gpDynamicsSequenceLogLikelihood.m: Return the log likelihood of a given latent sequence.
xyzankurVisualise2.m:
mlpLogLikeHessian.m: Multi-layer perceptron Hessian.
getSymbols.m: Get a cell array of different plot symbols.
ggKernExtractParam.m: Extract parameters from the GG kernel structure.
rbfKernExtractParam.m: Extract parameters from the RBF kernel structure.
multiKernComputeBlock.m:
sdlfmKernMeanCovPartial.m: Helper function for derivatives in SDLFM kernel
demSwissRollLle3.m: Demonstrate LLE on the oil data.
dist1.m:	Calculates absolute distance (i.e. L1 norm) between two sets of
dexpKernDisplay.m: Display parameters of the double exponential kernel.
rbfperiodic2KernDisplay.m: Display parameters of the RBFPERIODIC2 kernel.
rbfinfwhiteKernCompute.m: Compute the RBF-WHITE kernel (with integration limits
multiKernParamInit.m: MULTI kernel parameter initialisation.
xyzankurError.m: Computes the error between two poses in xyz format 
xyzhumanevaGenerateMovie.m:
multimodelLogLikeGradients.m: Gradient of MULTIMODEL model log likelihood with respect to parameters.
demCmu35SequenceOptimise.m: 
cmpndNoiseLikelihood.m: Likelihood of the data under the CMPND noise model.
gpOptions.m: Return default options for GP model.
fgplvmOptimiseSequence.m: Optimise the postion of a latent sequence.
bvhConnectionMatrix.m: Compute the connection matrix for the structure.
modelSetOutputWeights.m: Wrapper function to return set output weight and bias matrices.
multiKernDisplay.m: Display parameters of the MULTI kernel.
gpCreate.m: Create a GP model with inducing varibles/pseudo-inputs.
lleEmbed.m: Embed data set with LLE.
lfmKernExpandParam.m: Create kernel structure from LFM kernel's parameters.
gpTimeDynamicsLogLikeGradients.m: Gradients of the GP dynamics wrt parameters.
mlpardKernExtractParam.m: Extract parameters from the MLPARD kernel structure.
gaussianKernExtractParam.m: Extract parameters from the gaussian kernel structure.
gaussianNoiseNuG.m: Compute nu and g for GAUSSIAN noise model.
lfmXlfmKernGradient.m: Compute a cross gradient between two LFM kernels.
xlogy.m: z = x*log(y) returns zero if x=y=0
ratquadKernDiagGradX.m: Gradient of RATQUAD kernel's diagonal with respect to X.
ardKernDiagCompute.m: Compute diagonal of ARD kernel.
ndsimKernExtractParamTransformSettings.m: Extract parameter transform settings from the SIM kernel structure.
robTwoDynamicsLogLikelihood.m: Give the log likelihood of the robot one dynamics part.
lfmGradientH41.m: Gradient of the function h_i(z) with respect to some of the
dnetEstep.m: Do an E-step (update importance weights) on an Density Network model.
ggKernCompute.m: Compute the GG kernel given the parameters and X.
dnetOptions.m: Options for a density network.
ouKernCompute.m: Compute the Ornstein-Uhlenbeck (OU) kernel arising from the
gibbsKernParamInit.m: GIBBS kernel parameter initialisation.
simwhiteKernParamInit.m: SIM-WHITE kernel parameter initialisation.
mvuReconstruct.m: Reconstruct an MVU form component parts.
sdlfmaXsdlfmKernCompute.m: Cross kernel between a SDLFMA and a SDLFM kernels.
rbfwhiteKernGradX.m: Gradient of RBF-WHITE kernel with respect to a point t.
linard2KernExtractParam.m: Extract parameters from the LINARD2 kernel structure.
demRobotWirelessFgplvm2.m: Wireless Robot data from University of Washington, without dynamics and without back constraints.
rbfhKernExtractParam.m: Extract parameters from the RBFH kernel structure.
logisticNormalPriorGradient.m: Gradient wrt x of the logistic-normal prior.
fgplvmReadFromFID.m: Load from a FID produced by the C++ implementation.
ouKernGradX.m: Gradient of OU kernel with respect to a point x (see
ncnmNoiseNuG.m: Update nu and g parameters associated with null category noise model.
demSwissRollLle4.m: Demonstrate LLE on the oil data.
mlpExtractParam.m: Extract weights and biases from an MLP.
bvhModify.m: Helper code for visualisation of bvh data.
optimiMinimize.m: Wrapper for Carl Rasmussen's minimize function.
lfmComputeUpsilonDiagVector.m: 
cmpndNoiseExpandParam.m: Create noise structure from CMPND noise's parameters.
leOptimise.m: Optimise an LE model.
kldivGaussian.m: Give the KL divergence between two Gaussians.
kernWriteParamsToFID.m: Write the kernel parameters to a stream.
vectorModify.m: Helper code for visualisation of vectorial data.
orderedNoiseOut.m: Compute the output of the ORDERED noise given the input mean and variance.
mgaussianNoiseParamInit.m: MGAUSSIAN noise parameter initialisation.
probitNoiseExpandParam.m: Create noise structure from PROBIT noise's parameters.
linardKernParamInit.m: LINARD kernel parameter initialisation.
sdlfmaXsdlfmKernComputeBlock.m: Computes SDLFM kernel matrix for block i,j
acclaimLoadChannels.m: Load the channels from an AMC file.
noiseCreate.m: Initialise a noise structure.
lfmGradientH42AP.m: Gradient of the function h_i(z) with respect to some of the
biasKernCompute.m: Compute the BIAS kernel given the parameters and X.
lfmComputeH4VP.m: Helper function for computing part of the LFMVXLFM kernel.
demOilFgplvm8.m: Oil data with variational sparse approximation.
dnetLogLikelihood.m: Density network log likelihood.
nddisimKernExtractParamTransformSettings.m: Extract parameter transform settings from the NDDISIM kernel structure.
demSwissRollFullLle2.m: Demonstrate LLE on the oil data.
lfmGradientSigmaH4AV.m: Gradient of the function h_i(z) with respect \sigma.
noiseExpectationLogLikelihood.m: Return the expectation of the log likelihood.
whiteXrbfKernGradient.m: Compute a cross gradient between WHITE and RBF kernels.
gpPointLogLikelihood.m: Log-likelihood of a test point for a GP.
modelAddDynamics.m: Add a dynamics kernel to the model.
sdlfmKernGradientConstant.m: Gradients for constants for the SDLFM kernel
rbfperiodicKernDiagGradient.m: Compute the gradient of the RBFPERIODIC kernel's diagonal wrt parameters.
ivmPrintPlot.m: Make a 3-D or contour plot of the IVM.
normuniPriorExpandParam.m: Expand Normal uniform prior structure from param vector.
demUnlabelledIvm2.m: Test IVM code on a toy crescent data.
lfmjaComputeUpsilonMatrix.m: Upsilon matrix jolt. accel. with t1, t2 limits
cmpndKernDiagGradX.m: Gradient of CMPND kernel's diagonal with respect to X.
demOil100Fgplvm1.m: Oil100 data with fully independent training conditional.
gibbsKernDiagGradX.m: Gradient of GIBBS kernel's diagonal with respect to X.
rbfperiodicKernDiagGradX.m: Gradient of RBFPERIODIC kernel's diagonal with respect to X.
gpComputeAlpha.m: Update the vector `alpha' for computing posterior mean quickly.
smallrandEmbed.m: Embed data set with small random values.
wienerKernParamInit.m: WIENER kernel parameter initialisation.
noiseExpandParam.m: Expand the noise model's parameters from params vector.
lvmClassClickVisualise.m: Callback function for visualising data in 2-D with clicks.
gpReversibleDynamicsExpandParam.m: Place the parameters vector into the model for GP dynamics.
treeFindParents.m: Given a tree that lists only children, add parents.
lfmGradientSigmaUpsilon.m: Gradient of the function \upsilon(z) with respect
kernGradient.m: Compute the gradient wrt the kernel's parameters.
gaussianNoiseOut.m: Compute the output of the GAUSSIAN noise given the input mean and variance.
polyKernDiagCompute.m: Compute diagonal of POLY kernel.
dexpKernGradient.m: Gradient of the double exponential kernel's parameters.
wienerKernExpandParam.m: Create kernel structure from WIENER kernel's parameters.
gammaPriorGradient.m: Gradient wrt x of the gamma prior.
vivmRunDataSet.m: Try the virtual IVM on a data set and save the results.
sparseKernDisplay.m: Display parameters of the SPARSE kernel.
multiKernGradientBlockX.m:
demOilLle3.m: Demonstrate LLE on the oil data.
ngaussNoiseLogLikelihood.m: Log likelihood of the data under the NGAUSS noise model.
lfmKernDiagGradient.m: Compute the gradient of the LFM kernel's diagonal wrt parameters.
demOilFgplvm2.m: Oil data with fully independent training conditional, and MLP back constraints.
doubleMatrixReadFromFID.m: Read a full matrix from an FID.
gpUpdateAD.m: Update the representations of A and D associated with the model.
rbfinfwhiteXwhiteKernGradient.m: Compute gradient between the RBF-WHITE kernel
modelGradient.m: Gradient of error function to minimise for given model.
indexKernParamInit.m: INDEX kernel parameter initialisation.
invcmpndKernSetIndex.m: Set the indices in the inv. compound kernel.
whiteblockKernGradX.m: Gradient of WHITEBLOCK kernel wrt input locations.
lfmLogLikelihood.m: Compute the log likelihood of a LFM model.
kernPriorGradient.m: Compute gradient terms associated with kernel priors.
linKernExtractParam.m: Extract parameters from the LIN kernel structure.
lfmExtractParam.m: Extract the parameters of an LFM model.
gpDynamicsSetLatentValues.m: Set the latent values inside the model.
mgaussianNoise3dPlot.m: Draws a 3D or contour plot for the MGAUSSIAN noise model.
biasKernGradient.m: Gradient of BIAS kernel's parameters.
lfmGradientH.m: Gradient of the function h_i(z) with respect to some of the
demClassificationTwoIvm1.m: IVM for classification on a data-set sampled from a GP
lfmXrbfKernGradient.m: Compute gradient between the LFM and RBF kernels.
demGpSample.m: Simple demonstration of sampling from a covariance function.
ardKernParamInit.m: ARD kernel parameter initialisation.
stack.m: Return column stacked vector of given matrix.
mlpardKernParamInit.m: MLPARD kernel parameter initialisation.
sheffieldmlToolboxes.m: External dependencies for Sheffield ML toolboxes.
mgaussianNoiseLikelihood.m: Likelihood of the data under the MGAUSSIAN noise model.
noisePointPlot.m: Plot the data-points for the given noise model.
gpReversibleDynamicsLogLikelihood.m: Give the log likelihood of the dynamics part.
rbfperiodic2KernExpandParam.m: Create kernel structure from RBFPERIODIC2 kernel's parameters.
fgplvmLoadResult.m: Load a previously saved result.
fgplvmDynamicsSample.m: Sample a field from the GP.
robTwoDynamicsLatentGradients.m: Gradients of the X vector given the dynamics model.
mogEstep.m: Do an E-step on an MOG model.
wangPriorLogProb.m: Log probability of Wang prior.
ivmUpdateSites.m: Update site parameters.
linearCreate.m: Create a linear model.
demProbit1.m: Test IVM code on a toy crescent data.
rbfardKernExtractParam.m: Extract parameters from the RBFARD kernel structure.
rocholBackSub.m: Backsubstitute the representation of the rank one Cholesky.
priorParamInit.m: Prior model's parameter initialisation.
demSilhouetteGp2.m: Model silhouette data with independent MLP GPs.
ratquadKernExtractParam.m: Extract parameters from the RATQUAD kernel structure.
lfmwhiteKernExpandParam.m: Create kernel structure from LFM-WHITE kernel's
mgaussianNoiseGradientParam.m: Gradient of MGAUSSIAN noise's parameters.
lfmwhiteKernDiagGradX.m: Gradient of LFM-WHITE kernel's diagonal w.r.t. t.
fgplvmSequenceObjectiveGradient.m: Wrapper function for objective
kernExtractParam.m: Extract parameters from kernel structure.
fileKernDisplay.m: Display parameters of the FILE kernel.
ouKernGradient.m: Gradient of OU kernel's parameters (see ouKernCompute or
gpReversibleDynamicsDisplay.m: Display a GP dynamics model.
simwhiteXrbfwhiteKernGradient.m: Compute a cross gradient between a SIM-WHITE
sheatKernGradient.m: Gradient of the parameters of a SHEAT kernel.
lfmGradientSigmaH4AA.m: Gradient of the function h_i(z) with respect \sigma.
escapeText.m: Add back slashes to escape existing backslashes in a text.
modelPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
gpScaleBiasGradient.m: Compute the log likelihood gradient wrt the scales.
tensorKernGradX.m: Gradient of TENSOR kernel with respect to a point x.
demCmu35Animate.m: Animate reconstructed right leg and body.
gpComputeTranslationLogLikelihood.m:  
demOil100Fgplvm2.m: Oil100 data with FGPLVM.
gammaPdf.m: PDF for the Gamma distribution.
ggwhiteKernExtractParam.m: Extract parameters from the GG WHITE kernel structure.
multimodelCreate.m: Create a MULTIMODEL model.
gammaPriorExtractParam.m: Extract params from gamma prior structure.
diagKernExpandParam.m: Create kernel structure from DIAG kernel's parameters.
multiKernExtractParamTransformSettings.m: Extract parameter transform settings 
demSwissRollFullLle1.m: Demonstrate LLE on the oil data.
dexpKernExpandParam.m: Create a kernel structure from the double exponential
fgplvmSequenceGradient.m: Wrapper function for gradient of a latent sequence.
noiseExtractParam.m: Extract the noise model's parameters.
lfmKernDiagCompute.m: Compute diagonal of LFM kernel.
rbfperiodicDisplay.m: Display parameters of the RBFPERIODIC model.
heatKernParamInit.m: HEAT kernel parameter initialisation.
lfmvpComputeUpsilonMatrix.m: Upsilon matrix vel. pos. with t1, t2 limits
distanceWarp.m: Dynamic Time Warping Algorithm
sdlfmaXsdlfmKernGradientBlock.m: Gradients of the parameters in block i,j
orderedNoiseParamInit.m: ORDERED noise parameter initialisation.
gpReadFromFile.m: Load a file produced by the C++ implementation.
rbfardKernGradient.m: Gradient of RBFARD kernel's parameters.
gpComputeM.m: Compute the matrix m given the model.
disimKernExtractParam.m: Extract parameters from the DISIM kernel structure.
lfmapGradientUpsilonVector.m: Gradient upsilon vector accel. pos.
gibbsperiodicKernExtractParam.m: Extract parameters from the GIBBSPERIODIC kernel structure.
whiteblockKernExpandParam.m: Fill WHITEBLOCK kernel structure with params.
dexpKernDiagGradX.m: Gradient of the double exponential kernel's diagonal
heatXrbfhKernGradient.m: Gradient wrt parameters between a HEAT and a RBFH.
lvmLoadData.m: Load a latent variable model dataset.
sqexpKernDiagGradient.m: Compute the gradient of the SQEXP kernel's diagonal wrt parameters.
kernWriteToFID.m: Load from an FID written by the C++ implementation.
lfmwhiteKernDiagCompute.m: Compute diagonal of LFM-WHITE kernel.
bvh2xyz.m: Compute XYZ values given structure and channels.
priorLogProb.m: Log probability of given prior.
demRegressionTwoIvm1.m: The data-set is sampled from a GP with known parameters.
gpObjectiveGradient.m: Wrapper function for GP objective and gradient.
modelOptimise.m: Optimise the given model.
laplacePriorExpandParam.m: Expand Laplace prior structure from param vector.
lvmClassVisualisePath.m: Latent variable model path drawing in latent space.
cmpndTieParameters.m: Tie parameters together.
polyKernExpandParam.m: Create kernel structure from POLY kernel's parameters.
spectralUpdateX.m: Update the latent representation for spectral model.
fgplvmResultsCpp.m: Load a results file and visualise them.
gibbsperiodicKernDiagGradX.m: Gradient of GIBBSPERIODIC kernel's diagonal with respect to X.
dnetOptimise.m: Optimise an DNET model.
fgplvmCreate.m: Create a GPLVM model with inducing variables.
