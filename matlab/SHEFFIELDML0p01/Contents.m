% GPMAT toolbox
% Version 0.01		20-Feb-2013
% Copyright (c) 2013, Neil D. Lawrence
% 
, Neil D. Lawrence
% SQEXPKERNEXTRACTPARAM Extract parameters from the SQEXP kernel structure.
% NONEKERNEXPANDPARAM Create kernel structure from NONE kernel's parameters.
% DISIMKERNGRADX Gradient of DISIM kernel with respect to a point x.
% DEMCMU35GPLVMRECONSTRUCT Reconstruct right leg and body of CMU 35.
% RBFPERIODICKERNEXPANDPARAM Create kernel structure from RBFPERIODIC kernel's parameters.
% DNETOUTPUTGRADX Evaluate derivatives of DNET model outputs with respect to inputs.
% FGPLVMADDCONSTRAINT Add latent constraints to FGPLVM model
% MODELLATENTGRADIENTS Gradients of the latent variables for dynamics models in the GPLVM.
% READDOUBLEFROMFID Read a double from an FID.
% MATERN32KERNCOMPUTE Compute the MATERN32 kernel given the parameters and X.
% LINARDKERNGRADIENT Gradient of LINARD kernel's parameters.
% TRANSLATEKERNPARAMINIT TRANSLATE kernel parameter initialisation.
% WHITEFIXEDXWHITEFIXEDKERNCOMPUTE Compute a cross kernel between two WHITEFIXED kernels.
% SDLFMAXSDRBFKERNCOMPUTE Cross kernel between a SDLFMA and a SDRBF kernels.
% LLEOPTIONS Options for a locally linear embedding.
% POLYARDKERNEXTRACTPARAM Extract parameters from the POLYARD kernel structure.
% GPREVERSIBLEDYNAMICSLOGLIKEGRADIENTS Gradients of the GP reversible dynamics wrt parameters.
% WHITEHKERNDIAGGRADX Gradient of WHITEH kernel's diagonal with respect to X.
% GGXGAUSSIANKERNGRADIENT Compute gradient between the GG and GAUSSIAN kernels.
% SPECTRALUPDATELAPLACIAN Update the Laplacian using graph connections.
% LMCKERNPARAMINIT LMC kernel parameter initialisation.
% NONEKERNDIAGGRADX Gradient of NONE kernel's diagonal with respect to X.
% RBFWHITEKERNEXPANDPARAM Create kernel structure from RBF-WHITE kernel's
% READBINARYDOUBLES Read information from a binary file in as doubles.
% GAUSSIANPRIOREXTRACTPARAM Extract params from Gaussian prior structure.
% ACCLAIMSEMILOADCHANNELS Load the channels from an AMC file for a subset
% FGPLVMEXTRACTPARAM Extract a parameter vector from a GP-LVM model.
% LMCKERNDIAGGRADIENT Gradient of the LMC kernel's diagonal wrt parameters.
% SIMWHITEXWHITEKERNGRADIENT Compute gradient between the SIM-WHITE and WHITE kernels.
% MODELOUT Give the output of a model for given X.
% RBFPERIODICOUTPUTGRADX Evaluate derivatives of a RBFPERIODIC model's output with respect to inputs.
% NOISEPARAMINIT Noise model's parameter initialisation.
% CUMGAUSSIAN Cumulative distribution for Gaussian.
% MLPOUTPUTGRAD Evaluate derivatives of mlp model outputs with respect to parameters.
% LFMGRADIENTSIGMAH3AA Gradient of the function h_i(z) with respect \sigma.
% PLOTMATRIXOPTIONS Default options for plot matrix.
% LMCKERNEXPANDPARAM Expands parameters into a LMC kernel structure.
% LEDECONSTRUCT break LE in pieces for saving.
% GIBBSKERNGRADX Gradient of GIBBS kernel with respect to input locations.
% ROBONEDYNAMICSCREATE Create the dynamics model. 
% GPGRADIENT Gradient wrapper for a GP model.
% IVM3DPLOT Make a 3-D or contour plot of the IVM.
% SPECTRUMMODIFY Helper code for visualisation of spectrum data.
% TRACEPRODUCT Returns the trace of the product of two matrices.
% FILEKERNGRADIENT Gradient of FILE kernel's parameters.
% DEMREGRESSIONONEIVM1 The data-set is sampled from a GP with known parameters.
% FGPLVMTESTMISSING Make sure missing data likelihood match full ones.
% VELOTRANSKERNDISPLAY Display parameters of the VELOTRANS kernel.
% TENSORKERNEXTRACTPARAM Extract parameters from the TENSOR kernel structure.
% IVMADDPOINT Add a point into the IVM representation.
% PROBITNOISELIKELIHOOD Likelihood of the data under the PROBIT noise model.
% SQEXPKERNDIAGCOMPUTE Compute diagonal of SQEXP kernel.
% NGAUSSNOISEEXPANDPARAM Create noise structure from NGAUSS noise's parameters.
% ORDEREDGRADX Gradient wrt x of log-likelihood for Ordered categorical model.
% WIENERKERNDIAGGRADX Gradient of WIENER kernel's diagonal with respect to X.
% RATQUADKERNDIAGGRADIENT Compute the gradient of the RATQUAD kernel's diagonal wrt parameters.
% COMPONENTKERNREADPARAMSFROMFID Read a component based kernel from a C++ file.
% GPREVERSIBLEDYNAMICSEXTRACTPARAM Extract parameters from the GP reversible dynamics model.
% CMPNDKERNEXPANDPARAMTRANSFORMSETTINGS Create kernel structure from CMPND kernel's parameter transformation settings.
% NDDISIMKERNGRADIENT Gradient of NDDISIM kernel's parameters.
% SKELVISUALISE For drawing a skel representation of 3-D data.
% TRANSLATEKERNEXTRACTPARAM Extract parameters from the TRANSLATE kernel structure.
% IVMAPPROXGRADX Returns the gradient of the approxmate log-likelihood wrt x.
% SDLFMKERNCOMPUTECONSTANT Compute constants for the SDLFM kernel
% MGAUSSIANNOISEEXPANDPARAM Create noise structure from MGAUSSIAN noise's parameters.
% RBFPERIODICKERNGRADX Gradient of RBFPERIODIC kernel with respect to a point x.
% GAUSSIANNOISESITES Update the site parameters for the GAUSSIAN noise mode.
% ROBONEDYNAMICSEXTRACTPARAM Extract parameters from the robot one dynamics model.
% LINARD2KERNDISPLAY Display parameters of the LINARD2 kernel.
% KERNDISPLAY Display the parameters of the kernel.
% MLPARDKERNGRADIENT Gradient of MLPARD kernel's parameters.
% RBFINFWHITEKERNEXTRACTPARAM Extract parameters from the RBF-WHITE kernel
% RBFARDJITKERNEXPANDPARAM Create kernel structure from RBFARDJIT kernel's parameters.
% STRINGSPLIT Return separate parts of a string.
% SDLFMXSDLFMVKERNCOMPUTEBLOCK Computes SDLFM kernel matrix for block i,j
% BIASKERNPARAMINIT BIAS kernel parameter initialisation.
% KERNGETVARIANCE Get the signal associated with a the kernel.
% DEMVOWELSISOMAP Model the vowels data with a 2-D FGPLVM using RBF kernel.
% WHITEKERNEXTRACTPARAM Extract parameters from the WHITE kernel structure.
% SHEATKERNCOMPUTE Compute a cross kernel between two SHEAT kernels.
% DISIMKERNGRADIENT Gradient of DISIM kernel's parameters.
% SQEXPKERNEXPANDPARAM Create kernel structure from SQEXP kernel's parameters.
% MLPARDKERNEXPANDPARAM Create kernel structure from MLPARD kernel's parameters.
% UNIFORMPRIORGRADIENT Gradient wrt x of the uniform prior.
% RBFARD2KERNEXTRACTPARAM Extract parameters from the RBFARD2 kernel structure.
% LAPLACEPRIORLOGPROB Log probability of Laplace prior.
% PPCAEMBED Embed data set with probabilistic PCA.
% IVMGRADX Returns the gradient of the log-likelihood wrt x.
% DEMROBOTWIRELESSFGPLVM4 Wireless Robot data from University of Washington with dynamics and back constraints.
% NOISEWRITETOFID Load from an FID written by the C++ implementation.
% LNDIFFERFS Helper function for computing the log of difference
% MULTIKERNGRADIENT Gradient of MULTI kernel's parameters.
% PPCAOUT Output of an PPCA model.
% KERNCORRELATION Compute the correlation matrix kernel given the parameters and X.
% LAPLACEPRIOREXTRACTPARAM Extract params from Laplace prior structure.
% GPRECONSTRUCT Reconstruct an GP form component parts.
% XYZHUMANEVAJOINT2POS
% BVHWRITEFILE Write a bvh file from a given structure and channels.
% MLPPARAMINIT Initialise the parameters of an MLP model.
% FGPLVMTAYLORANGLEERRORS Helper function for computing angle errors for CMU 35 data.
% SKEL2XYZ Compute XYZ values given skeleton structure and channels.
% HEATKERNDISPLAY Display parameters of the HEAT kernel.
% RBFWHITEXWHITEKERNGRADIENT Compute gradient between the RBF-WHITE and
% GPDYNAMICSEXTRACTPARAM Extract parameters from the GP dynamics model.
% TREEGETWIDTHS give width of each level of tree.
% SIMKERNCOMPUTE Compute the SIM kernel given the parameters and X.
% FGPLVMPRINTPLOT Print latent space for learnt model.
% LFMCOMPUTEH4AV Helper function for computing part of the LFMAV kernel.
% FINDDIRECTEDNEIGHBOURS find the k nearest neighbours for each point in Y preventing cycles in the graph.
% DISIMSAMPLE Sample from SIM kernel
% FGPLVMSCATTERPLOTCOLOR 2-D scatter plot of the latent points with color - for Swiss Roll data.
% GPDISPLAY Display a Gaussian process model.
% FGPLVMPOINTGRADIENT Wrapper function for gradient of a single point.
% RBFKERNDIAGGRADX Gradient of RBF kernel's diagonal with respect to X.
% CMPNDNOISEGRADIENTPARAM Gradient of CMPND noise's parameters.
% DNETRECONSTRUCT Reconstruct an DNET form component parts.
% GAUSSIANKERNEXPANDPARAM Create kernel structure from gaussian kernel's parameters.
% DEMSTICKFGPLVM1 Model the stick man using an RBF kernel.
% PRINTPLOT Print a plot to eps and png files.
% WANGPRIOREXTRACTPARAM Extract params from Wang prior structure.
% RBFPERIODICLOGLIKEGRADIENTS Gradient of RBFPERIODIC model log likelihood with respect to parameters.
% DIAGKERNGRADIENT Gradient of DIAG kernel's parameters.
% MATRIXREADFROMFID Read a matrix from an FID.
% RBFKERNDIAGCOMPUTE Compute diagonal of RBF kernel.
% SRBFHKERNCOMPUTE Compute an SRBFH kernel.
% GGXGAUSSIANKERNCOMPUTE Compute a cross kernel between the GG and GAUSSIAN kernels.
% SIMXSIMKERNDIAGCOMPUTE Diagonal of a cross kernel between two SIM kernels.
% GAUSSIANNOISEGRADVALS Gradient of GAUSSIAN noise log Z with respect to input mean and variance.
% WHITEKERNDIAGGRADIENT Compute the gradient of the WHITE kernel's diagonal wrt parameters.
% NGAUSSIAN Compute a Gaussian with mean 0 and variance 1.
% NGAUSSNOISELIKELIHOOD Likelihood of the data under the NGAUSS noise model.
% XYZPOPPEANIM Animate point cloud of stick man from Poppe dataset.
% SIMKERNDISPLAY Display parameters of the SIM kernel.
% PCAEMBED Embed data set with PCA.
% LFMLOGLIKEGRADIENTS Compute the gradients of the log likelihood of a LFM model.
% DEMREGRESSIONGP Demonstrate Gaussian processes for regression.
% IMAGEVISUALISE Helper code for showing an image during 2-D visualisation.
% RBFDISPLAY Display an RBF network.
% INVCMPNDKERNEXPANDPARAM Create kernel structure from INVCMPND kernel's parameters.
% DEMOILFGPLVM5 Oil data with partially independent training conditional.
% ROBTHREEDYNAMICSCREATE Create the dynamics model. 
% GPDYNAMICSSEQUENCELOGLIKEGRADIENT Log-likelihood gradient for of a sequence of the GP-LVM dynamics.
% RATQUADKERNGRADX Gradient of RATQUAD kernel with respect to input locations.
% DISIMKERNDISPLAY Display parameters of the DISIM kernel.
% FGPLVMPOSTERIORVAR Variances of the posterior at points given by X.
% NORMUNIPRIORLOGPROB Log probability of a normal uniform.
% IVMLOGLIKELIHOOD Return the log-likelihood for the IVM.
% LFMAXRBFKERNGRADIENT Compute gradient between the LFMA and RBF kernels.
% SIMWHITEKERNDIAGCOMPUTE Compute the diagonal of the SIM-WHITE kernel.
% SCALENOISEPARAMINIT Scale noise model's parameter initialisation.
% KBREXPANDPARAM Create model structure from KBR model's parameters.
% EXPKERNPARAMINIT EXP kernel parameter initialisation.
% PROBIT3DPLOT Draw a 3D or contour plot for the probit.
% DEMBRENDANFGPLVM1 Use the GP-LVM to model the Frey face data with FITC.
% SIMWHITEKERNDISPLAY Display parameters of the SIM-WHITE kernel.
% LINARD2KERNPARAMINIT LINARD2 kernel parameter initialisation.
% LVMVISUALISEGENERAL Visualise the manifold.
% MULTIMODELDISPLAY Display parameters of the MULTIMODEL model.
% LFMGRADIENTSIGMAH3AV Gradient of the function h_i(z) with respect \sigma.
% SDLFMAXSDLFMAKERNGRADIENTBLOCK Gradients of the parameters in block i,j
% MLPOPTIMISE Optimise MLP for given inputs and outputs.
% MLPOPTIONS Options for the multi-layered perceptron.
% RBFOUTPUTGRAD Evaluate derivatives of rbf model outputs with respect to parameters.
% FILEKERNREAD Read kernel values from file or cache.
% LMCKERNDIAGGRADX Gradient of LMC kernel's diagonal with respect to X.
% SIMWHITEXRBFINFWHITEKERNGRADIENT Compute a cross gradient between a
% DNETDECONSTRUCT break DNET in pieces for saving.
% ROBTHREEDYNAMICSSETLATENTVALUES Set the latent values inside the model.
% DEMFOURWALKSRECONSTRUCT Reconstruct right leg of CMU 35.
% DEMOILFGPLVM7 Oil data with variational sparse approximation.
% BIASKERNDIAGGRADIENT Compute the gradient of the BIAS kernel's diagonal wrt parameters.
% LFMJXLFMVKERNCOMPUTE Jolt and velocity LFM kernel  
% DEMSILHOUETTEGP1 Model silhouette data with independent RBF GPs.
% GAUSSIANPRIORLOGPROB Log probability of Gaussian prior.
% LAPLACEPRIORGRADIENT Gradient wrt x of the log Laplace prior.
% WHITEBLOCKKERNDIAGGRADIENT WHITEBLOCK kernel's diagonal gradient wrt par.
% IVMRUN Run the IVM on a given data set.
% XYZPOPPEDRAW Helper function for drawing data from Poppe.
% LFMGRADIENTUPSILON Gradient of the function \upsilon(z) with respect to
% HEATKERNGRADIENT Gradient of HEAT kernel's parameters.
% MATERN32KERNDISPLAY Display parameters of the MATERN32 kernel.
% FGPLVMOPTIMISE Optimise the FGPLVM.
% LFMVISUALISE Visualise the outputs in a latent force model
% GPREVERSIBLEDYNAMICSCREATE Create a reversible dynamics model. 
% DEMSPGP1DGP2 Do a simple 1-D regression after Snelson & Ghahramani's example.
% NGAUSSNOISE3DPLOT Draws a 3D or contour plot for the NGAUSS noise model.
% MLPCREATE Multi-layer peceptron model.
% MODELTEST Run some tests on the specified model.
% FGPLVMVISUALISE Visualise the manifold.
% LFMAXLFMKERNGRADIENT Compute a cross gradient between a LFMA and a LFM.
% IDENTITYTRANSFORM 
% CMDSROADDATA This script uses classical MDS to visualise some road distance data.
% GAUSSIANKERNPARAMINIT Gaussian kernel parameter initialisation.
% LFMJXLFMAKERNCOMPUTE Jolt and acceleration LFM kernel  
% MODELSAMP Give a sample from a model for given X.
% LFMCOMPUTEH3 Helper function for computing part of the LFM kernel.
% DEMCMU35GPLVM4 Learn a GPLVM on CMU 35 data set.
% WANGPRIORGRADIENT Gradient wrt x of the Wang prior.
% MULTIKERNEXPANDPARAM Create kernel structure from MULTI kernel's parameters.
% RBFARD2KERNDIAGGRADX Gradient of RBFARD2 kernel's diagonal with respect to X.
% SCG2 Scaled conjugate gradient optimization like netlab's scg, with a slight modification for speeding it up.
% ARDKERNGRADIENT Gradient of ARD kernel's parameters.
% GPCOVGRADSTEST Test the gradients of the likelihood wrt the covariance.
% WHITEKERNEXPANDPARAM Create kernel structure from WHITE kernel's parameters.
% XGAMRND Draw a sample from the gamma distribution.
% LINARDKERNDIAGCOMPUTE Compute diagonal of LINARD kernel.
% MATERN32KERNPARAMINIT MATERN32 kernel parameter initialisation.
% MODELREADFROMFID Load from a FID produced by C++ code.
% SDLFMKERNCOMPUTE Compute the SDLFM kernel given the parameters and X.
% GPREADFROMFID Load from a FID produced by the C++ implementation.
% RBFARDJITKERNGRADX Gradient of RBFARDJIT kernel with respect to input locations.
% DEMFOURWALKS1 Model four seperate walsk using an RBF kernel and dynamics.
% RATQUADKERNDIAGCOMPUTE Compute diagonal of RATQUAD kernel.
% RBFXNONEKERNCOMPUTE Compute a cross kernel between RBF and NONE kernels.
% LFMVXLFMKERNGRADIENT Compute a cross gradient for a LFMVXLFM.
% DEMSWISSROLLFULLLLE3 Demonstrate LLE on the oil data.
% LINARD2KERNGRADX Gradient of LINARD2 kernel with respect to input locations.
% LFMVKERNEXTRACTPARAM Extract parameters from the LFMV kernel structure.
% KERNPRIORLOGPROB Compute penalty terms associated with kernel priors.
% GAUSSIANNOISEPOINTPLOT Plot the data-points for the GAUSSIAN noise model.
% SDLFMVXSDLFMKERNCOMPUTE Cross kernel between a SDLFMV and a SDLFM kernels.
% NCNMNOISEEXPANDPARAM Expand null category noise model's structure from param vector.
% MOGPRINTPLOT Print projection of MOG into two dimensions.
% DISIMXDISIMKERNCOMPUTE Compute a cross kernel between two DISIM kernels.
% GAUSSIANKERNDIAGCOMPUTE Compute diagonal of gaussian kernel.
% RBFWHITEKERNPARAMINIT RBF-WHITE kernel parameter initialisation. The RBF-
% DEMSILHOUETTEPLOT
% CMPNDNOISEGRADVALS Gradient of CMPND noise log Z with respect to input mean and variance.
% DEMSTICKGP1 Demonstrate Gaussian processes for regression on stick man data.
% RBFARD2KERNDIAGCOMPUTE Compute diagonal of RBFARD2 kernel.
% RATQUADKERNPARAMINIT RATQUAD kernel parameter initialisation.
% LFMCOMPUTEH4JP Helper function for computing part of the LFMJP kernel.
% MGAUSSIANNOISEOUT Compute the output of the MGAUSSIAN noise given the input mean and variance.
% MODELCREATE Create a model of the specified type.
% LNCUMGAUSSSUM The log of the weighted sum of two cumulative Gaussians.
% GPSUBSPACEEXPANDPARAM 
% MLPARDKERNDIAGGRADX Gradient of MLPARD kernel's diagonal with respect to X.
% WRITEINTTOFID Writes an integer to an FID.
% SDLFMVKERNDIAGCOMPUTEBLOCK Diagonal of a SDLFM kernel matrix for block i
% XYZHUMANEVAHEADINGANGLE 
% DNETUPDATEOUTPUTWEIGHTS Do an M-step (update parameters) on an Density Network model.
% SDLFMXSDLFMKERNGRADIENT Gradients of cross kernel between 2 SDLFM kernels.
% LERECONSTRUCT Reconstruct an LE form component parts.
% RBFARDKERNDIAGGRADX Gradient of RBFARD kernel's diagonal with respect to X.
% GIBBSKERNDISPLAY Display parameters of the GIBBS kernel.
% MLPARDKERNCOMPUTE Compute the MLPARD kernel given the parameters and X.
% RBFINFWHITEXRBFINFWHITEKERNGRADIENT Compute a cross gradient between two
% KERNFACTORS Extract factors associated with transformed
% IVMUPDATEM Update matrix M, L, v and mu.
% LFMKERNDIAGGRADX Gradient of LFM kernel's diagonal with respect to X.
% GPPOSTERIORGRADMEANCOVAR Gadient of the mean and variances of the posterior at points given by X.
% LFMAACOMPUTEUPSILONMATRIX Upsilon matrix acce. accel. with t1, t2 limits
% PRIORWRITETOFID Write a prior to a C++ stream.
% DEMEP1 Demonstrate Expectation propagation on a toy data set..
% ROBTWODYNAMICSEXPANDPARAM Place the parameters vector into the model for first robot dynamics.
% SDLFMVXSDRBFKERNGRADIENT Gradients cross kernel between a SDLFM and SDRBF
% PRINTLATEXOPTIONS Options for printing a plot to LaTeX.
% LFMJXLFMKERNCOMPUTE Jolt and position LFM kernel  
% STICKVISUALISE For drawing a stick representation of 3-D data.
% MULTIKERNDIAGGRADIENT Compute the gradient of the MULTI kernel's diagonal wrt parameters.
% SIMWHITEXRBFWHITEKERNCOMPUTE Compute a cross kernel between a SIM-WHITE
% FGPLVMPOINTSAMPLELOGLIKELIHOOD
% MODELTIEPARAM Tie parameters of a model together.
% WIENERKERNCOMPUTE Compute the WIENER kernel given the parameters and X.
% GIBBSPERIODICKERNGRADIENT Gradient of GIBBSPERIODIC kernel's parameters.
% LVMRESULTSDYNAMIC Load a results file and visualise them.
% NGAUSSNOISEGRADIENTPARAM Gradient of NGAUSS noise's parameters.
% LFMAPCOMPUTEUPSILONVECTOR Upsilon vector for acce. pos. with t1 limit
% GAUSSIANPRIOREXPANDPARAM Expand Gaussian prior structure from param vector.
% GPDYNAMICSCREATE Create the dynamics model. 
% NDDISIMKERNEXPANDPARAM Create kernel structure from NDDISIM kernel's parameters.
% ORDEREDNOISEEXTRACTPARAM Extract parameters from the ORDERED noise structure.
% SDLFMKERNEXPANDPARAM Pass parameters from params to SDLFM kernel
% TENSORKERNDIAGGRADX Gradient of TENSOR kernel's diagonal with respect to X.
% NOISEUPDATENUG Update nu and g for a given noise model.
% GGWHITEKERNEXPANDPARAM Create kernel structure from GG white kernel's parameters.
% NGAUSSNOISENUG Update nu and g parameters associated with noiseless Gaussian noise model.
% ACCLAIM2XYZ Compute XYZ values given skeleton structure and channels.
% LFMCOMPUTEH4AP Helper function for computing part of the LFMAP kernel.
% SDLFMXSDLFMKERNGRADIENTBLOCKIGJ 
% ROCHOLFORESUB Foreward substitute the representation of the rank one Cholesky.
% DISIMXSIMKERNCOMPUTE Compute a cross kernel between DISIM and SIM kernels.
% SDRBFKERNCOMPUTE Compute the SDRBF kernel given the parameters and t1.
% LVMSCATTERPLOTNOVAR 2-D scatter plot of the latent points.
% TREEFINDLEAVES Return indices of all leaf nodes in a tree structure.
% LEOPTIONS Options for a Laplacian eigenmaps.
% RBFPERIODICEXTRACTPARAM Extract parameters from the RBFPERIODIC model structure.
% LMCKERNEXTRACTPARAM Extract parameters from the LMC kernel struc.
% IVMVIRTUAL Create virtual data points with the specified invariance.
% WHITEKERNDISPLAY Display parameters of the WHITE kernel.
% SIMCOMPUTETEST Test the file simComputeH.
% RBFPERIODICKERNPARAMINIT RBFPERIODIC kernel parameter initialisation.
% SDLFMVKERNDISPLAY Display parameters of the SDLFMV kernel.
% ARDKERNCOMPUTE Compute the ARD kernel given the parameters and X.
% SPRINGDAMPERSMODIFY Helper code for visualisation of springDamper data.
% WHITEFIXEDKERNDIAGCOMPUTE Compute diagonal of WHITEFIXED kernel.
% PARSEWIRELESSDATA Load wireless strength data.
% FGPLVMBACKCONSTRAINTGRAD Gradient with respect to back constraints if present.
% NCNMNOISEEXTRACTPARAM Extract parameters from null category noise model.
% RBFARD2KERNDISPLAY Display parameters of the RBFARD2 kernel.
% XYZHUMANEVAALIGN
% MOGOPTIMISE Optimise an MOG model.
% CMPNDNOISESITES Site updates for compound noise model.
% ARDKERNEXTRACTPARAM Extract parameters from the ARD kernel structure.
% INDEXARDKERNGRADIENT Gradient of INDEXARD kernel's parameters.
% RBFPERIODICLOGLIKELIHOOD Log likelihood of RBFPERIODIC model.
% LFMGRADIENTH41VP Gradient of the function h_i(z) with respect to some of the
% LFMAXLFMAKERNCOMPUTE Acceleration and acceleration LFM kernel  
% SDLFMKERNEXTRACTPARAM Extract parameters from the SDLFM kernel structure.
% DEMROBOTWIRELESSNAVIGATE Take some test data for the robot and navigate with it.
% MOCAPPARSETEXT Parse a motion capture text file.
% RBFPERIODICKERNEXTRACTPARAM Extract parameters from the RBFPERIODIC kernel structure.
% RBFKERNDIAGGRADIENT Compute the gradient of the RBF kernel's diagonal wrt parameters.
% DEMCMU35GPLVMFGPLVM2 Learn a GPLVM on CMU 35 data set.
% SIMKERNDIAGCOMPUTE Compute diagonal of SIM kernel.
% DEFAULTOPTIONS The default options for optimisation.
% GGKERNPARAMINIT GG kernel parameter initialisation.
% FGPLVMPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% LFMCOMPUTEH3JA Helper function for computing part of the LFMJA kernel.
% CMPNDKERNGRADX Gradient of CMPND kernel with respect to a point x.
% DISIMCOMPUTEHPRIME Helper function for comptuing part of the DISIM kernel.
% STICKMODIFY Helper code for visualisation of a stick man.
% MATERN52KERNGRADX Gradient of MATERN52 kernel with respect to input locations.
% LVMVISUALISE Visualise the manifold.
% MULTIMODELEXTRACTPARAM Extract parameters from the MULTIMODEL model structure.
% NORMUNIPRIORGRADIENT Gradient wrt x of the log normal uniform prior.
% MULTIKERNEXPANDPARAMTRANSFORMSETTINGS Create kernel structure from MULTI kernel's parameter transformation settings.
% GPOPTIMISE Optimise the inducing variable based kernel.
% PRIORTEST Run some tests on the specified prior.
% NOISEOUT Give the output of the noise model given the mean and variance.
% SDLFMXSDRBFKERNGRADIENTBLOCKIEJ 
% RBFPERIODIC2KERNEXTRACTPARAM Extract parameters from the RBFPERIODIC2 kernel structure.
% SDLFMAXSDRBFKERNGRADIENT Gradients cross kernel between a SDLFM and SDRBF
% SMOOTHANGLECHANNELS Try and remove artificial discontinuities associated with angles.
% ORDEREDNOISEPOINTPLOT Plot the data-points for the ORDERED noise model.
% RBFARDKERNDISPLAY Display parameters of the RBFARD kernel.
% HESSIANCHECK Check Hessian of objective function.
% INVCMPNDKERNGRADIENT Gradient of INVERSE-PRESICION-CMPND kernel's parameters.
% SCALENOISEEXPANDPARAM Expand Scale noise structure from param vector.
% LFMAVCOMPUTEUPSILONMATRIX Upsilon matrix acce. vel. with t1, t2 limits
% LFMCOMPUTEH3VV Helper function for computing part of the LFMVXLFMV kernel.
% WHITEBLOCKKERNEXTRACTPARAM Extract parameters from WHITEBLOCK kernel str.
% FGPLVMPOINTOBJECTIVEGRADIENT Wrapper function for objective and gradient of a single point in latent space and the output location..
% GPCOVGRADS Sparse objective function gradients wrt Covariance functions for inducing variables.
% MOCAPRESULTSCPPBVH Load results from cpp file and visualise as a bvh format.
% SIMWHITEKERNCOMPUTE Compute the SIM-WHITE kernel given the parameters, t1
% IVMCOMPUTELANDM Compute the L and M matrix.
% GPEXTRACTPARAM Extract a parameter vector from a GP model.
% POLYKERNEXTRACTPARAM Extract parameters from the POLY kernel structure.
% ROBTWODYNAMICSCREATE Create the dynamics model. 
% ORDEREDNOISEGRADVALS Gradient of ORDERED noise log Z with respect to input mean and variance.
% PROBITNOISEGRADVALS Gradient of PROBIT noise log Z with respect to input mean and variance.
% POLYKERNCOMPUTE Compute the POLY kernel given the parameters and X.
% RBFPERIODICOPTIONS Create a default options structure for the RBFPERIODIC model.
% WHITEHKERNDISPLAY Display parameters of the WHITEH kernel.
% DEMCMU35GPLVMRECONSTRUCTTAYLOR Reconstruct right leg of CMU 35.
% XYZMODIFY Update visualisation of skeleton data.
% MATERN32KERNDIAGGRADX Gradient of MATERN32 kernel's diagonal with respect to X.
% RBFINFWHITEKERNGRADX Gradient of RBF-WHITE kernel (with integration limits
% LFMJVCOMPUTEUPSILONMATRIX Upsilon matrix jolt. vel. with t1, t2 limits
% IVMOPTIMISEKERNEL Optimise the kernel parameters.
% SDLFMJMEANCOMPUTE Jolt mean for the switching dynamical LFM model.
% ROBONEDYNAMICSSETLATENTVALUES Set the latent values inside the model.
% IVMEPUPDATEM Update matrix M, L, varSigma and mu for EP.
% MULTIKERNEXTRACTPARAM Extract parameters from the MULTI kernel structure.
% KERNGRADX Compute the gradient of the kernel wrt X.
% LINEAROUT Obtain the output of the linear model.
% GGWHITEXWHITEKERNGRADIENT Compute gradient between the GGWHITE and WHITE kernels.
% IVMCOMPUTEINFOCHANGE Compute the information change associated with each point.
% SHEATKERNDIAGGRADIENT Gradient of the parameters of diagonal of a SHEAT kernel.
% LFMWHITEXRBFWHITEKERNGRADIENT Compute a cross gradient between an LFM-WHITE
% LFMAAGRADIENTSIGMAUPSILONMATRIX Gradient of upsilon matrix aa wrt sigma
% IVMCONTOUR Special contour plot showing decision boundary.
% TRANSLATEKERNGRADX Gradient of TRANSLATE kernel with respect to a point x.
% SQEXPKERNPARAMINIT SQEXP kernel parameter initialisation.
% FGPLVMADDDYNAMICS Add a dynamics kernel to the model.
% RBFARDJITKERNCOMPUTE Compute the RBFARD kernel given the parameters and X.
% IVMDOWNDATEM Remove point from M, L, mu and varSigma.
% LMCKERNDISPLAY Display parameters of the LMC kernel.
% WRITESTRINGTOFID Writes a string to an FID.
% LMCKERNGRADIENT Gradient of LMC kernel's parameters.
% DEMOILFGPLVM9 Oil data with three dimensions and variational sparse approximation.
% XYZANKURANIMCOMPAREMULTIPLE Animate many predictions and ground truth for
% MULTIKERNTEST Run some tests on the multiple output block kernel.
% WHITEKERNDIAGGRADX Gradient of WHITE kernel's diagonal with respect to X.
% DEMOILLLE1 Demonstrate LLE on the oil data.
% LINARDKERNDISPLAY Display parameters of the LINARD kernel.
% SDLFMXSDLFMVKERNGRADIENTBLOCKIGJ 
% WHITEHKERNGRADIENT Gradient of WHITEH kernel's parameters.
% ROCHOLMULTIPLY Multiply by the rank one Cholesky.
% GRADLNDIFFERFS Compute the gradient of the log difference of two erfs.
% BIASKERNDISPLAY Display parameters of the BIASkernel.
% LINEAROUTPUTGRAD Evaluate derivatives of linear model outputs with respect to parameters.
% MATERN32KERNDIAGGRADIENT Compute the gradient of the MATERN32 kernel's diagonal wrt parameters.
% GAMMAPRIORPARAMINIT Gamma prior model's parameter initialisation.
% SWISSROLLSCATTER 3-D scatter plot with colors.
% SDRBFKERNEXPANDPARAM Pass parameters from params to SDRBF kernel
% FGPLVMDYNAMICSPLOT 2-D scatter plot of the latent points.
% DEMSWISSROLLLLE2 Demonstrate LLE on the oil data.
% DEMMPPCA1 Demonstrate MPPCA on a artificial dataset.
% DEMOILFGPLVM3 Oil data with deterministic training conditional.
% GENERATECRESCENTDATA Generate crescent data.
% RBFARDJITKERNDIAGCOMPUTE Compute diagonal of RBFARDJIT kernel.
% GIBBSKERNDIAGCOMPUTE Compute diagonal of GIBBS kernel.
% DEG2RAD Transform degrees to radians.
% LINEARDISPLAY Display a linear model.
% DEMOILLLE4 Demonstrate LLE on the oil data.
% LFMCOMPUTEH3AA Helper function for computing part of the LFMAA kernel.
% WIENERKERNGRADIENT Gradient of WIENER kernel's parameters.
% LFMAKERNPARAMINIT LFMA kernel parameter initialisation. 
% LFMGRADIENTH42 Gradient of the function h_i(z) with respect to some of the
% GAUSSIANWHITEKERNDISPLAY Display parameters of the GAUSSIAN white kernel.
% FGPLVMRESULTSDYNAMIC Load a results file and visualise them.
% NEGLOGLOGIT Function which returns the negative log of the logistic function.
% NCNMLOADDATA Load a dataset.
% RBFPERIODIC2KERNCOMPUTE Compute the RBFPERIODIC2 kernel given the parameters and X.
% RBFPERIODICOUTPUTGRAD Evaluate derivatives of RBFPERIODIC model outputs with respect to parameters.
% GGWHITEKERNDIAGCOMPUTE Compute diagonal of GG WHITE kernel.
% BVHVISUALISE For updating a bvh representation of 3-D data.
% GAUSSIANWHITEKERNEXPANDPARAM Create kernel structure from gaussian white 
% WHITEFIXEDKERNPARAMINIT WHITEFIXED kernel parameter initialisation.
% LFMAPGRADIENTSIGMAUPSILONVECTOR Gradient of upsilon vector ap wrt sigma
% MODELOUTPUTGRADX Compute derivatives with respect to model inputs of model outputs.
% RBFWHITEKERNDIAGCOMPUTE Compute diagonal of RBF-WHITE kernel.
% FGPLVMDECONSTRUCT break FGPLVM in pieces for saving.
% GGWHITEXGAUSSIANWHITEKERNGRADX Compute gradient between the GG white and
% SDLFMXSDRBFKERNCOMPUTE Cross kernel between a SDLFM and a SDRBF kernels.
% NGAUSSNOISEPARAMINIT NGAUSS noise parameter initialisation.
% MOGCREATE Create a mixtures of Gaussians model.
% RBFARDKERNEXPANDPARAM Create kernel structure from RBFARD kernel's parameters.
% SIMWHITEKERNGRADX Gradient of SIM-WHITE kernel with respect to a point t.
% LFMWHITEXWHITEKERNCOMPUTE Compute a cross kernel between the LFM-WHITE
% CMPNDNOISEOUT Compute the output of the CMPND noise given the input mean and variance.
% SDLFMAXSDLFMKERNGRADIENTBLOCKIEJ 
% EXPKERNGRADIENT Gradient of EXP kernel's parameters.
% NDDISIMKERNDISPLAY Display parameters of the NDDISIM kernel.
% LFMVXRBFVKERNCOMPUTE Compute a cross kernel between the LFMV and RBFV kernels.
% LFMGRADIENTSIGMAH4AP Gradient of the function h_i(z) with respect \sigma.
% RBFPERIODICKERNDISPLAY Display parameters of the RBFPERIODIC kernel.
% LFMVXLFMKERNCOMPUTE Velocity and position LFM kernel  
% WHITEBLOCKKERNDIAGCOMPUTE Compute diagonal of WHITEBLOCK kernel.
% PROBITNOISEGRADIENTPARAM Gradient of PROBIT noise's parameters.
% INDEXKERNCOMPUTE Compute the INDEX kernel given the parameters and X.
% LFMGRADIENTSIGMAH3VP Gradient of the function h_i(z) with respect \sigma.
% KERNSETWHITE Helper function to set the white noise in a kernel if it exists.
% LVMLOADRESULT Load a previously saved result.
% DNETWRITERESULT Write a DNET result.
% GAUSSIANNOISEPARAMINIT GAUSSIAN noise parameter initialisation.
% LLEDECONSTRUCT break LLE in pieces for saving.
% DIAGKERNDIAGGRADX Gradient of DIAG kernel's diagonal with respect to X.
% CMPNDKERNCOMPUTE Compute the CMPND kernel given the parameters and X.
% DEMSPGP1DGP1 Do a simple 1-D regression after Snelson & Ghahramani's example.
% DEMSILHOUETTEPLOTTRUE Plot the true poses for the silhouette data.
% PROBITNOISEPOINTPLOT Plot the data-points for the PROBIT noise model.
% MODELEXTRACTPARAM Extract the parameters of a model.
% MLPKERNGRADX Gradient of MLP kernel with respect to input locations.
% SDLFMVKERNGRADIENT Gradient of SDLFM kernel's parameters.
% LFMWHITEXLFMWHITEKERNCOMPUTE Compute a cross kernel between two LFM-WHITE
% GPLOGLIKEGRADIENTS Compute the gradients for the parameters and X.
% RBFKERNCOMPUTE Compute the RBF kernel given the parameters and X.
% MODELGRADIENTCHECK Check gradients of given model.
% WHITEFIXEDKERNCOMPUTE Compute the WHITEFIXED kernel given the parameters and X.
% RBFARD2KERNCOMPUTE Compute the RBFARD kernel given the parameters and X.
% LFMVXRBFKERNGRADIENT Compute gradient between the LFMV and RBF kernels.
% LFMWHITECOMPUTEGRADTHETAH2 computes a portion of the LFM-WHITE kernel's gradient w.r.t. theta.
% WHITEFIXEDKERNEXTRACTPARAM Extract parameters from the WHITEFIXED kernel structure.
% GGWHITEKERNDISPLAY Display parameters of the GG WHITE kernel.
% MLPKERNCOMPUTE Compute the MLP kernel given the parameters and X.
% PDINV Invert a positive definite matrix.
% XYZANKURMODIFY  Helper function for modifying the point cloud from Agarwal and Triggs data.
% DEMSWISSROLLFULLLLE5 Demonstrate LLE on the oil data.
% MULTIKERNCOMPUTE Compute the MULTI kernel given the parameters and X.
% VIVMUSPSRESULTS Summarise the USPS result files in LaTeX.
% LFMCOMPUTEH3JV Helper function for computing part of the LFMJV kernel.
% SIMXSIMKERNGRADIENT Compute a cross gradient between two SIM kernels.
% ACCLAIMNUMBEROFFRAMES Extract the number of frames.
% SDLFMXSDLFMKERNCOMPUTE Compute a cross kernel between two SDLFM kernels.
% INVGAMMAPRIOREXTRACTPARAM Extract params from inverse gamma prior structure.
% GPTIMEDYNAMICSCREATE Create the time dynamics model. 
% GPREVERSIBLEDYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.
% HEATKERNCOMPUTE Compute a kernel matrix for a HEAT kernel.
% PROBITNOISEOUT Compute the output of the PROBIT noise given the input mean and variance.
% LFMWHITEXLFMWHITEKERNGRADIENT Compute a cross gradient between two
% NGAUSSNOISEOUT Compute the output of the NGAUSS noise given the input mean and variance.
% GIBBSPERIODICKERNPARAMINIT GIBBSPERIODIC kernel parameter initialisation.
% NGAUSSNOISEPOINTPLOT Plot the data-points for the NGAUSS noise model.
% LVMTWODPLOT Helper function for plotting the labels in 2-D.
% MODELREADFROMFILE Read model from a file FID produced by the C++ implementation.
% MVUDECONSTRUCT break MVU in pieces for saving.
% MVUOPTIONS Options for a MVU.
% MATERN52KERNDIAGCOMPUTE Compute diagonal of MATERN52 kernel.
% XYZANKURANIM Animate point cloud of stick man from Agarwal & Triggs dataset.
% XYZHUMANEVAVISUALISE2
% MLPLOGLIKEGRADIENTS Multi-layer perceptron gradients.
% LOGISTICNORMALPRIOREXPANDPARAM Expand logistic-normal prior structure from params.
% CMPNDNOISENUG  Update nu and g parameters associated with compound noise model.
% MGAUSSIANNOISEDISPLAY Display parameters of the MGAUSSIAN noise.
% LMCKERNGRADX Gradient of LMC kernel with respect to input locations.
% DEMVOWELSFGPLVM3 Model the vowels data with a 2-D FGPLVM using RBF kernel and back constraints, but without PCA initialisation.
% NCNMNOISELOGLIKELIHOOD Log-likelihood of data under null category noise model.
% SIMXSIMKERNDIAGGRADIENT Gradient for the diagonal between two SIM kernels.
% PREPAREPLOT Helper function for tidying up the plot before printing.
% DEMROBOTWIRELESSFGPLVM3 Wireless Robot data from University of Washington with dynamics and no back constraints.
% NONEKERNGRADX Gradient of NONE kernel with respect to a point x.
% GAUSSIANNOISEEXPANDPARAM Create noise structure from GAUSSIAN noise's parameters.
% XYZHUMANEVA2JOINT
% GAUSSIANKERNDIAGGRADIENT Compute the gradient of the gaussian kernel's diagonal wrt parameters.
% FGPLVMOPTIMISEPOINT Optimise the postion of a latent point.
% KBROUT Compute the output of a KBR model given the structure and input X.
% SKELREVERSELOOKUP Return the number associated with the joint name.
% GAUSSIANKERNGRADIENT Gradient of gaussian kernel's parameters.
% LFMVXLFMVKERNCOMPUTE Velocity and velocity LFM kernel  
% MATERN52KERNDIAGGRADIENT Compute the gradient of the MATERN52 kernel's diagonal wrt parameters.
% RBFHKERNPARAMINIT RBFH kernel parameter initialisation.
% IVMEPLOGLIKELIHOOD Return the EP approximation to the log-likelihood.
% LFMAXLFMVKERNGRADIENT Compute a cross gradient between a LFMA and a LFMV.
% DEMSTICKFGPLVM3 Model the stick man using an RBF kernel and RBF kernel based back constraints.
% MODELWRITETOFID Write to a stream a given model.
% MVUCREATE Maximum variance unfolding embedding model.
% INVGAMMAPRIORGRADIENT Gradient wrt x of the log Gaussian prior.
% GPDYNAMICSSAMP Sample from the dynamics for a given input.
% SDLFMAXSDLFMVKERNGRADIENTBLOCK Gradients of the parameters in block i,j
% PLOTMATRIX Fill a given axis with a matrix plot.
% LFMCOMPUTEH4JA Helper function for computing part of the LFMJA kernel.
% XYZVISUALISE For drawing an xyz representation of 3-D data.
% GPOBJECTIVE Wrapper function for GP objective.
% GPDYNAMICSLOGLIKELIHOOD Give the log likelihood of GP dynamics.
% RATQUADKERNEXPANDPARAM Create kernel structure from RATQUAD kernel's parameters.
% FGPLVMOBJECTIVE Wrapper function for GP-LVM objective.
% LFMAXRBFKERNCOMPUTE Compute cross kernel between the LFMA and RBF kernels.
% DATASETSDIRECTORY Returns directory where data is stored.
% LOGISTICNORMALPRIORLOGPROB Log probability of logistic-normal prior.
% LFMTEST Test the gradients of the LFM model.
% MULTIKERNDIAGGRADIENTBLOCK
% GGWHITEXGGWHITEKERNGRADIENT Compute a cross gradient between two GG WHITE kernels.
% LVMSETPLOTNOVAR A copy of lvmSetPlot where the variance in the input
% LFMWHITEKERNPARAMINIT LFM-WHITE kernel parameter initialisation.
% NDDISIMKERNPARAMINIT NDDISIM kernel parameter initialisation.
% IVMOPTIONS Return default options for IVM model.
% LINEARLOGLIKEGRADIENTS Linear model gradients.
% RBFXNONEKERNGRADIENT Compute a cross gradient between RBF and DUMMY
% RBFARDKERNDIAGGRADIENT Compute the gradient of the RBFARD kernel's diagonal wrt parameters.
% NOISEDISPLAY Display the parameters of the noise model.
% FGPLVMRECONSTRUCT Reconstruct an FGPLVM from component parts.
% MULTIKERNCACHEBLOCK
% HEATKERNDIAGGRADIENT Gradient of the HEAT kernel's diagonal wrt parameters.
% FGPLVMSEQUENCEOBJECTIVE Wrapper function for objective of a single sequence in latent space and the corresponding output sequence.
% SDLFMXSDLFMKERNGRADIENTBLOCK Gradients of the parameters in block i,j
% MOGLOWERBOUND Computes lower bound on log likelihood for an MOG model.
% DEMCMU35GPLVMFGPLVM1 Learn a GPLVM on CMU 35 data set.
% GPSUBSPACECREATE 
% DISIMXSIMKERNGRADIENT Compute gradient between the DISIM and SIM kernels.
% GAUSSIANWHITEKERNGRADX Gradient of gaussian white kernel with respect 
% TOKENISE Split a string into separate tokens.
% MLPKERNEXPANDPARAM Create kernel structure from MLP kernel's parameters.
% SDLFMVXSDRBFKERNCOMPUTEBLOCK Cross kernel between SDLFM and SDRBF for i,j
% POLYKERNPARAMINIT POLY kernel parameter initialisation.
% LFMWHITEKERNGRADX Gradient of LFM-WHITE kernel with respect to a point t.
% POLYARDKERNDISPLAY Display parameters of the POLYARD kernel.
% DNETTEST Test some settings for the density network.
% SDLFMKERNDIAGCOMPUTEBLOCK Diagonal of a SDLFM kernel matrix for block i
% WHITEKERNGRADX Gradient of WHITE kernel with respect to input locations.
% DEMTHREEFIVEIVM1 Try the IVM & NCNM on 3 vs 5.
% PRIOREXTRACTPARAM Extract the prior model's parameters.
% LFMVPGRADIENTUPSILONMATRIX Gradient upsilon matrix vel. pos.
% MATERN32KERNGRADIENT Gradient of MATERN32 kernel's parameters.
% ORDEREDNOISEUPDATEPARAMS Update parameters for ordered categorical noise model.
% SIMWHITEXRBFINFWHITEKERNCOMPUTE Compute a cross kernel between a SIM-WHITE
% IVMREMOVEPOINT Removes a given point from the IVM.
% SIMWHITEXWHITEKERNCOMPUTE Compute a cross kernel between the SIM-WHITE
% LOGISTICNORMALPRIORPARAMINIT Logistic-normal prior model's parameter initialisation.
% LFMGRADIENTSIGMAH4 Gradient of the function h_i(z) with respect \sigma.
% CMPNDKERNDISPLAY Display parameters of the CMPND kernel.
% CMPNDKERNEXTRACTPARAMTRANSFORMSETTINGS Extract parameter transform settings 
% NDDISIMXNDSIMKERNCOMPUTE Compute a cross kernel between DISIM and SIM kernels with no decay in the SIM part.
% GAUSSIANKERNGRADX Gradient of gaussian kernel with respect to input locations.
% VELOTRANSKERNPARAMINIT VELOTRANS kernel parameter initialisation.
% NOISEWRITEPARAMSTOFID Write the noise parameters to a stream.
% GPTIMEDYNAMICSSETLATENTVALUES Set the latent values inside the model.
% INVCMPNDKERNCOMPUTE Compute the INVERSE-PRECISION-CMPND kernel given the parameters and X.
% NOISETEST Run some tests on the specified noise model.
% OPTIMISEPARAMS Optimise parameters.
% SDLFMAXSDRBFKERNGRADIENTBLOCK Gradients of the parameters in block i,j
% SIMKERNGRADX Gradient of SIM kernel with respect to each time point in t1.
% WANGPRIORPARAMINIT Wang prior model's parameter initialisation.
% HEATKERNEXPANDPARAM Create kernel structure from HEAT kernel's parameters.
% IVMPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% LFMVXRBFKERNCOMPUTE Compute a cross kernel between the LFMV and RBF kernels.
% MULTIMODELPARAMINIT MULTIMODEL model parameter initialisation.
% MODELLOGLIKELIHOOD Compute a model log likelihood.
% GIBBSKERNEXTRACTPARAM Extract parameters from the GIBBS kernel structure.
% WHITEKERNPARAMINIT WHITE kernel parameter initialisation.
% MOGLOGLIKELIHOOD Mixture of Gaussian's log likelihood.
% LVMNEARESTNEIGHBOUR Give the number of errors in latent space for 1 nearest neighbour.
% SDLFMAXSDLFMAKERNGRADIENT Gradients of cross kernel between 2 SDLFM kernels.
% NGAUSSNOISEEXTRACTPARAM Extract parameters from the NGAUSS noise structure.
% LINEAROPTIONS Options for learning a linear model.
% NONEKERNDISPLAY Display parameters of the NONE kernel.
% READSTRINGFROMFID Read an boolean from an FID.
% CMPNDNOISEEXTRACTPARAM Extract parameters from the CMPND noise structure.
% MULTIKERNFIXBLOCKS
% KBRDISPLAY Display parameters of the KBR model.
% POLYARDKERNEXPANDPARAM Create kernel structure from POLYARD kernel's parameters.
% MODELEXPANDPARAM Update a model structure with parameters.
% XYZPOPPEMODIFY
% LFMGRADIENTSIGMAH3AP Gradient of the function h_i(z) with respect \sigma.
% UNIFORMPRIOREXPANDPARAM Expand uniform prior structure from params.
% ARDKERNGRADX Gradient of ARD kernel with respect to a point x.
% NCNMNOISEGRADIENTPARAM Gradient of parameters for NCNM.
% LVMCLICKVISUALISE Visualise the manifold using clicks.
% SIMCOMPUTEHSTAT Helper function for computing part of the stationary version
% NDSIMKERNGRADIENT Gradient of SIM kernel's parameters.
% TENSORKERNDIAGCOMPUTE Compute diagonal of TENSOR kernel.
% OPTIMIDEFAULTOPTIMISER Returns the default optimiser to be used.
% POLYARDKERNGRADX Gradient of POLYARD kernel with respect to input locations.
% DEMBRENDANFGPLVM3 Use the GP-LVM to model the Frey face data with DTCVAR.
% FGPLVMDYNAMICSPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% DEMOILFGPLVM4 Oil data with deterministic training conditional, and MLP back constraints.
% FGPLVMDYNAMICSRUN Runs auto regressive dynamics in a forward manner.
% WHITEKERNCOMPUTE Compute the white-noise (WHITE) kernel between
% LFMVKERNEXPANDPARAM Create kernel structure from LFMV kernel's parameters.
% LFMVPGRADIENTUPSILONVECTOR Gradient upsilon vector vel. pos.
% FGPLVMPOINTLOGLIKELIHOOD Log-likelihood of a point for the GP-LVM.
% IMAGEMODIFY Helper code for visualisation of image data.
% FGPLVMSEQUENCELOGLIKEGRADIENT Log-likelihood gradient for of a sequence of the GP-LVM.
% RBFINFWHITEKERNDISPLAY Display parameters of the RBF-WHITE kernel (with
% SIMXSIMKERNCOMPUTE Compute a cross kernel between two SIM kernels.
% IVMRECONSTRUCT Reconstruct an IVM form component parts.
% CMPNDKERNEXTRACTPARAM Extract parameters from the CMPND kernel structure.
% RBFHKERNDIAGCOMPUTE Compute diagonal of RBFH kernel.
% DEMUNLABELLEDONEIVM1 Test IVM code on a toy crescent data.
% DIAGKERNDIAGCOMPUTE Compute diagonal of DIAG kernel.
% IVMPOSTERIORGRADMEANVAR Gradient of mean and variances of the posterior wrt X.
% SDLFMXSDLFMKERNGRADIENTICBLOCK Partial derivatives initial conditions
% DEMSILHOUETTEAVERAGE Shows the average of the poses.
% LFMAPGRADIENTUPSILONMATRIX Gradient upsilon matrix accel. pos.
% GRADIENTCHECK Check gradients of objective function.
% MOCAPLOADTEXTDATA Load a motion capture data set from a text file.
% DOUBLEMATRIXWRITETOFID Writes a double matrix to an FID.
% HEATXHEATKERNCOMPUTE Compute a cross kernel between two HEAT kernels.
% INDEXARDKERNEXTRACTPARAM Extract parameters from the INDEXARD kernel structure.
% MULTIKERNGRADX Gradient of MULTI kernel with respect to a point x.
% LFMAXRBFVKERNCOMPUTE Compute cross kernel between the LFMA and RBFV kernels.
% PPCAOPTIONS Options for probabilistic PCA.
% INVCMPNDKERNDIAGGRADIENT Compute the gradient of the INVCMPND kernel's diagonal wrt parameters.
% NONEKERNCOMPUTE Compute the NONE kernel given the parameters and X.
% GPDYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.
% MODELLOADRESULT Load a previously saved result.
% GETLOCALGRADALPHAOMEGA Gradients of parameters in alpha and omega
% DEMGPTWOSAMPLE Test GP two sample code.
% LFMVVCOMPUTEUPSILONMATRIX Upsilon matrix vel. vel. with t1, t2 limits
% SIMWHITEXSIMWHITEKERNGRADIENT Compute a cross gradient between two
% GPWRITERESULT Write a GP result.
% RBFHKERNDISPLAY Display parameters of the RBFH kernel.
% GPDYNAMICSLOGLIKEGRADIENTS Gradients of the GP dynamics wrt parameters.
% NDSIMKERNDISPLAY Display parameters of the NDSIM kernel.
% GAUSSIANWHITEKERNGRADIENT Gradient of gaussian white kernel's parameters.
% DEMSWISSROLLLLE5 Demonstrate LLE on the oil data.
% SDLFMXSDRBFKERNGRADIENTBLOCKIGJ 
% XYZHUMANEVAVISUALISEMODES
% DEMVOWELSFGPLVM2 Model the vowels data with a 2-D FGPLVM using RBF kernel.
% LFMWHITECOMPUTEH Helper function for computing part of the LFM-WHITE
% LOGISTICNORMALPRIORSETBOUNDS Set logistic-normal prior bounds.
% XYZHUMANEVAVISUALISE3D
% PRINTLATEXPLOT Print a plot to LaTeX.
% FILEKERNGRADX Gradient of FILE kernel with respect to a point x.
% SDLFMVXSDLFMKERNGRADIENTBLOCK Gradients of the parameters in block i,j
% IVMSELECTPOINT Choose a point for inclusion or removal.
% DEXPKERNGRADX Gradient of the double exponential kernel with respect to a
% XYZHUMANEVAMODIFY2
% INVGAMMAPRIOREXPANDPARAM Expand inverse gamma prior structure from params.
% LVMPRINTPLOT Print latent space for learnt model.
% LFMVVGRADIENTUPSILONMATRIX Gradient upsilon matrix vel. vel.
% NORMUNIPRIOREXTRACTPARAM Extract params from normal uniform prior structure.
% OUKERNDIAGCOMPUTE Compute diagonal of OU kernel (see ouKernCompute or
% SIMWHITEKERNDIAGGRADX Gradient of SIM-WHITE kernel's diagonal w.r.t. t.
% LINKERNGRADIENT Gradient of LIN kernel's parameters.
% NDSIMKERNEXPANDPARAMTRANSFORMSETTINGS Create kernel structure from SIM kernel's parameters' transform settings.
% PPCAPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% LFMCOMPUTEUPSILON Helper function for comptuing part of the LFM kernel.
% NORMUNIPRIORPARAMINIT Normal uniform prior model's parameter initialisation.
% GAUSSIANNOISEGRADIENTPARAM Gradient of GAUSSIAN noise's parameters.
% NCNMNOISEDISPLAY Display  parameters from null category noise model.
% DEMSWISSROLLLLE1 Demonstrate LLE on the oil data.
% ROCHOLFACTORISE Rank one Cholesky factorise.
% LFMAKERNCOMPUTE Compute the LFMA kernel given the parameters and X.
% RBFARDKERNCOMPUTE Compute the RBFARD kernel given the parameters and X.
% XYZHUMANEVAMODIFY
% LFMJPCOMPUTEUPSILONMATRIX Upsilon matrix jolt. pos. with t1, t2 limits
% SDLFMVKERNPARAMINIT SDLFMV kernel initialization
% LNDIFFCUMGAUSSIAN Log of the difference between two cumulative Gaussians.
% ORDEREDNOISE3DPLOT Draws a 3D or contour plot for the ORDERED noise model.
% GAUSSIANNOISELIKELIHOOD Likelihood of the data under the GAUSSIAN noise model.
% DEMOILFGPLVM1 Oil data with fully independent training conditional.
% RBFWHITEKERNEXTRACTPARAM Extract parameters from the RBF-WHITE kernel
% SDLFMVKERNDIAGCOMPUTE Compute diagonal of a SDLFMV kernel.
% NOISEGRADVALS Gradient of noise model wrt mu and varsigma.
% IVMDOWNDATENUG Downdate nu and g parameters associated with noise model.
% DISIMCOMPUTEH Helper function for comptuing part of the DISIM kernel.
% MATERN32KERNEXTRACTPARAM Extract parameters from the MATERN32 kernel structure.
% WHITEKERNEXPANDPARAMTRANSFORMSETTINGS Create kernel structure from WHITE kernel's
% GAUSSIANWHITEKERNDIAGGRADIENT Compute the gradient of the gaussian white 
% GPDYNAMICSPOINTLOGLIKELIHOOD Compute the log likelihood of a point under the GP dynamics model.
% LFMKERNGRADX Gradient of LFM kernel with respect to a point x.
% NONEKERNEXTRACTPARAM Extract parameters from the NONE kernel structure.
% LFMCOMPUTEH4JV Helper function for computing part of the LFMJV kernel.
% EXPKERNCOMPUTE Compute the EXP kernel given the parameters and X.
% RBFPERIODIC2KERNGRADIENT Gradient of RBFPERIODIC2 kernel's parameters.
% WHITEKERNDIAGCOMPUTE Compute diagonal of WHITE kernel.
% GAUSSIANPRIORPARAMINIT Gaussian prior model's parameter initialisation.
% TENSORKERNGRADIENT Gradient of TENSOR kernel's parameters.
% SDLFMAXSDLFMVKERNGRADIENT Gradients of cross kernel between 2 SDLFM kernels.
% SDLFMKERNDIAGCOMPUTE Compute diagonal of a SDLFM kernel.
% PROBITNOISELOGLIKELIHOOD Log likelihood of the data under the PROBIT noise model.
% PROBITNOISEEXTRACTPARAM Extract parameters from the PROBIT noise structure.
% GPUPDATEKERNELS Update the kernels that are needed.
% READBOOLFROMFID Read a boolean from an FID.
% MLPKERNPARAMINIT MLP kernel parameter initialisation.
% LLEOPTIMISE Optimise an LLE model.
% LFMAXLFMAKERNGRADIENT Compute a cross gradient between a LFMA and a LFMA.
% GGWHITEKERNCOMPUTE Compute the GG white kernel given the parameters and X.
% GPTIMEDYNAMICSEXPANDPARAM Place the parameters vector into the model for GP time dynamics.
% LFMAPCOMPUTEUPSILONMATRIX Upsilon matrix acce. pos. with t1, t2 limits
% NOISEREADPARAMSFROMFID Read the noise parameters from C++ file FID.
% LFMVPGRADIENTSIGMAUPSILONVECTOR Gradient of upsilon vector vp wrt sigma
% SDLFMAXSDLFMKERNGRADIENTBLOCKILJ 
% DEMINTERPOLATIONGP Demonstrate Gaussian processes for interpolation.
% SKELPLAYDATA Play skel motion capture data.
% KERNDIAGGRADX Compute the gradient of the  kernel wrt X.
% LFMGRADIENTSIGMAH3 Gradient of the function h_i(z) with respect \sigma.
% GPREVERSIBLEDYNAMICSSAMP Sample from the dynamics for a given input.
% KERNEXPANDPARAM Expand parameters to form a kernel structure.
% BVHREADFILE Reads a bvh file into a tree structure.
% SKELCONNECTIONMATRIX Compute the connection matrix for the structure.
% FILEKERNPARAMINIT FILE kernel parameter initialisation.
% VIVMRUNDATASETLEARN Try the virtual IVM on a data set and save the results.
% MODELLOGLIKEGRADIENTS Compute a model's gradients wrt log likelihood.
% LFMKERNEXTRACTPARAM Extract parameters from the LFM kernel structure.
% IVMSELECTPOINTS Selects the point for an IVM.
% GPSUBSPACEEXTRACTPARAM
% KERNSETINDEX Set the indices on a compound kernel.
% FGPLVMLOGLIKELIHOOD Log-likelihood for a GP-LVM.
% PRIORCREATE Create a prior structure given a type.
% LINARDKERNEXPANDPARAM Create kernel structure from LINARD kernel's parameters.
% SCALENOISEOUT A simple noise model that scales and centres the data.
% FGPLVMDISPLAY Display an FGPLVM model.
% MODELSETLATENTVALUES Set the latent variables for dynamics models in the GPLVM.
% PRIORSETBOUNDS Set the bounded prior model's bounds from bounds vector.
% GPTIMEDYNAMICSEXTRACTPARAM Extract parameters from the GP time dynamics model.
% LFMKERNPARAMINIT LFM kernel parameter initialisation. The latent force
% GPTIMEDYNAMICSLOGLIKELIHOOD Give the log likelihood of GP time dynamics.
% RBFPERIODIC2KERNGRADX Gradient of RBFPERIODIC2 kernel with respect to a point x.
% RBFINFWHITEKERNEXPANDPARAM Create kernel structure from RBF-WHITE kernel's
% TABLEREAD Read in data which has column titles in the first line and separated values in each other line.
% RBFCREATE Wrapper for NETLAB's rbf `net'.
% WHITEFIXEDKERNEXPANDPARAM Create kernel structure from WHITEFIXED kernel's parameters.
% ROBTHREEDYNAMICSEXPANDPARAM Place the parameters vector into the model for first robot dynamics.
% LFMAKERNEXPANDPARAM Create kernel structure from LFMA kernel's parameters.
% LVMSCATTERPLOTCOLOR 2-D scatter plot of the latent points with color.
% LINKERNEXPANDPARAM Create kernel structure from LIN kernel's parameters.
% PRIORWRITEPARAMSTOFID Write prior params from C++ written FID.
% MODELPOINTLOGLIKELIHOOD Compute the log likelihood of a given point.
% SDLFMVXSDLFMVKERNCOMPUTE Compute a cross kernel between two SDLFMV kernels.
% MULTIMODELLOGLIKELIHOOD Log likelihood of MULTIMODEL model.
% XYZHUMANEVADRAW
% SIMKERNDIAGGRADX Gradient of SIM kernel's diagonal with respect to the
% CMPNDKERNPARAMINIT CMPND kernel parameter initialisation.
% XYZHUMANEVAVISUALISE
% RBFKERNDISPLAY Display parameters of the RBF kernel.
% MATERN52KERNDIAGGRADX Gradient of MATERN52 kernel's diagonal with respect to X.
% SDLFMXSDLFMKERNGRADIENTBLOCKIEJ 
% MATERN32KERNGRADX Gradient of MATERN32 kernel with respect to input locations.
% SDLFMAMEANGRADIENT Gradients wrt parameters of the accel. mean SDLFM.
% DEMCMU35TAYLORNEARESTNEIGHBOUR Recreate the Nearest Neighbour result from Taylor et al, NIPS 2006.
% POLYARDKERNPARAMINIT POLYARD kernel parameter initialisation.
% XYZHUMANEVAANIM
% LFMAXLFMVKERNCOMPUTE Acceleration and velocity LFM kernel  
% LVMSCOREMODEL Score model with a GP log likelihood.
% SCALENOISESITES Site updates for Scale noise model.
% HEATKERNDIAGCOMPUTE Diagonal of a kernel matrix for a HEAT kernel.
% ISOMAPEMBED Embed data set with Isomap.
% XYZPOPPEVISUALISE Draw the Poppe figure return the graphics handle.
% RBFHKERNDIAGGRADX Gradient of RBFH kernel's diagonal with respect to X.
% POLYARDKERNCOMPUTE Compute the POLYARD kernel given the parameters and X.
% INDEXKERNDIAGCOMPUTE Compute diagonal of INDEX kernel.
% DEMBRENDANFGPLVM4 Use the GP-LVM to model the Frey face data with DTCVAR and back constraints.
% INVCUMGAUSSIAN Computes inverse of the cumulative Gaussian.
% GGWHITEKERNGRADIENT Gradient of GG WHITE kernel's parameters.
% KERNEXPANDPARAMTRANSFORMSETTINGS Expand parameters' transform settings to form a kernel structure.
% SCATTERPLOT 2-D scatter plot of labelled points.
% LFMCLASSVISUALISE Callback function to visualize LFM in 2D
% SIGMOIDTRANSFORM Constrains a parameter to be between 0 and 1.
% MATERN52KERNEXTRACTPARAM Extract parameters from the MATERN52 kernel structure.
% INVCMPNDKERNDIAGCOMPUTE Compute diagonal of INVCMPND kernel.
% VELOTRANSKERNDIAGGRADX Gradient of VELOTRANS kernel's diagonal with respect to X.
% NGAUSSNOISEGRADVALS Gradient of NGAUSS noise log Z with respect to input mean and variance.
% IVMGUNNARDATA Script for running experiments on Gunnar data.
% IVMNEGGRADIENTNOISE Wrapper function for calling noise param gradients.
% DEMTHREEFIVERESULTS Plot results from the three vs five experiments.
% GPTIMEDYNAMICSDISPLAY Display a GP time dynamics model.
% RBFHKERNEXPANDPARAM Create kernel structure from RBFH kernel's parameters.
% RBFINFWHITEXRBFINFWHITEKERNCOMPUTE Compute a cross kernel between two
% WHITEHKERNEXPANDPARAM Create kernel structure from WHITEH kernel's parameters.
% MODELDISPLAY Display a text output of a model.
% GPPOSTERIORGRADMEANVAR Gadient of the mean and variances of the posterior at points given by X.
% ARDKERNDISPLAY Display parameters of the ARD kernel.
% SDLFMXSDLFMKERNGRADIENTIC Computes partial derivative for init. const.
% LFMCOMPUTEH3AP Helper function for computing part of the LFMAP kernel.
% WHITEFIXEDKERNGRADX Gradient of WHITEFIXED kernel with respect to a point x.
% FGPLVMPOINTLOGLIKEGRADIENT Log-likelihood gradient for of a point of the GP-LVM.
% LFMWHITEKERNGRADIENT Gradient of LFM-WHITE kernel's parameters.
% GPDYNAMICSEXPANDPARAM Place the parameters vector into the model for GP dynamics.
% IVMOPTIMISEIVM Selects the points for an IVM model.
% WHITEHKERNDIAGCOMPUTE Compute diagonal of WHITEH kernel.
% GGKERNGRADIENT Gradient of GG kernel's parameters.
% ROBTHREEDYNAMICSLOGLIKELIHOOD Give the log likelihood of the robot three dynamics part.
% DEMUSPSIVM3 Try the ARD IVM on some digits data.
% NOISEREADFROMFID Load from an FID written by the C++ implementation.
% GAUSSIANNOISELOGLIKELIHOOD Log likelihood of the data under the GAUSSIAN noise model.
% ROBTHREEDYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.
% NCNMNOISEPARAMINIT null category noise model's parameter initialisation.
% DISIMKERNCOMPUTE Compute the DISIM kernel given the parameters and X.
% RAYLEIGHSAMP Sample from a Rayleigh with a given sigma.
% EXPKERNDIAGCOMPUTE Compute diagonal of EXP kernel.
% MATERN32KERNDIAGCOMPUTE Compute diagonal of MATERN32 kernel.
% IVMUSPSRESULTS Summarise the USPS result files in LaTeX.
% MLPLOGLIKELIHOOD Multi-layer perceptron log likelihood.
% LLERECONSTRUCT Reconstruct an LLE form component parts.
% BIASKERNGRADX Gradient of BIAS kernel with respect to input locations.
% RBFPERIODICOUT Compute the output of a RBFPERIODIC model given the structure and input X.
% INVCMPNDKERNGRADX 
% MVUEMBED Embed data set with MVU.
% MATERN32KERNEXPANDPARAM Create kernel structure from MATERN32 kernel's parameters.
% DEMOPTIMISEGPTUTORIAL Shows that there is an optimum for the covariance function length scale.
% LFMGRADIENTSIGMAH3VV Gradient of the function h_i(z) with respect \sigma.
% ORDEREDNOISEGRADIENTPARAM Gradient of ORDERED noise's parameters.
% GPTIMEDYNAMICSLATENTGRADIENTS Gradients of the X vector given the time dynamics model.
% SCALENOISEDISPLAY Display the parameters of the scaled noise model.
% RBFHKERNDIAGGRADIENT Gradient of the RBFH kernel's diagonal wrt parameters.
% DEMUSPSIVM1 Try the IVM on the USPS digits data with RBF kernel.
% INDEXARDKERNPARAMINIT INDEXARD kernel parameter initialisation.
% LVMSETPLOT Sets up the plot for visualization of the latent space.
% SDLFMKERNGRADIENTMEAN Gradients of the parameters of mean function cov.
% KBROPTIONS Create a default options structure for the KBR model.
% TRANSLATEKERNDISPLAY Display parameters of the TRANSLATE kernel.
% INVCMPNDKERNEXTRACTPARAM Extract parameters from the INVCMPND kernel structure.
% SDLFMAXSDLFMAKERNCOMPUTEBLOCK Computes SDLFM kernel matrix for block i,j
% SPARSEDIAG Create a diagonal matrix that is sparse from a vector.
% ISOMAPOPTIONS Options for a isomap.
% KERNREADPARAMSFROMFID Read the kernel parameters from C++ file FID.
% LVMCLASSVISUALISE Callback function for visualising data.
% MAPLOADDATA Load a mapping model dataset (e.g. classification, regression).
% IVMUPDATENUG Update nu and g parameters associated with noise model.
% FGPLVMCOVGRADSTEST Test the gradients of the covariance.
% TENSORKERNCOMPUTE Compute the TENSOR kernel given the parameters and X.
% INDEXARDKERNEXPANDPARAM Create kernel structure from INDEXARD kernel's parameters.
% DEMSILHOUETTELINEAR1 Model silhouette data with independent linear models.
% RBFPERIODIC2KERNPARAMINIT RBFPERIODIC2 kernel parameter initialisation.
% GRADLOGCUMGAUSSIAN Gradient of the log of the cumulative Gaussian.
% XLSLOADDATA Wrapper function for xlsread to get files from the datasets directory.
% DISIMKERNEXPANDPARAM Create kernel structure from DISIM kernel's parameters.
% MGAUSSIANNOISEGRADVALS Gradient of MGAUSSIAN noise log Z with respect to input mean and variance.
% RBFHKERNGRADIENT Gradient of RBFH kernel's parameters.
% PARAMNAMEREGULAREXPRESSIONLOOKUP Returns the indices of the parameter containing the given regular expression.
% TRANSLATEKERNDIAGGRADX Gradient of TRANSLATE kernel's diagonal with respect to X.
% SDLFMVXSDLFMKERNCOMPUTEBLOCK Computes SDLFM kernel matrix for block i,j
% FGPLVMDYNAMICSFIELDPLOT 2-D field plot of the dynamics.
% BVHPLAYFILE Play motion capture data from a bvh file.
% DISIMXRBFKERNGRADIENT Compute gradient between the DISIM and RBF kernels.
% LFMAPGRADIENTSIGMAUPSILONMATRIX Gradient of upsilon matrix ap wrt sigma
% DISIMKERNDIAGCOMPUTE Compute diagonal of DISIM kernel.
% SDLFMVXSDLFMVKERNGRADIENTICBLOCK Partial derivatives initial conditions
% WHITEXWHITEKERNGRADIENT Compute a cross gradient between two WHITE kernels.
% DISIMKERNDIAGGRADX Gradient of DISIM kernel's diagonal with respect to X.
% RBFARDKERNGRADX Gradient of RBFARD kernel with respect to input locations.
% MODELOPTIONS Returns a default options structure for the given model.
% RBFWHITEKERNDISPLAY Display parameters of the RBF-WHITE kernel.
% KERNTEST Run some tests on the specified kernel.
% KERNCREATE Initialise a kernel structure.
% NDDISIMKERNEXTRACTPARAM Extract parameters from the NDDISIM kernel structure.
% GIBBSPERIODICKERNDISPLAY Display parameters of the GIBBSPERIODIC kernel.
% LINARDKERNCOMPUTE Compute the LINARD kernel given the parameters and X.
% RATQUADKERNGRADIENT Gradient of RATQUAD kernel's parameters.
% ARDKERNDIAGGRADIENT Compute the gradient of the ARD kernel's diagonal wrt parameters.
% OPTIMIDEFAULTCONSTRAINT Returns function for parameter constraint.
% DIAGKERNDISPLAY Display parameters of the DIAG kernel.
% PLOT3VISUALISE  Helper code for plotting a plot3 visualisation.
% SIMKERNPARAMINIT SIM kernel parameter initialisation.
% RBFARDJITKERNPARAMINIT RBFARD2 kernel parameter initialisation.
% NCNMNOISEGRADVALS Compute gradient with respect to inputs to noise model.
% TRANSLATEKERNGRADIENT Gradient of TRANSLATE kernel's parameters.
% INDEXKERNEXPANDPARAM Create kernel structure from INDEX kernel's parameters.
% LFMVKERNDISPLAY Display parameters of the LFMV kernel.
% SIMCOMPUTEH Helper function for comptuing part of the SIM kernel.
% WHITEXWHITEKERNCOMPUTE Compute a cross kernel between two WHITE kernels.
% KERNDIAGGRADIENT Compute the gradient of the kernel's parameters for the diagonal.
% GPMEANFUNCTIONGRADIENT Compute the log likelihood gradient wrt the scales.
% DEMOILFGPLVM6 Oil data with partially independent training conditional, and MLP back constraints.
% OUKERNDISPLAY Display parameters of the OU kernel (see ouKernCompute or
% DNETLOWERBOUND Computes lower bound on log likelihood for an DNET model.
% CELL2NUM Converts a cell array of strings to numbers.
% READINTFROMFID Read an integer from an FID.
% SDLFMAXSDRBFKERNCOMPUTEBLOCK Cross kernel between SDLFM and SDRBF for i,j
% DEMSWISSROLLFULLLLE4 Demonstrate LLE on the oil data.
% RBFARD2KERNDIAGGRADIENT Compute the gradient of the RBFARD2 kernel's diagonal wrt parameters.
% RBFKERNGRADIENT Gradient of RBF kernel's parameters.
% LFMWHITEXRBFWHITEKERNCOMPUTE Compute a cross kernel between an LFM-WHITE
% LFMAVGRADIENTSIGMAUPSILONMATRIX Gradient of upsilon matrix av wrt sigma
% GPPOSTERIORVAR Variances of the posterior at points given by X.
% ROTATIONMATRIX Compute the rotation matrix for an angle in each direction.
% SHEATKERNDIAGCOMPUTE Compute a diagonal for the SHEAT kernel matrix.
% LFMCOMPUTEH4 Helper function for computing part of the LFM kernel.
% UNIFORMPRIORLOGPROB Log probability of uniform prior.
% GAUSSIANKERNDIAGGRADX Gradient of gaussian kernel's diagonal with respect to X.
% SDLFMAXSDLFMVKERNCOMPUTEBLOCK Computes SDLFM kernel matrix for block i,j
% ROBONEDYNAMICSEXPANDPARAM Place the parameters vector into the model for first robot dynamics.
% NCNMNOISEPOINTPLOT Plot the data-points for null category noise model.
% TRANSLATEKERNDIAGCOMPUTE Compute diagonal of TRANSLATE kernel.
% MULTIMODELEXPANDPARAM Create model structure from MULTIMODEL model's parameters.
% MOGUPDATECOVARIANCE Update the covariances of an MOG model.
% COMPUTEKERNEL Compute the kernel given the parameters and X.
% DIAGKERNCOMPUTE Compute the DIAG kernel given the parameters and X.
% LFMVPGRADIENTSIGMAUPSILONMATRIX Gradient of upsilon matrix vp wrt sigma
% MAPMODELREADFROMFID Load from a FID produced by C++ code.
% CMPNDNOISE3DPLOT Draws a 3D or contour plot for the CMPND noise model.
% FGPLVMFIELDPLOT 2-D field plot of the dynamics.
% RBFPERIODIC2KERNDIAGCOMPUTE Compute diagonal of RBFPERIODIC2 kernel.
% SIMKERNDIAGGRADIENT Compute the gradient of the SIM kernel's diagonal wrt parameters.
% SDRBFKERNGRADIENT Gradient of SDRBF kernel's parameters.
% XYZANKUR2JOINT Converts data to xyz positions for each joint.
% RBFARDJITKERNDISPLAY Display parameters of the RBFARDJIT kernel.
% SQEXPKERNGRADIENT Gradient of SQEXP kernel's parameters.
% MGAUSSIANNOISEEXTRACTPARAM Extract parameters from the MGAUSSIAN noise structure.
% NEGNOISELOGLIKELIHOOD Wrapper function for calling noise likelihoods.
% DEMBRENDANPPCA1 Use PPCA to model the Frey face data with five latent dimensions.
% GPPOSTERIORSAMPLE Create a plot of samples from a posterior covariance.
% CMPNDNOISELOGLIKELIHOOD Log likelihood of the data under the CMPND noise model.
% PRINTLATEXTEXT Print a text string to file for latex input.
% LFMVVGRADIENTSIGMAUPSILONMATRIX Gradient of upsilon matrix vv wrt sigma
% SDLFMXSDRBFKERNGRADIENTBLOCK Gradients of the parameters in block i,j
% MLPKERNDIAGGRADX Gradient of MLP kernel's diagonal with respect to X.
% RBFARD2KERNGRADIENT Gradient of RBFARD2 kernel's parameters.
% WHITEBLOCKKERNDISPLAY Display parameters of the WHITEBLOCK kernel.
% WIENERKERNDIAGCOMPUTE Compute diagonal of WIENER kernel.
% SIMWHITEKERNEXTRACTPARAM Extract parameters from the SIM-WHITE kernel
% OUKERNPARAMINIT Ornstein-Uhlenbeck (OU) kernel parameter initialisation.
% WHITEFIXEDKERNDIAGGRADX Gradient of WHITEFIXED kernel's diagonal with respect to X.
% SAFESAVE Safe save
% SIMXSIMCOMPUTEDIAGHSTAT Helper function for computing part of the stationary version
% FILEKERNDIAGCOMPUTE Compute diagonal of FILE kernel.
% DEMVOWELSLLE Model the vowels data with a 2-D FGPLVM using RBF kernel.
% DIAGKERNEXTRACTPARAM Extract parameters from the DIAG kernel structure.
% NONEKERNDIAGCOMPUTE Compute diagonal of NONE kernel.
% MODELOBJECTIVE Objective function to minimise for given model.
% GPDECONSTRUCT break GP in pieces for saving.
% GAUSSIANNOISEDISPLAY Display parameters of the GAUSSIAN noise.
% LAPLACEPRIORPARAMINIT Laplace prior model's parameter initialisation.
% FILEKERNEXPANDPARAM Create kernel structure from FILE kernel's parameters.
% SDRBFKERNPARAMINIT SDRBF kernel initialization
% IVMCOVARIANCEGRADIENT The gradient of the likelihood approximation wrt the covariance.
% INVGAMMAPRIORLOGPROB Log probability of inverse gamma prior.
% GPEXPANDPARAM Expand a parameter vector into a GP model.
% LINEAROPTIMISE Optimise a linear model.
% PROBITNOISE3DPLOT Draws a 3D or contour plot for the PROBIT noise model.
% NDSIMKERNDIAGCOMPUTE Compute diagonal of NDSIM kernel.
% PRIOREXPANDPARAM Expand the prior model's parameters from params vector.
% POLYARDKERNGRADIENT Gradient of POLYARD kernel's parameters.
% RBFPERIODICEXPANDPARAM Create model structure from RBFPERIODIC model's parameters.
% ROBONEDYNAMICSLOGLIKEGRADIENTS Gradients of the robot one dynamics wrt parameters.
% LFMWHITECOMPUTEGRADTHETAH1 computes a portion of the LFM-WHITE kernel's gradient w.r.t. theta.
% KERNPCA performs KPCA.
% POLYKERNDIAGGRADX Gradient of POLY kernel's diagonal with respect to X.
% RBFINFWHITEKERNPARAMINIT The RBF-WHITE-INF kernel is a convolutional
% IVMREADFROMFILE Load a file produced by the C++ implementation.
% OUKERNEXTRACTPARAM Extract parameters from the OU kernel structure (see
% WHITEHKERNGRADX Gradient of WHITEH kernel with respect to input locations.
% NCNMNOISESITES Site updates for null category model.
% ROBTHREEDYNAMICSEXTRACTPARAM Extract parameters from the robot three dynamics model.
% GAMMAPRIOREXPANDPARAM Expand gamma prior structure from params.
% GAUSSIANWHITEKERNEXTRACTPARAM Extract parameters from the gaussian white 
% SDLFMVXSDLFMVKERNCOMPUTEBLOCK Computes SDLFM kernel matrix for block i,j
% GIBBSKERNDIAGGRADIENT Compute the gradient of the GIBBS kernel's diagonal wrt parameters.
% SIMXSIMCOMPUTEDIAGH Helper function for comptuing part of the SIM kernel.
% GGWHITEKERNPARAMINIT GG WHITE kernel parameter initialisation.
% RBFOUTPUTGRADX Evaluate derivatives of a RBF model's output with respect to inputs.
% LVMSCATTERPLOT 2-D scatter plot of the latent points.
% LFMCOMPUTEH Helper function for computing part of the LFM kernel.
% MULTIMODELOPTIONS Create a default options structure for the MULTIMODEL model.
% LVMRESULTSCLICK Load a results file and visualise them with clicks
% DEXPKERNDIAGCOMPUTE Compute diagonal of the double exponential kernel.
% RBFPERIODICKERNDIAGCOMPUTE Compute diagonal of RBFPERIODIC kernel.
% PARAMNAMEREVERSELOOKUP Returns the index of the parameter with the given name.
% EXPKERNEXPANDPARAM Create kernel structure from EXP kernel's parameters.
% GPLOGLIKELIHOOD Compute the log likelihood of a GP.
% LINKERNDIAGCOMPUTE Compute diagonal of LIN kernel.
% NOISEGRADX Returns the gradient of the log-likelihood wrt x.
% WHITEXWHITEKERNGRADX
% GPTIMEDYNAMICSSEQUENCELOGLIKELIHOOD Return the log likelihood of a given latent sequence.
% RBFPERIODICKERNCOMPUTE Compute the RBFPERIODIC kernel given the parameters and X.
% MLPEXPANDPARAM Update mlp model with new vector of parameters.
% DEMBRENDANFGPLVM5 Use the GP-LVM to model the Frey face data with DTCVAR and five latent dimensions..
% RBFWHITEXRBFWHITEKERNGRADIENT Compute a cross gradient between two
% SDLFMAXSDLFMKERNGRADIENT Gradients of cross kernel between 2 SDLFM kernels.
% PLOT3MODIFY Helper code for visualisation of 3-d data.
% WHITEHKERNCOMPUTE Compute the WHITEH kernel given the parameters and X.
% KBREXTRACTPARAM Extract parameters from the KBR model structure.
% INDEXKERNEXTRACTPARAM Extract parameters from the INDEX kernel structure.
% RBFWHITEXWHITEKERNCOMPUTE Compute a cross kernel between the RBF-WHITE
% NDSIMKERNPARAMINIT SIM kernel parameter initialisation.
% DNETUPDATEBETA Do an M-step (update parameters) on an Density Network model.
% COMPONENTKERNWRITEPARAMSTOFID Write a component based kernel to a stream.
% SDLFMXSDLFMVKERNGRADIENTICBLOCK Partial derivatives initial conditions
% SDRBFKERNDISPLAY Display parameters of the SDRBF kernel.
% FGPLVMREADFROMFILE Load a file produced by the C++ implementation.
% GAUSSIANWHITEKERNCOMPUTE Compute the covariance of the output samples 
% INDEXARDKERNCOMPUTE Compute the INDEXARD kernel given the parameters and X.
% BIASKERNEXTRACTPARAM Extract parameters from the BIAS kernel structure.
% LFMEXPANDPARAM Expand the given parameters into a LFM structure.
% ISOMAPDECONSTRUCT break isomap in pieces for saving.
% WHITEFIXEDXWHITEFIXEDKERNGRADIENT Compute a cross gradient between two WHITEFIXED kernels.
% RBFARDJITKERNEXTRACTPARAM Extract parameters from the RBFARD2 kernel structure.
% RBFPERIODIC2KERNDIAGGRADX Gradient of RBFPERIODIC2 kernel's diagonal with respect to X.
% CMPNDNOISEDISPLAY Display parameters of the CMPND noise.
% SIGMOIDABTRANSFORM Constrains a parameter to be between A and B
% NOISEGRADIENTPARAM Gradient wrt the noise model's parameters.
% GPPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% MOGOPTIONS Sets the default options structure for MOG models.
% TIMESERIESDATA make a time series data set with the given window length.
% GPDYNAMICSDISPLAY Display a GP dynamics model.
% IVMDOWNDATESITES Downdate site parameters.
% MGAUSSIANNOISEPOINTPLOT Plot the data-points for the MGAUSSIAN noise model.
% DNETOUTPUTGRAD Evaluate derivatives of dnet model outputs with respect to parameters.
% DEMSPGP1DGP5 Do a simple 1-D regression after Snelson & Ghahramani's example.
% MATERN52KERNPARAMINIT MATERN52 kernel parameter initialisation.
% FINDACYCLICNEIGHBOURS2 find the k nearest neighbours for each point in Y preventing cycles in the graph.
% TREEFINDROOTS Return indices of all root nodes in a tree structure.
% MLPARDKERNGRADX Gradient of MLPARD kernel with respect to input locations.
% IVMKERNELOBJECTIVE Compute the negative of the IVM log likelihood approximation.
% LFMKERNGRADIENT Gradient of LFM kernel's parameters.
% LFMVPCOMPUTEUPSILONDIAGVECTOR Upsilon diag vector vel. pos. with t1, t2 limits
% CUMGAMMA Cumulative distribution for gamma.
% WRITEVERSIONTOFID Writes a version to an FID.
% MOCAPCONNECTIONS Give a connection matrix for the motion capture data.
% MATERN52KERNDISPLAY Display parameters of the MATERN52 kernel.
% VECTORVISUALISE  Helper code for plotting a vector during 2-D visualisation.
% DISIMKERNPARAMINIT DISIM kernel parameter initialisation.
% DNETPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% SDLFMAXSDLFMAKERNCOMPUTE Cross kernel between a SDLFMA and a SDLFMA kernels.
% MODELOUTPUTGRAD Compute derivatives with respect to params of model outputs.
% ISOMAPOPTIMISE Optimise an ISOMAP model.
% HEATXRBFHKERNCOMPUTE Cross kernel between a HEAT and a RBF kernels.
% GPTWOSAMPLE Do Oliver Stegles simple two sample test on a data set.
% DEXPKERNPARAMINIT The double exponential kernel is usually called
% GPDATAINDICES Return indices of present data.
% WHITEBLOCKKERNGRADIENT Gradient of WHITEBLOCK kernel's parameters.
% NCNMNOISEOUT Ouput from null category noise model.
% FGPLVMWRITERESULT Write a FGPLVM result.
% FGPLVMGRADIENT GP-LVM gradient wrapper.
% TREESWAPNODE Swap two nodes in the tree structure array.
% PRIORREADFROMFID Read a prior from a C++ written FID.
% LFMGRADIENTSIGMAH4VP Gradient of the function h_i(z) with respect \sigma.
% GPSUBSPACEOPTIMISE
% MATERN52KERNCOMPUTE Compute the MATERN52 kernel given the parameters and X.
% SIMXRBFKERNGRADIENT Compute gradient between the SIM and RBF kernels.
% INDEXARDKERNDISPLAY Display parameters of the INDEXARD kernel.
% OPTOPTIONS Give optimisation options for NETLAB.
% DNETCREATE Density network model.
% FILEKERNCOMPUTE Compute the FILE kernel given the parameters and X.
% PRIORGRADIENT Gradient of the prior with respect to its variables
% LINARDKERNGRADX Gradient of LINARD kernel with respect to input locations.
% CMPNDNOISEPARAMINIT CMPND noise parameter initialisation.
% HEATKERNEXTRACTPARAM Extract parameters from the HEAT kernel structure.
% SIMWHITEXSIMWHITEKERNCOMPUTE Compute a cross kernel between two SIM-WHITE
% PROBITNOISEPARAMINIT PROBIT noise parameter initialisation.
% READVERSIONFROMFID Read version number from an FID.
% MODELPARAMINIT Initialise the parameters of the model.
% JITCHOL Do a Cholesky decomposition with jitter.
% LFMWHITEKERNEXTRACTPARAM Extract parameters from the LFM-WHITE kernel
% IVMAPPROXLOGLIKELIHOOD Return the approximate log-likelihood for the IVM.
% MLPOUTPUTGRADX Evaluate derivatives of mlp model outputs with respect to inputs.
% LMVUEMBED Embed data set with landmark MVU
% PPCACREATE Density network model.
% XYZANKURANIMCOMPARE Animate a prediction and ground truth for stick man from Agarwal & Triggs dataset.
% GGWHITEXWHITEKERNCOMPUTE Compute a cross kernel between a GG white and
% SDLFMVXSDRBFKERNCOMPUTE Cross kernel between a SDLFMV and a SDRBF kernels.
% KBRCREATE Create a KBR model.
% GAUSSOVERDIFFCUMGAUSSIAN A Gaussian over difference of cumulative Gaussians.
% ROBTWODYNAMICSSETLATENTVALUES Set the latent values inside the model.
% GIBBSPERIODICKERNDIAGCOMPUTE Compute diagonal of GIBBSPERIODIC kernel.
% DEMGPTWOSAMPLELIFSH Run GP two sample code on LifSh.
% NEGNOISEGRADIENTPARAM Wrapper function for calling noise gradients.
% RBFWHITEXRBFWHITEKERNCOMPUTE Compute a cross kernel between two RBF-WHITE
% RBFARDJITKERNGRADIENT Gradient of RBFARD2 kernel's parameters.
% CMPNDKERNEXPANDPARAM Create kernel structure from CMPND kernel's parameters.
% HTKLOADMMF File for loading synthesis data from HTK files.
% EXPKERNEXTRACTPARAM Extract parameters from the EXP kernel structure.
% LFMCOMPUTEH4VV Helper function for computing part of the LFMVXLFMV kernel.
% TRANSLATEKERNCOMPUTE Compute the TRANSLATE kernel given the parameters and X.
% MATERN52KERNEXPANDPARAM Create kernel structure from MATERN52 kernel's parameters.
% MATERN52KERNGRADIENT Gradient of MATERN52 kernel's parameters.
% SDLFMVXSDLFMVKERNGRADIENTBLOCK Gradients of the parameters in block i,j
% ARDKERNEXPANDPARAM Create kernel structure from ARD kernel's parameters.
% MODELHESSIAN Hessian of error function to minimise for given model.
% HEATXHEATKERNGRADIENT Gradient wrt parameters between two HEAT kernels.
% LNCUMGAUSSIAN log cumulative distribution for the normalised Gaussian.
% VELOTRANSKERNDIAGCOMPUTE Compute diagonal of VELOTRANS kernel.
% ROBTWODYNAMICSDISPLAY Display the robot dynamics model. 
% RBFOPTIMISE Optimise RBF for given inputs and outputs.
% GRADFUNCWRAPPER Wrapper function to enable use of Carl Rasmussen's minimze function.
% NOISEUPDATESITES Update site parameters for a given noise model.
% BIASKERNDIAGCOMPUTE Compute diagonal of BIAS kernel.
% GGWHITEXGAUSSIANWHITEKERNCOMPUTE Compute a cross kernel between the GG white and GAUSSIAN white kernels.
% WHITEBLOCKKERNCOMPUTE Compute the WHITEBLOCK kernel. 
% LFMGRADIENTH32 Gradient of the function h_i(z) with respect to some of the
% FGPLVMOPTIONS Return default options for FGPLVM model.
% TENSORKERNSETINDEX Set the indices in the tensor kernel.
% ROCHOLHFACTORISE Rank one Cholesky factorise.
% VELOTRANSKERNGRADIENT Gradient of VELOTRANS kernel's parameters.
% SDLFMVKERNCOMPUTE Compute the SDLFMV kernel given the parameters and X.
% ORDEREDNOISEEXPANDPARAM Create noise structure from ORDERED noise's parameters.
% CGCARL Wrapper for Carl Rasmussen's conjugate gradient implemntation.
% VELOTRANSKERNGRADX Gradient of VELOTRANS kernel with respect to a point x.
% SDLFMXSDLFMVKERNGRADIENTBLOCK Gradients of the parameters in block i,j
% MULTIKERNGRADIENTBLOCK
% WHITEFIXEDKERNDISPLAY Display parameters of the WHITEFIXED kernel.
% RBFINFWHITEKERNDIAGGRADX Gradient of RBF-WHITE kernel's (with integration
% LINKERNGRADX Gradient of LIN kernel with respect to input locations.
% NONEKERNGRADIENT Gradient of NONE kernel's parameters.
% GPLVMCMU35ANIMATE Animate the test data jointly with predictions.
% LFMVVCOMPUTEUPSILONDIAGVECTOR Upsilon vector vel. vel. with t1 = t2
% RBFINFWHITEKERNDIAGCOMPUTE Compute diagonal of RBF-WHITE kernel (with
% ORDEREDNOISEDISPLAY Display parameters of the ORDERED noise.
% WIENERKERNGRADX Gradient of WIENER kernel with respect to a point x.
% LFMCOMPUTETEST Test the file lfmComputeH.
% DNETEXPANDPARAM Update dnet model with new vector of parameters.
% LFMVPCOMPUTEUPSILONVECTOR Upsilon vector for vel. pos. with t1 limit
% PROBITNOISEDISPLAY Display parameters of the PROBIT noise.
% UNIFORMPRIOREXTRACTPARAM Extract params from uniform prior structure.
% NEGLOGLOGITTRANSFORM Constrains a parameter to be positive.
% GPTEST Test the gradients of the gpCovGrads function and the gp models.
% PSKERNELGRADIENT Gradient on likelihood approximation for point set IVM.
% MLPOUT Output of an MLP model.
% GIBBSPERIODICKERNGRADX Gradient of GIBBSPERIODIC kernel with respect to a point x.
% PPCADECONSTRUCT break PPCA in pieces for saving.
% GAUSSIANNOISE3DPLOT Draws a 3D or contour plot for the GAUSSIAN noise model.
% SDLFMKERNGRADIENT Gradient of SDLFM kernel's parameters.
% LINEAREXTRACTPARAM Extract weights from a linear model.
% LFMAKERNDISPLAY Display parameters of the LFMA kernel.
% LVMSCATTERPLOTNEIGHBOURS 2-D scatter plot of the latent points with neighbourhood.
% ARDKERNDIAGGRADX Gradient of ARD kernel's diagonal with respect to X.
% SDLFMKERNPARAMINIT SDLFM kernel initialization
% WHITEHKERNDIAGGRADIENT Compute the gradient of the WHITEH kernel's diagonal wrt parameters.
% SIMKERNGRADIENT Gradient of SIM kernel's parameters.
% DEMUSPSIVM2 Try the IVM on the USPS digits data with MLP kernel.
% IVMOUT Evaluate the output of an IVM model.
% ACCLAIMREADSKEL Reads an ASF file into a skeleton structure.
% DEMGPTWOSAMPLELIF Run GP two sample code on LIF.
% MLPKERNEXTRACTPARAM Extract parameters from the MLP kernel structure.
% UNIFORMPRIORPARAMINIT Uniform prior model's parameter initialisation.
% RBFOPTIONS Default options for RBF network.
% LINEARPARAMINIT Initialise the parameters of an LINEAR model.
% CMPNDKERNDIAGGRADIENT Compute the gradient of the CMPND kernel's diagonal wrt parameters.
% DEMSTICKFGPLVM5 Model the stick man using an RBF kernel and regressive dynamics.
% ROBONEDYNAMICSDISPLAY Display the robot dynamics model. 
% RBFOUT Output of an RBF model.
% FGPLVMVITERBISEQUENCE Viterbi align a latent sequence.
% GAUSSIANWHITEKERNDIAGGRADX Gradient of gaussian white kernel's diagonal with respect to X.
% SDLFMAXSDLFMAKERNGRADIENTBLOCKIGJ 
% LFMXLFMKERNCOMPUTE Compute a cross kernel between two LFM kernels.
% POLYKERNGRADX Gradient of POLY kernel with respect to input locations.
% ROBTHREEDYNAMICSLOGLIKEGRADIENTS Gradients of the robot three dynamics wrt parameters.
% CMPNDKERNDIAGCOMPUTE Compute diagonal of CMPND kernel.
% GGXGGKERNCOMPUTE Compute a cross kernel between two GG kernels.
% DEXPKERNEXTRACTPARAM Extract parameters from the double exponential's
% GPREVERSIBLEDYNAMICSSETLATENTVALUES Set the latent values inside the model.
% SIMXRBFKERNCOMPUTE Compute a cross kernel between the SIM and RBF kernels.
% IVMDISPLAY Display parameters of an IVM model.
% LFMGRADIENTSIGMAH4VV Gradient of the function h_i(z) with respect \sigma.
% MOGSAMPLE Sample from a mixture of Gaussians model.
% KBROPTIMISE Optimise a KBR model.
% GPWRITETOFILE Write a file to be read by the C++ implementation.
% NDSIMKERNEXTRACTPARAM Extract parameters from the SIM kernel structure.
% FINDNEIGHBOURS find the k nearest neighbours for each point in Y.
% RBFPERIODICCREATE Create a RBFPERIODIC model.
% IVMREADFROMFID Load from a FID produced by the C++ implementation.
% LFMWHITEKERNDISPLAY Display parameters of the LFM-WHITE kernel.
% SDLFMVXSDLFMKERNGRADIENT Gradients of cross kernel between 2 SDLFM kernels.
% NOISELOGLIKELIHOOD Return the log-likelihood under the noise model.
% VELOTRANSKERNCOMPUTE Compute the VELOTRANS kernel given the parameters and X.
% SDLFMVXSDLFMKERNGRADIENTICBLOCK Partial derivatives initial conditions
% NDDISIMXNDSIMKERNGRADIENT Compute a cross gradient between NDDISIM and NDSIM
% SPRINGDAMPERSVISUALISE Helper code for showing an spring dampers during 2-D visualisation.
% ROBTWODYNAMICSEXTRACTPARAM Extract parameters from the robot two dynamics model.
% DEMORDEREDONEIVM1 Run a demonstration of the ordered categories noise model (linear data).
% GGXGGKERNGRADIENT Compute a cross gradient between two GG kernels.
% GPPOSTERIORMEANCOVARTEST Test the gradients of the mean and covariance.
% LFMKERNCOMPUTE Compute the LFM kernel given the parameters and X.
% DEMGPCOV2D Simple demonstration of sampling from a covariance function.
% SIGMOID The sigmoid function
% MGAUSSIANNOISELOGLIKELIHOOD Log likelihood of the data under the MGAUSSIAN noise model.
% PRIORREADPARAMSFROMFID Read prior params from C++ written FID.
% XYZANKURVISUALISE Draw the Agarwal & Triggs figure return the graphics handle.
% GGKERNEXPANDPARAM Create kernel structure from GG kernel's parameters.
% DEMGPTWOSAMPLEEB Run GP two sample code on EB.
% GPSAMPLE Create a plot of samples from a GP.
% GPTIMEDYNAMICSSEQUENCELOGLIKEGRADIENT Log-likelihood gradient for of a sequence of the GP-LVM time dynamics.
% GIBBSPERIODICKERNEXPANDPARAM Create kernel structure from GIBBSPERIODIC kernel's parameters.
% SKELMODIFY Update visualisation of skeleton data.
% LFMUPDATEKERNELS Updates the kernel representations in the LFM structure.
% INVGAMMAPRIORPARAMINIT Inverse gamma prior model's parameter initialisation.
% NGAUSSNOISESITES Site updates for noiseless Gaussian noise model.
% GAUSSIANNOISEEXTRACTPARAM Extract parameters from the GAUSSIAN noise structure.
% GPBLOCKINDICES Return indices of given block.
% MAPPINGOPTIMISE Optimise the given model.
% LINKERNCOMPUTE Compute the LIN kernel given the parameters and X.
% SDRBFKERNEXTRACTPARAM Extract parameters from the SDRBF kernel structure.
% MLPARDKERNDIAGCOMPUTE Compute diagonal of MLPARD kernel.
% DNETLOADRESULT Load a previously saved result.
% DNETOBJECTIVE Wrapper function for Density Network objective.
% BIASKERNEXPANDPARAM Create kernel structure from BIAS kernel's parameters.
% PARSENOBLEKERNELFILE Parse a kernel file from Bill Stafford Noble.
% LFMGRADIENTH42VP Gradient of the function h_i(z) with respect to some of the
% DEMGPCOVFUNCSAMPLE Sample from some different covariance functions.
% WHITEXNONEKERNGRADIENT Compute a cross gradient between WHITE and DUMMY kernels.
% DEMCMU35GPLVMFGPLVM3 Learn a GPLVM on CMU 35 data set.
% RBFPERIODICPARAMINIT RBFPERIODIC model parameter initialisation.
% XYZANKURDRAW Helper function for drawing the point cloud from Agarwal and Triggs data.
% WIENERKERNEXTRACTPARAM Extract parameters from the WIENER kernel structure.
% ISOCTAVE Returns true if the software running is Octave.
% LFMCOMPUTEH3AV Helper function for computing part of the LFMAV kernel.
% SDLFMKERNDISPLAY Display parameters of the SDLFM kernel.
% LFMJPCOMPUTEUPSILONVECTOR Upsilon vector jolt. pos. with t1, t2 limits
% DIAGKERNGRADX Gradient of DIAG kernel with respect to a point x.
% INVCMPNDKERNPARAMINIT INV.PRECISION-CMPND kernel parameter initialisation.
% ORDEREDNOISELOGLIKELIHOOD Log likelihood of the data under the ORDERED noise model.
% RBFHKERNCOMPUTE Compute the RBFH kernel given the parameters and X.
% GIBBSKERNCOMPUTE Compute the GIBBS kernel given the parameters and X.
% DNETOUT Output of an DNET model.
% RBFARDJITKERNDIAGGRADX Gradient of RBFARDJIT kernel's diagonal with respect to X.
% DEMUNLABELLEDONEIVM2 Test IVM code on a toy crescent data.
% VELOTRANSKERNEXTRACTPARAM Extract parameters from the VELOTRANS kernel structure.
% DNETGRADIENT Density Network gradient wrapper.
% SPECTRUMVISUALISE Helper code for showing an spectrum during 2-D visualisation.
% DEMOILLLE2 Demonstrate LLE on the oil data.
% DNETLOGLIKEGRADIENTS Density network gradients.
% FGPLVMCLASSVISUALISE Callback function for visualising data in 2-D.
% SDLFMVKERNEXTRACTPARAM Extract parameters from the SDLFMV kernel structure.
% GIBBSKERNSETLENGTHSCALEFUNC Set the length scale function of the GIBBS kernel.
% MOGTWODPLOT Helper function for plotting the labels in 2-D.
% SDLFMAXSDLFMVKERNCOMPUTE Cross kernel between a SDLFMA and a SDLFMV kernels.
% INDEXKERNDISPLAY Display parameters of the INDEX kernel.
% LFMCOMPUTEH3JP Helper function for computing part of the LFMJP kernel.
% PPCARECONSTRUCT Reconstruct an PPCA form component parts.
% ORDEREDNOISELIKELIHOOD Likelihood of the data under the ORDERED noise model.
% WHITEFIXEDKERNDIAGGRADIENT Compute the gradient of the WHITEFIXED kernel's diagonal wrt parameters.
% EXPKERNDISPLAY Display parameters of the EXP kernel.
% SDLFMXSDLFMKERNGRADIENTBLOCKILJ 
% STRINGSIGFIGS Convert number to a string with a number of significant digits.
% DEMSPGP1DGP4 Do a simple 1-D regression after Snelson & Ghahramani's example.
% TENSORKERNDISPLAY Display parameters of the TENSOR kernel.
% ACCLAIMGRADIENT computes the gradient of x,y,z locations wrt angles.
% LFMAACOMPUTEUPSILONDIAGVECTOR Diag. of Upsilon matrix acce. accel. 
% RBFARD2KERNPARAMINIT RBFARD2 kernel parameter initialisation.
% XYZPOPPE2JOINT
% POLYKERNGRADIENT Gradient of POLY kernel's parameters.
% IVMINIT Initialise the IVM model.
% ROBONEDYNAMICSLOGLIKELIHOOD Give the log likelihood of the robot one dynamics part.
% IVMOPTIMISE Optimise the IVM.
% GPPOSTERIORMEANCOVAR Mean and covariances of the posterior at points given by X.
% GGKERNDISPLAY Display parameters of the GG kernel.
% LFMRESULTSDYNAMICWALKING Load a results file and visualise them.
% GAUSSSAMP Sample from a Gaussian with a given covariance.
% DEMORDEREDTWOIVM1 Run a demonstration of the ordered categories noise model (circular data).
% KERNPARAMINIT Kernel parameter initialisation.
% LINARD2KERNDIAGGRADX Gradient of LINARD2 kernel's diagonal with respect to X.
% GAUSSIANKERNCOMPUTE Compute the Gaussian kernel given the parameters and X.
% IVMAPPROXLOGLIKEKERNGRAD Gradient of the approximate likelihood wrt kernel parameters.
% ROBONEDYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.
% KPCAEMBED Embed data set with kernel PCA.
% NDSIMKERNCOMPUTE Compute the NDSIM kernel with no decay given the parameters and X.
% RBFWHITEKERNGRADIENT Gradient of RBF-WHITE kernel's parameters.
% LFMOPTIONS Creates a set of default options for a LFM model.
% LFMWHITECOMPUTEPSI Helper function for comptuing part of the LFM-WHITE
% SIMSAMPLE Sample from SIM kernel
% LINARDKERNEXTRACTPARAM Extract parameters from the LINARD kernel structure.
% WHITEKERNGRADIENT Gradient of white-noise (WHITE) kernel's
% LFMAAGRADIENTUPSILONMATRIX Gradient upsilon matrix accel. accel.
% GAUSSIANWHITEKERNDIAGCOMPUTE Compute diagonal of gaussian white kernel.
% UNIFORMPRIORSETBOUNDS Set uniform prior bounds.
% SQEXPKERNCOMPUTE Compute the SQEXP kernel given the parameters and X.
% TENSORKERNEXPANDPARAM Create kernel structure from TENSOR kernel's parameters.
% SIMWHITEKERNGRADIENT Gradient of SIM-WHITE kernel's parameters.
% INVCMPNDKERNDIAGGRADX Gradient of INVCMPND kernel's diagonal with respect to X.
% IVMCREATE Create a IVM model with the IVM sparse approximaiton.
% KERNREADFROMFID Load from an FID written by the C++ implementation.
% BIASKERNDIAGGRADX Gradient of BIAS kernel's diagonal with respect to X.
% GPSEQUENCELOGLIKEGRADIENT Log-likelihood gradient for of a sequence of the GP-LVM.
% FGPLVMKERNDYNAMICSSAMPLE Sample a field from a given kernel.
% LINARD2KERNGRADIENT Gradient of LINARD2 kernel's parameters.
% WANGPRIOREXPANDPARAM Expand wang prior structure from params.
% POLYARDKERNDIAGGRADX Gradient of POLYARD kernel's diagonal with respect to X.
% WRITEDOUBLETOFID Writes a double to an FID.
% KERNCOMPUTE Compute the kernel given the parameters and X.
% GPREVERSIBLEDYNAMICSOPTIONS Return default options for GP reversible dynamics model.
% IVMSELECTVISUALISE Visualise the selected point.
% ROCHOLEXTRACT Extract the lower triangular matrix from the Cholesky structure.
% GGKERNDIAGCOMPUTE Compute diagonal of GG kernel.
% RBFKERNEXPANDPARAM Create kernel structure from RBF kernel's parameters.
% MODELPOSTERIORVAR variances of the posterior at points given by X.
% PPCAPOSTERIORVAR Mean and variances of the posterior at points given by X.
% CMPNDKERNGRADIENT Gradient of CMPND kernel's parameters.
% MLPKERNGRADIENT Gradient of MLP kernel's parameters.
% ISOMAPCREATE isomap embedding model.
% GAUSSIANWHITEKERNPARAMINIT Gaussian white kernel parameter initialisation.
% FGPLVMPOINTOBJECTIVE Wrapper function for objective of a single point in latent space and the output location..
% FGPLVMCMU35ANIMATE Animate the test data jointly with predictions.
% TENSORKERNPARAMINIT TENSOR kernel parameter initialisation.
% LFMCOMPUTEH3VP Helper function for computing part of the LFMVXLFM kernel.
% LVMTHREEDPLOT Helper function for plotting the labels in 3-D.
% IVM Initialise an IVM model.
% INVCMPNDKERNDISPLAY Display parameters of the INVCMPND kernel.
% GPOUT Evaluate the output of an Gaussian process model.
% LECREATE Laplacian eigenmap model.
% MOGUPDATEMEAN Update the means of an MOG model.
% RBFARDJITKERNDIAGGRADIENT Compute the gradient of the RBFARD2 kernel's diagonal wrt parameters.
% IVMDECONSTRUCT break IVM in pieces for saving.
% MVUOPTIMISE Optimise an MVU model.
% GPTIMEDYNAMICSOUT Evaluate the output of GPTIMEDYNAMICS.
% SIMKERNEXPANDPARAM Create kernel structure from SIM kernel's parameters.
% NONEKERNPARAMINIT NONE kernel parameter initialisation.  
% SDLFMXSDLFMKERNCOMPUTEBLOCK Computes SDLFM kernel matrix for block i,j
% GPCOMPUTEOBSERVATIONLOGLIKELIHOOD  
% KBRPARAMINIT KBR model parameter initialisation.
% SQEXPKERNDIAGGRADX Gradient of SQEXP kernel's diagonal with respect to X.
% LINKERNDIAGGRADX Gradient of LIN kernel's diagonal with respect to X.
% LINEAREXPANDPARAM Update linear model with vector of parameters.
% LFMVKERNDIAGCOMPUTE Compute diagonal of LFMV kernel.
% CENTERINGMATRIX returns the centering matrix for the given dimensionality.
% FGPLVMNEARESTNEIGHBOUR Give the number of errors in latent space for 1 nearest neighbour.
% DEMSTICKFGPLVM4 Model the stick man using an RBF kernel and 3-D latent space.
% INDEXKERNGRADIENT Gradient of INDEX kernel's parameters.
% LLECREATE Locally linear embedding model.
% OUKERNDIAGGRADX Gradient of OU kernel's diagonal with respect to t (see
% BVHPLAYDATA Play bvh motion capture data.
% DEMOPTIMISEGP Shows that there is an optimum for the covariance function length scale.
% LINKERNDISPLAY Display parameters of the LIN kernel.
% TRANSLATEKERNEXPANDPARAM Create kernel structure from TRANSLATE kernel's parameters.
% WHITEBLOCKKERNDIAGGRADX Gradient of WHITEBLOCK kernel's diagonal wrt X.
% LMCKERNDIAGCOMPUTE Compute the diagonal of the LMC kernel.
% INVSIGMOID The inverse of the sigmoid function.
% LFMVISUALISEWALKING Visualise the outputs in a latent force model
% LINARDKERNDIAGGRADX Gradient of LINARD kernel's diagonal with respect to X.
% LFMWHITEXWHITEKERNGRADIENT Compute gradient between the LFM-WHITE and
% RBFINFWHITEXWHITEKERNCOMPUTE Compute a cross kernel between the RBF-WHITE
% RBFWHITEKERNDIAGGRADX Gradient of RBF-WHITE kernel's diagonal w.r.t. t.
% INDEXARDKERNDIAGCOMPUTE Compute diagonal of INDEXARD kernel.
% WHITEHKERNEXTRACTPARAM Extract parameters from the WHITEH kernel structure.
% IVMOPTIMISENOISE Optimise the noise parameters.
% EXPKERNDIAGGRADX Gradient of EXP kernel's diagonal with respect to X.
% LINEAROUTPUTGRADX Evaluate derivatives of linear model outputs with respect to inputs.
% SDLFMVXSDLFMVKERNGRADIENT Gradients of cross kernel between 2 SDLFM kernels.
% DEMBRENDANFGPLVM2 Use the GP-LVM to model the Frey face data with FITC and back constraints.
% GPLOADRESULT Load a previously saved result.
% SDLFMXSDRBFKERNGRADIENT Cross gradient between a SDLFM and a SDRBF kernels.
% LFMAKERNEXTRACTPARAM Extract parameters from the LFMA kernel structure.
% DEMVOWELSFGPLVM1 Model the vowels data with a 2-D FGPLVM using RBF kernel and back constraints.
% MLPKERNDISPLAY Display parameters of the MLP kernel.
% FGPLVMEXPANDPARAM Expand a parameter vector into a GP-LVM model.
% LFMAXLFMKERNCOMPUTE Acceleration and position LFM kernel  
% DEMSTICKFGPLVM2 Model the stick man using an RBF kernel and dynamics.
% RBFARDKERNDIAGCOMPUTE Compute diagonal of RBFARD kernel.
% RBFWHITEKERNCOMPUTE Compute the RBF-WHITE kernel given the parameters, t1
% SDLFMVKERNEXPANDPARAM Pass parameters from params to SDLFMV kernel
% NCNMNOISELIKELIHOOD Likelihood of data under null category noise model.
% LFMCOMPUTEH4AA Helper function for computing part of the LFMAA kernel.
% RBFKERNGRADX Gradient of RBF kernel with respect to input locations.
% NUMSF2STR Convert number to a string with a number of significant digits.
% TENSORKERNSLASH Tensor kernel created by removing ith component.
% WHITEBLOCKKERNPARAMINIT WHITE BLOCK kernel parameter initialisation.
% WHITEHKERNPARAMINIT WHITEH kernel parameter initialisation.
% EXPTRANSFORM Constrains a parameter to be positive through exponentiation.
% MULTIKERNDIAGCOMPUTE Compute diagonal of MULTI kernel.
% WHITEXNONEKERNCOMPUTE Compute a cross kernel between WHITE and NONE kernels.
% MLPDISPLAY Display the multi-layer perceptron model.
% DEMWALKSITJOGDYNAMICSLEARN Learn the stick man dynamics.
% TREEFINDCHILDREN Given a tree that lists only parents, add children.
% RBFEXTRACTPARAM Wrapper for NETLAB's rbfpak.
% DEXPKERNCOMPUTE Compute the double exponential kernel,
% GAUSSIANKERNDISPLAY Display parameters of the GAUSSIAN kernel.
% WIENERKERNDISPLAY Display parameters of the WIENER kernel.
% LFMVKERNPARAMINIT LFMV kernel parameter initialisation. 
% LFMRESULTSDYNAMIC Load a results file and visualise them.
% CMPNDNOISEPOINTPLOT Plot the data-points for the CMPND noise model.
% SQEXPKERNDISPLAY Display parameters of the SQEXP kernel.
% MLPARDKERNDISPLAY Display parameters of the MLPARD kernel.
% SDLFMVXSDRBFKERNGRADIENTBLOCK Gradients of the parameters in block i,j
% LOGISTICNORMALPRIOREXTRACTPARAM Extract params from logistic-normal prior structure.
% DISIMXRBFKERNCOMPUTE Compute a cross kernel between the DISIM and RBF kernels.
% DISIMKERNDIAGGRADIENT Compute the gradient of the DISIM kernel's diagonal wrt parameters.
% MULTIKERNDIAGGRADX Gradient of MULTI kernel's diagonal with respect to X.
% LMCKERNCOMPUTE Compute the LMC kernel given the parameters and X.
% LINEARLOGLIKELIHOOD Linear model log likelihood.
% LFMJXRBFKERNCOMPUTE Compute cross kernel between the LFMJ and RBF kernels.
% FGPLVMLOGLIKEGRADIENTS Compute the gradients for the FGPLVM.
% RBFPERIODICKERNGRADIENT Gradient of RBFPERIODIC kernel's parameters.
% NCNMNOISE3DPLOT Draw a 3D or contour plot for the NCNM noise model.
% LINARD2KERNEXPANDPARAM Create kernel structure from LINARD2 kernel's parameters.
% RATQUADKERNDISPLAY Display parameters of the RATQUAD kernel.
% FGPLVMOBJECTIVEGRADIENT Wrapper function for FGPLVM objective and gradient.
% MODELWRITERESULT Write a model to file.
% LFMAVGRADIENTUPSILONMATRIX Gradient upsilon matrix accel. vel.
% NGAUSSNOISEDISPLAY Display parameters of the NGAUSS noise.
% DEMSPGP1DPLOT Plot results from 1-D sparse GP.
% POLYKERNDISPLAY Display parameters of the POLY kernel.
% VITERBIALIGN Compute the Viterbi alignment.
% RBFARDKERNPARAMINIT RBFARD kernel parameter initialisation.
% IVMMESHVALS Give the output of the IVM for contour plot display.
% NOISELIKELIHOOD Return the likelihood for each point under the noise model.
% SQEXPKERNGRADX Gradient of SQEXP kernel with respect to a point x.
% MOGMEANCOV Project a mixture of Gaussians to a low dimensional space.
% LINARD2KERNDIAGCOMPUTE Compute diagonal of LINARD2 kernel.
% RBFINFWHITEKERNGRADIENT Gradient of the parameters of the RBF-WHITE kernel
% RBFKERNPARAMINIT RBF kernel parameter initialisation.
% RBFARD2KERNEXPANDPARAM Create kernel structure from RBFARD2 kernel's parameters.
% XYZHUMANEVAREMOVEPART
% VELOTRANSKERNEXPANDPARAM Create kernel structure from VELOTRANS kernel's parameters.
% LFMXRBFKERNCOMPUTE Compute a cross kernel between the LFM and RBF kernels.
% IVMNEGLOGLIKELIHOOD Wrapper function for calling IVM likelihood.
% KBROUTPUTGRAD Evaluate derivatives of KBR model outputs with respect to parameters.
% GGWHITEXGAUSSIANWHITEKERNGRADIENT Compute gradient between the GG white
% DEMROBOTTRACES1 Wireless Robot data from University of Washington, with tailored dynamics.
% NCNMTWODPLOT Make a 2-D plot of the null category noise model.
% RBFEXPANDPARAM Update rbf model with new vector of parameters.
% GAMMAPRIORLOGPROB Log probability of Gamma prior.
% NOISE3DPLOT Draw a 3D or contour plot for the relevant noise model.
% RATQUADKERNCOMPUTE Compute the RATQUAD kernel given the parameters and X.
% GPSUBSPACEOUT
% DNETEXTRACTPARAM Extract weights and biases from an DNET.
% IVMEPUPDATEPOINT Do an EP update of a point.
% ROBTHREEDYNAMICSDISPLAY Display the robot dynamics model. 
% ROCCURVE Draw ROC curve and return labels.
% LINKERNPARAMINIT LIN kernel parameter initialisation.
% SIMWHITEKERNEXPANDPARAM Create kernel structure from SIM-WHITE kernel's
% KERNELCENTER Attempts to Center Kernel Matrix
% LFMXRBFVKERNCOMPUTE Compute a cross kernel between the LFM and RBFV kernels.
% MOGPROJECT Project a mixture of Gaussians to a low dimensional space.
% FGPLVMTEST Test the gradients of the gpCovGrads function and the fgplvm models.
% DEMROBOTWIRELESSFGPLVM1 Wireless Robot data from University of Washington, without dynamics and without back constraints.
% GIBBSKERNEXPANDPARAM Create kernel structure from GIBBS kernel's parameters.
% ACCLAIMPLAYFILE Play motion capture data from a asf and amc file.
% ROBTWODYNAMICSLOGLIKEGRADIENTS Gradients of the robot two dynamics wrt parameters.
% POLYARDKERNDIAGCOMPUTE Compute diagonal of POLYARD kernel.
% GETLINE Get a line from a file.
% FGPLVMSEQUENCELOGLIKELIHOOD Log-likelihood of a sequence for the GP-LVM.
% CMPNDKERNSETINDEX Set the indices in the compound kernel.
% LOGDET The log of the determinant when argument is positive definite.
% MLPKERNDIAGCOMPUTE Compute diagonal of MLP kernel.
% ROCHOLTRANSMULTIPLY Multiply by the transposed version of the rank one Cholesky.
% RBFARD2KERNGRADX Gradient of RBFARD2 kernel with respect to input locations.
% KERNDIAGCOMPUTE Compute the kernel given the parameters and X.
% NDSIMKERNEXPANDPARAM Create kernel structure from NDSIM kernel's parameters.
% NDDISIMKERNEXPANDPARAMTRANSFORMSETTINGS Create kernel structure from DISIM kernel's parameters.
% DYNAMICSTEST Run some tests on the specified dynamics model.
% LFMGRADIENTH31 Gradient of the function h_i(z) with respect to some of the
% LFMKERNDISPLAY Display parameters of the LFM kernel.
% LINARD2KERNCOMPUTE Compute the LINARD2 kernel given the parameters and X.
% ISOMAPRECONSTRUCT Reconstruct an isomap form component parts.
% SIMKERNEXTRACTPARAM Extract parameters from the SIM kernel structure.
% DEMTWOCLUSTERS1
% IVMKERNELGRADIENT Gradient of likelihood approximation wrt kernel parameters.
% SRBFHKERNGRADIENT Gradient of the parameters of a SRBFH kernel.
% OUKERNEXPANDPARAM Create kernel structure from OU kernel's parameters
% EXPKERNGRADX Gradient of EXP kernel with respect to a point x.
% LFMSAMPLE Sample from LFM kernel
% MODELGETOUTPUTWEIGHTS Wrapper function to return output weight and bias matrices.
% DIAGKERNPARAMINIT DIAG kernel parameter initialisation.
% LFMGRADIENTSIGMAH Gradient of the function h_i(z) with respect \sigma.
% GIBBSPERIODICKERNCOMPUTE Compute the GIBBSPERIODIC kernel given the parameters and X.
% PSKERNELOBJECTIVE Likelihood approximation for point set IVM.
% GAUSSIANPRIORGRADIENT Gradient wrt x of the log Gaussian prior.
% GIBBSPERIODICKERNDIAGGRADIENT Compute the gradient of the GIBBSPERIODIC kernel's diagonal wrt parameters.
% LFMCREATE Create a LFM model.
% GGWHITEXGGWHITEKERNCOMPUTE Compute a cross kernel between two GG white kernels.
% NCNMCONTOUR Special contour plot showing null category region.
% KERNEXTRACTPARAMTRANSFORMSETTINGS Extract parameter transform settings from kernel structure.
% GGKERNDIAGGRADIENT Compute gradient of the diagonal of GG kernel.
% GIBBSKERNGRADIENT Gradient of GIBBS kernel's parameters.
% WHITEFIXEDKERNGRADIENT Gradient of WHITEFIXED kernel's parameters.
% DEMSPGP1DGP3 Do a simple 1-D regression after Snelson & Ghahramani's example.
% LFMVXLFMVKERNGRADIENT Compute a cross gradient between a LFMV and a LFMV.
% DISIMXDISIMKERNGRADIENT Compute a cross gradient between two DISIM kernels.
% FINDACYCLICNEIGHBOURS find the k nearest neighbours for each point in Y preventing cycles in the graph.
% LFMWHITEKERNCOMPUTE Compute the LFM-WHITE kernel given the parameters, t1
% FILEKERNEXTRACTPARAM Extract parameters from the FILE kernel structure.
% SDLFMXSDRBFKERNCOMPUTEBLOCK Cross kernel between SDLFM and SDRBF for i,j
% ZEROAXES A function to move the axes crossing point to the origin.
% WRITEBOOLTOFID Writes a boolean to an FID.
% GGXGAUSSIANKERNGRADX Compute gradient between the GG and GAUSSIAN
% LFMAKERNDIAGCOMPUTE Compute diagonal of LFMAXLFMA kernel.
% GPDYNAMICSSEQUENCELOGLIKELIHOOD Return the log likelihood of a given latent sequence.
% XYZANKURVISUALISE2
% MLPLOGLIKEHESSIAN Multi-layer perceptron Hessian.
% GETSYMBOLS Get a cell array of different plot symbols.
% GGKERNEXTRACTPARAM Extract parameters from the GG kernel structure.
% RBFKERNEXTRACTPARAM Extract parameters from the RBF kernel structure.
% MULTIKERNCOMPUTEBLOCK
% SDLFMKERNMEANCOVPARTIAL Helper function for derivatives in SDLFM kernel
% DEMSWISSROLLLLE3 Demonstrate LLE on the oil data.
% DIST1	Calculates absolute distance (i.e. L1 norm) between two sets of
% DEXPKERNDISPLAY Display parameters of the double exponential kernel.
% RBFPERIODIC2KERNDISPLAY Display parameters of the RBFPERIODIC2 kernel.
% RBFINFWHITEKERNCOMPUTE Compute the RBF-WHITE kernel (with integration limits
% MULTIKERNPARAMINIT MULTI kernel parameter initialisation.
% XYZANKURERROR Computes the error between two poses in xyz format 
% XYZHUMANEVAGENERATEMOVIE
% MULTIMODELLOGLIKEGRADIENTS Gradient of MULTIMODEL model log likelihood with respect to parameters.
% DEMCMU35SEQUENCEOPTIMISE 
% CMPNDNOISELIKELIHOOD Likelihood of the data under the CMPND noise model.
% GPOPTIONS Return default options for GP model.
% FGPLVMOPTIMISESEQUENCE Optimise the postion of a latent sequence.
% BVHCONNECTIONMATRIX Compute the connection matrix for the structure.
% MODELSETOUTPUTWEIGHTS Wrapper function to return set output weight and bias matrices.
% MULTIKERNDISPLAY Display parameters of the MULTI kernel.
% GPCREATE Create a GP model with inducing varibles/pseudo-inputs.
% LLEEMBED Embed data set with LLE.
% LFMKERNEXPANDPARAM Create kernel structure from LFM kernel's parameters.
% GPTIMEDYNAMICSLOGLIKEGRADIENTS Gradients of the GP dynamics wrt parameters.
% MLPARDKERNEXTRACTPARAM Extract parameters from the MLPARD kernel structure.
% GAUSSIANKERNEXTRACTPARAM Extract parameters from the gaussian kernel structure.
% GAUSSIANNOISENUG Compute nu and g for GAUSSIAN noise model.
% LFMXLFMKERNGRADIENT Compute a cross gradient between two LFM kernels.
% XLOGY z = x*log(y) returns zero if x=y=0
% RATQUADKERNDIAGGRADX Gradient of RATQUAD kernel's diagonal with respect to X.
% ARDKERNDIAGCOMPUTE Compute diagonal of ARD kernel.
% NDSIMKERNEXTRACTPARAMTRANSFORMSETTINGS Extract parameter transform settings from the SIM kernel structure.
% ROBTWODYNAMICSLOGLIKELIHOOD Give the log likelihood of the robot one dynamics part.
% LFMGRADIENTH41 Gradient of the function h_i(z) with respect to some of the
% DNETESTEP Do an E-step (update importance weights) on an Density Network model.
% GGKERNCOMPUTE Compute the GG kernel given the parameters and X.
% DNETOPTIONS Options for a density network.
% OUKERNCOMPUTE Compute the Ornstein-Uhlenbeck (OU) kernel arising from the
% GIBBSKERNPARAMINIT GIBBS kernel parameter initialisation.
% SIMWHITEKERNPARAMINIT SIM-WHITE kernel parameter initialisation.
% MVURECONSTRUCT Reconstruct an MVU form component parts.
% SDLFMAXSDLFMKERNCOMPUTE Cross kernel between a SDLFMA and a SDLFM kernels.
% RBFWHITEKERNGRADX Gradient of RBF-WHITE kernel with respect to a point t.
% LINARD2KERNEXTRACTPARAM Extract parameters from the LINARD2 kernel structure.
% DEMROBOTWIRELESSFGPLVM2 Wireless Robot data from University of Washington, without dynamics and without back constraints.
% RBFHKERNEXTRACTPARAM Extract parameters from the RBFH kernel structure.
% LOGISTICNORMALPRIORGRADIENT Gradient wrt x of the logistic-normal prior.
% FGPLVMREADFROMFID Load from a FID produced by the C++ implementation.
% OUKERNGRADX Gradient of OU kernel with respect to a point x (see
% NCNMNOISENUG Update nu and g parameters associated with null category noise model.
% DEMSWISSROLLLLE4 Demonstrate LLE on the oil data.
% MLPEXTRACTPARAM Extract weights and biases from an MLP.
% BVHMODIFY Helper code for visualisation of bvh data.
% OPTIMIMINIMIZE Wrapper for Carl Rasmussen's minimize function.
% LFMCOMPUTEUPSILONDIAGVECTOR 
% CMPNDNOISEEXPANDPARAM Create noise structure from CMPND noise's parameters.
% LEOPTIMISE Optimise an LE model.
% KLDIVGAUSSIAN Give the KL divergence between two Gaussians.
% KERNWRITEPARAMSTOFID Write the kernel parameters to a stream.
% VECTORMODIFY Helper code for visualisation of vectorial data.
% ORDEREDNOISEOUT Compute the output of the ORDERED noise given the input mean and variance.
% MGAUSSIANNOISEPARAMINIT MGAUSSIAN noise parameter initialisation.
% PROBITNOISEEXPANDPARAM Create noise structure from PROBIT noise's parameters.
% LINARDKERNPARAMINIT LINARD kernel parameter initialisation.
% SDLFMAXSDLFMKERNCOMPUTEBLOCK Computes SDLFM kernel matrix for block i,j
% ACCLAIMLOADCHANNELS Load the channels from an AMC file.
% NOISECREATE Initialise a noise structure.
% LFMGRADIENTH42AP Gradient of the function h_i(z) with respect to some of the
% BIASKERNCOMPUTE Compute the BIAS kernel given the parameters and X.
% LFMCOMPUTEH4VP Helper function for computing part of the LFMVXLFM kernel.
% DEMOILFGPLVM8 Oil data with variational sparse approximation.
% DNETLOGLIKELIHOOD Density network log likelihood.
% NDDISIMKERNEXTRACTPARAMTRANSFORMSETTINGS Extract parameter transform settings from the NDDISIM kernel structure.
% DEMSWISSROLLFULLLLE2 Demonstrate LLE on the oil data.
% LFMGRADIENTSIGMAH4AV Gradient of the function h_i(z) with respect \sigma.
% NOISEEXPECTATIONLOGLIKELIHOOD Return the expectation of the log likelihood.
% WHITEXRBFKERNGRADIENT Compute a cross gradient between WHITE and RBF kernels.
% GPPOINTLOGLIKELIHOOD Log-likelihood of a test point for a GP.
% MODELADDDYNAMICS Add a dynamics kernel to the model.
% SDLFMKERNGRADIENTCONSTANT Gradients for constants for the SDLFM kernel
% RBFPERIODICKERNDIAGGRADIENT Compute the gradient of the RBFPERIODIC kernel's diagonal wrt parameters.
% IVMPRINTPLOT Make a 3-D or contour plot of the IVM.
% NORMUNIPRIOREXPANDPARAM Expand Normal uniform prior structure from param vector.
% DEMUNLABELLEDIVM2 Test IVM code on a toy crescent data.
% LFMJACOMPUTEUPSILONMATRIX Upsilon matrix jolt. accel. with t1, t2 limits
% CMPNDKERNDIAGGRADX Gradient of CMPND kernel's diagonal with respect to X.
% DEMOIL100FGPLVM1 Oil100 data with fully independent training conditional.
% GIBBSKERNDIAGGRADX Gradient of GIBBS kernel's diagonal with respect to X.
% RBFPERIODICKERNDIAGGRADX Gradient of RBFPERIODIC kernel's diagonal with respect to X.
% GPCOMPUTEALPHA Update the vector `alpha' for computing posterior mean quickly.
% SMALLRANDEMBED Embed data set with small random values.
% WIENERKERNPARAMINIT WIENER kernel parameter initialisation.
% NOISEEXPANDPARAM Expand the noise model's parameters from params vector.
% LVMCLASSCLICKVISUALISE Callback function for visualising data in 2-D with clicks.
% GPREVERSIBLEDYNAMICSEXPANDPARAM Place the parameters vector into the model for GP dynamics.
% TREEFINDPARENTS Given a tree that lists only children, add parents.
% LFMGRADIENTSIGMAUPSILON Gradient of the function \upsilon(z) with respect
% KERNGRADIENT Compute the gradient wrt the kernel's parameters.
% GAUSSIANNOISEOUT Compute the output of the GAUSSIAN noise given the input mean and variance.
% POLYKERNDIAGCOMPUTE Compute diagonal of POLY kernel.
% DEXPKERNGRADIENT Gradient of the double exponential kernel's parameters.
% WIENERKERNEXPANDPARAM Create kernel structure from WIENER kernel's parameters.
% GAMMAPRIORGRADIENT Gradient wrt x of the gamma prior.
% VIVMRUNDATASET Try the virtual IVM on a data set and save the results.
% SPARSEKERNDISPLAY Display parameters of the SPARSE kernel.
% MULTIKERNGRADIENTBLOCKX
% DEMOILLLE3 Demonstrate LLE on the oil data.
% NGAUSSNOISELOGLIKELIHOOD Log likelihood of the data under the NGAUSS noise model.
% LFMKERNDIAGGRADIENT Compute the gradient of the LFM kernel's diagonal wrt parameters.
% DEMOILFGPLVM2 Oil data with fully independent training conditional, and MLP back constraints.
% DOUBLEMATRIXREADFROMFID Read a full matrix from an FID.
% GPUPDATEAD Update the representations of A and D associated with the model.
% RBFINFWHITEXWHITEKERNGRADIENT Compute gradient between the RBF-WHITE kernel
% MODELGRADIENT Gradient of error function to minimise for given model.
% INDEXKERNPARAMINIT INDEX kernel parameter initialisation.
% INVCMPNDKERNSETINDEX Set the indices in the inv. compound kernel.
% WHITEBLOCKKERNGRADX Gradient of WHITEBLOCK kernel wrt input locations.
% LFMLOGLIKELIHOOD Compute the log likelihood of a LFM model.
% KERNPRIORGRADIENT Compute gradient terms associated with kernel priors.
% LINKERNEXTRACTPARAM Extract parameters from the LIN kernel structure.
% LFMEXTRACTPARAM Extract the parameters of an LFM model.
% GPDYNAMICSSETLATENTVALUES Set the latent values inside the model.
% MGAUSSIANNOISE3DPLOT Draws a 3D or contour plot for the MGAUSSIAN noise model.
% BIASKERNGRADIENT Gradient of BIAS kernel's parameters.
% LFMGRADIENTH Gradient of the function h_i(z) with respect to some of the
% DEMCLASSIFICATIONTWOIVM1 IVM for classification on a data-set sampled from a GP
% LFMXRBFKERNGRADIENT Compute gradient between the LFM and RBF kernels.
% DEMGPSAMPLE Simple demonstration of sampling from a covariance function.
% ARDKERNPARAMINIT ARD kernel parameter initialisation.
% STACK Return column stacked vector of given matrix.
% MLPARDKERNPARAMINIT MLPARD kernel parameter initialisation.
% SHEFFIELDMLTOOLBOXES External dependencies for Sheffield ML toolboxes.
% MGAUSSIANNOISELIKELIHOOD Likelihood of the data under the MGAUSSIAN noise model.
% NOISEPOINTPLOT Plot the data-points for the given noise model.
% GPREVERSIBLEDYNAMICSLOGLIKELIHOOD Give the log likelihood of the dynamics part.
% RBFPERIODIC2KERNEXPANDPARAM Create kernel structure from RBFPERIODIC2 kernel's parameters.
% FGPLVMLOADRESULT Load a previously saved result.
% FGPLVMDYNAMICSSAMPLE Sample a field from the GP.
% ROBTWODYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.
% MOGESTEP Do an E-step on an MOG model.
% WANGPRIORLOGPROB Log probability of Wang prior.
% IVMUPDATESITES Update site parameters.
% LINEARCREATE Create a linear model.
% DEMPROBIT1 Test IVM code on a toy crescent data.
% RBFARDKERNEXTRACTPARAM Extract parameters from the RBFARD kernel structure.
% ROCHOLBACKSUB Backsubstitute the representation of the rank one Cholesky.
% PRIORPARAMINIT Prior model's parameter initialisation.
% DEMSILHOUETTEGP2 Model silhouette data with independent MLP GPs.
% RATQUADKERNEXTRACTPARAM Extract parameters from the RATQUAD kernel structure.
% LFMWHITEKERNEXPANDPARAM Create kernel structure from LFM-WHITE kernel's
% MGAUSSIANNOISEGRADIENTPARAM Gradient of MGAUSSIAN noise's parameters.
% LFMWHITEKERNDIAGGRADX Gradient of LFM-WHITE kernel's diagonal w.r.t. t.
% FGPLVMSEQUENCEOBJECTIVEGRADIENT Wrapper function for objective
% KERNEXTRACTPARAM Extract parameters from kernel structure.
% FILEKERNDISPLAY Display parameters of the FILE kernel.
% OUKERNGRADIENT Gradient of OU kernel's parameters (see ouKernCompute or
% GPREVERSIBLEDYNAMICSDISPLAY Display a GP dynamics model.
% SIMWHITEXRBFWHITEKERNGRADIENT Compute a cross gradient between a SIM-WHITE
% SHEATKERNGRADIENT Gradient of the parameters of a SHEAT kernel.
% LFMGRADIENTSIGMAH4AA Gradient of the function h_i(z) with respect \sigma.
% ESCAPETEXT Add back slashes to escape existing backslashes in a text.
% MODELPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% GPSCALEBIASGRADIENT Compute the log likelihood gradient wrt the scales.
% TENSORKERNGRADX Gradient of TENSOR kernel with respect to a point x.
% DEMCMU35ANIMATE Animate reconstructed right leg and body.
% GPCOMPUTETRANSLATIONLOGLIKELIHOOD  
% DEMOIL100FGPLVM2 Oil100 data with FGPLVM.
% GAMMAPDF PDF for the Gamma distribution.
% GGWHITEKERNEXTRACTPARAM Extract parameters from the GG WHITE kernel structure.
% MULTIMODELCREATE Create a MULTIMODEL model.
% GAMMAPRIOREXTRACTPARAM Extract params from gamma prior structure.
% DIAGKERNEXPANDPARAM Create kernel structure from DIAG kernel's parameters.
% MULTIKERNEXTRACTPARAMTRANSFORMSETTINGS Extract parameter transform settings 
% DEMSWISSROLLFULLLLE1 Demonstrate LLE on the oil data.
% DEXPKERNEXPANDPARAM Create a kernel structure from the double exponential
% FGPLVMSEQUENCEGRADIENT Wrapper function for gradient of a latent sequence.
% NOISEEXTRACTPARAM Extract the noise model's parameters.
% LFMKERNDIAGCOMPUTE Compute diagonal of LFM kernel.
% RBFPERIODICDISPLAY Display parameters of the RBFPERIODIC model.
% HEATKERNPARAMINIT HEAT kernel parameter initialisation.
% LFMVPCOMPUTEUPSILONMATRIX Upsilon matrix vel. pos. with t1, t2 limits
% DISTANCEWARP Dynamic Time Warping Algorithm
% SDLFMAXSDLFMKERNGRADIENTBLOCK Gradients of the parameters in block i,j
% ORDEREDNOISEPARAMINIT ORDERED noise parameter initialisation.
% GPREADFROMFILE Load a file produced by the C++ implementation.
% RBFARDKERNGRADIENT Gradient of RBFARD kernel's parameters.
% GPCOMPUTEM Compute the matrix m given the model.
% DISIMKERNEXTRACTPARAM Extract parameters from the DISIM kernel structure.
% LFMAPGRADIENTUPSILONVECTOR Gradient upsilon vector accel. pos.
% GIBBSPERIODICKERNEXTRACTPARAM Extract parameters from the GIBBSPERIODIC kernel structure.
% WHITEBLOCKKERNEXPANDPARAM Fill WHITEBLOCK kernel structure with params.
% DEXPKERNDIAGGRADX Gradient of the double exponential kernel's diagonal
% HEATXRBFHKERNGRADIENT Gradient wrt parameters between a HEAT and a RBFH.
% LVMLOADDATA Load a latent variable model dataset.
% SQEXPKERNDIAGGRADIENT Compute the gradient of the SQEXP kernel's diagonal wrt parameters.
% KERNWRITETOFID Load from an FID written by the C++ implementation.
% LFMWHITEKERNDIAGCOMPUTE Compute diagonal of LFM-WHITE kernel.
% BVH2XYZ Compute XYZ values given structure and channels.
% PRIORLOGPROB Log probability of given prior.
% DEMREGRESSIONTWOIVM1 The data-set is sampled from a GP with known parameters.
% GPOBJECTIVEGRADIENT Wrapper function for GP objective and gradient.
% MODELOPTIMISE Optimise the given model.
% LAPLACEPRIOREXPANDPARAM Expand Laplace prior structure from param vector.
% LVMCLASSVISUALISEPATH Latent variable model path drawing in latent space.
% CMPNDTIEPARAMETERS Tie parameters together.
% POLYKERNEXPANDPARAM Create kernel structure from POLY kernel's parameters.
% SPECTRALUPDATEX Update the latent representation for spectral model.
% FGPLVMRESULTSCPP Load a results file and visualise them.
% GIBBSPERIODICKERNDIAGGRADX Gradient of GIBBSPERIODIC kernel's diagonal with respect to X.
% DNETOPTIMISE Optimise an DNET model.
% FGPLVMCREATE Create a GPLVM model with inducing variables.
