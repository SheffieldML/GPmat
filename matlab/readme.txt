FGPLVM software
Version 0.132		Friday 24 Feb 2006 at 19:26
Copyright (c) 2006 Neil D. Lawrence

The FGPLVM toolbox is a new implementation of the GP-LVM that uses the Pseudo-Input method of Snelson and Ghahramani (NIPS 2005) for sparsification and efficiency improvements. 

Version 0.132
-------------

Learning with missing data fully implemented across all models. Two big speed improvements on the fitc approximation (thanks to Ed Snelson for pointing out how slow it was!).

Version 0.131
-------------

Added learning with missing data for the FTC and reversible dynamics model.

Version 0.13
------------

This version includes a much cleaner way of incorporating different dynamics models. It is released in line two imminent reports on learning large scale Gaussian processes and learning with back constraints.

Version 0.11
------------

This version now includes the Snelson-Ghahramani approximation (called FITC by Quinonero-Candela and Rasmussen) and the partially independent training criterion (PITC). Additionally the approximations can be used in standard Gaussian process regression.

Version 0.1
-----------

In the first release, only the projected latent variables approximation of Seeger et al is implemented. The toolbox also implements back-constraints as proposed by Lawrence and Quinonero-Candela

The first release containing a couple of examples on the oil data (demOil1.m and demOil2.m) and dynamics (demStick1.m and demStick2.m). The toolbox can load in the C++ files with dynamics associated and (through the mocapResultsCppBvh in the MOCAP toolbox) can run motion capture files with dynamics.

MATLAB Files
------------

Matlab files associated with the toolbox are:

cmdsRoadData.m: This script uses classical MDS to visualise some road distance data.
demOil1.m: Oil data with fully independent training conditional.
demOil2.m: Oil data with fully independent training conditional, and MLP back constraints.
demOil3.m: Oil data with deterministic training conditional.
demOil4.m: Oil data with deterministic training conditional, and MLP back constraints.
demOil5.m: Oil data with partially independent training conditional.
demOil6.m: Oil data with partially independent training conditional, and MLP back constraints.
demOilTest.m: Take some test data for the robot and navigate with it.
demRobotTraces1.m: Wireless Robot data from University of Washington, with tailored dynamics.
demRobotWireless1.m: Wireless Robot data from University of Washington, without dynamics and without back constraints.
demRobotWireless2.m: Wireless Robot data from University of Washington, without dynamics and without back constraints.
demRobotWireless3.m: Wireless Robot data from University of Washington with dynamics and no back constraints.
demRobotWireless4.m: Wireless Robot data from University of Washington with dynamics and back constraints.
demRobotWirelessNavigate.m: Take some test data for the robot and navigate with it.
demSpgp1d1.m: Do a simple 1-D regression after Snelson & Ghahramani's example.
demSpgp1d2.m: Do a simple 1-D regression after Snelson & Ghahramani's example.
demSpgp1d3.m: Do a simple 1-D regression after Snelson & Ghahramani's example.
demSpgp1d4.m: Do a simple 1-D regression after Snelson & Ghahramani's example.
demSpgp1dPlot.m: Plot results from 1-D sparse GP.
demStick1.m: Model the stick man using an RBF kernel.
demStick2.m: Model the stick man using an RBF kernel and dynamics.
demStick3.m: Model the stick man using an RBF kernel and mlp back constraints.
demStick4.m: Model the stick man using an RBF kernel and 3-D latent space.
demVowels1.m: Model the vowels data with a 2-D FGPLVM using RBF kernel.
demVowels2.m: Model the vowels data with a 2-D FGPLVM using RBF kernel and back constraints.
demVowels3.m: Model the vowels data with a 2-D FGPLVM using RBF kernel and back constraints, but without PCA initialisation.
demVowelsIsomap.m: Model the vowels data with a 2-D FGPLVM using RBF kernel.
demVowelsLle.m: Model the vowels data with a 2-D FGPLVM using RBF kernel.
demWalkSitJog1.m: Model the stick man using an RBF kernel.
demWalkSitJog2.m: Model the stick man using an RBF kernel.
demWalkSitJog3.m: Model the stick man using an RBF kernel.
demWalkSitJog4.m: Model the stick man using an RBF kernel.
demWalkSitJogDynamicsLearn.m: Learn the stick man dynamics.
dynamicsTest.m: Run some tests on the specified dynamics model.
fgplvmAddDynamics.m: Add a dynamics kernel to the model.
fgplvmClassVisualise.m: Callback function for visualising data in 2-D.
fgplvmCovGradsTest.m: Test the gradients of the covariance.
fgplvmCreate.m: Create a GPLVM model with inducing varibles.
fgplvmDisplay.m: Display an FGPLVM model.
fgplvmDynamicsFieldPlot.m: 2-D field plot of the dynamics.
fgplvmDynamicsPlot.m: 2-D scatter plot of the latent points.
fgplvmDynamicsPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
fgplvmDynamicsRun.m: Visualise the manifold.
fgplvmDynamicsSample.m: Sample a field from the GP.
fgplvmExpandParam.m: Expand a parameter vector into a GP-LVM model.
fgplvmExtractParam.m: Extract a parameter vector from a GP-LVM model.
fgplvmFieldPlot.m: 2-D field plot of the dynamics.
fgplvmGradient.m: GP-LVM gradient wrapper.
fgplvmKernDynamicsSample.m: Sample a field from a given kernel.
fgplvmLoadResult.m: Load a previously saved result.
fgplvmLogLikeGradients.m: Compute the gradients of the EZFT sparse covariance.
fgplvmLogLikelihood.m: Log-likelihood for a GP-LVM.
fgplvmNearestNeighbour.m: Give the number of errors in latent space for 1 nearest neighbour.
fgplvmObjective.m: Wrapper function for GPLVM objective.
fgplvmOptimise.m: Optimise the inducing variable based kernel.
fgplvmOptimisePoint.m: Optimise the postion of a point.
fgplvmOptions.m: Return default options for FGPLVM model.
fgplvmPointGradient.m: Wrapper function for gradient of a single point.
fgplvmPointLogLikeGradient.m: Log-likelihood gradient for of a point of the GP-LVM.
fgplvmPointLogLikelihood.m: Log-likelihood of a point for the GP-LVM.
fgplvmPointObjective.m: Wrapper function for objective of a single point.
fgplvmPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
fgplvmPrintPlot.m: Print latent space for learnt model.
fgplvmReadFromFID.m: Load from a FID produced by the C++ implementation.
fgplvmReadFromFile.m: Load a file produced by the c++ implementation.
fgplvmResultsDynamic.m: Load a results file and visualise them.
fgplvmTest.m: Test the gradients of the gpCovGrads function and the fgplvm models.
fgplvmVisualise.m: Visualise the manifold.
gpBlockIndices.m: Return indices of given block.
gpComputeAlpha.m: Update the vector `alpha' for computing posterior mean quickly.
gpComputeM.m: Compute the matrix m given the model.
gpCovGrads.m: Sparse objective function gradients wrt Covariance functions for inducing variables.
gpCovGradsTest.m: Test the gradients of the covariance.
gpCreate.m: Create a GP model with inducing varibles/pseudo-inputs.
gpDataIndices.m: Return indices of present data.
gpDisplay.m: Display a Gaussian process model.
gpDynamicsCreate.m: Create the dynamics model. 
gpDynamicsDisplay.m: Display a GP dynamics model.
gpDynamicsExpandParam.m: Place the parameters vector into the model for GP dynamics.
gpDynamicsExtractParam.m: Extract parameters from the GP dynamics model.
gpDynamicsLatentGradients.m: Gradients of the X vector given the dynamics model.
gpDynamicsLogLikeGradients.m: Gradients of the GP dynamics wrt parameters.
gpDynamicsLogLikelihood.m: Give the log likelihood of the dynamics part.
gpDynamicsSamp.m: Sample from the dynamics for a given input.
gpDynamicsSetLatentValues.m: Set the latent values inside the model.
gpExpandParam.m: Expand a parameter vector into a GP model.
gpExtractParam.m: Extract a parameter vector from a GP model.
gpGradient.m: Gradient wrapper for a GP model.
gpLogLikeGradients.m: Compute the gradients for the parameters and X.
gpLogLikelihood.m: Compute the log likelihood of a GP.
gpObjective.m: Wrapper function for GP objective.
gpOptimise.m: Optimise the inducing variable based kernel.
gpOptions.m: Return default options for GP model.
gpOut.m: Evaluate the output of an Gaussian process model.
gpPosteriorGradMeanVar.m: Gadient of the mean and variances of the posterior at points given by X.
gpPosteriorMeanCovar.m: Mean and covariances of the posterior at points given by X.
gpPosteriorMeanVar.m: Mean and variances of the posterior at points given by X.
gpReversibleDynamicsCreate.m: Create the dynamics model. 
gpReversibleDynamicsDisplay.m: Display a GP dynamics model.
gpReversibleDynamicsExpandParam.m: Place the parameters vector into the model for GP dynamics.
gpReversibleDynamicsExtractParam.m: Extract parameters from the GP reversible dynamics model.
gpReversibleDynamicsLatentGradients.m: Gradients of the X vector given the dynamics model.
gpReversibleDynamicsLogLikeGradients.m: Gradients of the GP reversible dynamics wrt parameters.
gpReversibleDynamicsLogLikelihood.m: Give the log likelihood of the dynamics part.
gpReversibleDynamicsOptions.m: Return default options for GP reversible dynamics model.
gpReversibleDynamicsSamp.m: Sample from the dynamics for a given input.
gpReversibleDynamicsSetLatentValues.m: Set the latent values inside the model.
gpScaleBiasGradient.m: Compute the gradient of the scale and bias.
gpUpdateKernels.m: Update the kernels that are needed.
modelLatentGradients.m: Gradients of the latent variables for dynamics models in the GPLVM.
modelSetLatentValues.m: Set the latent variables for dynamics models in the GPLVM.
robOneDynamicsCreate.m: Create the dynamics model. 
robOneDynamicsDisplay.m: Display the robot dynamics model. 
robOneDynamicsExpandParam.m: Place the parameters vector into the model for first robot dynamics.
robOneDynamicsExtractParam.m: Extract parameters from the robot one dynamics model.
robOneDynamicsLatentGradients.m: Gradients of the X vector given the dynamics model.
robOneDynamicsLogLikeGradients.m: Gradients of the robot one dynamics wrt parameters.
robOneDynamicsLogLikelihood.m: Give the log likelihood of the robot one dynamics part.
robOneDynamicsSetLatentValues.m: Set the latent values inside the model.
robThreeDynamicsCreate.m: Create the dynamics model. 
robThreeDynamicsDisplay.m: Display the robot dynamics model. 
robThreeDynamicsExpandParam.m: Place the parameters vector into the model for first robot dynamics.
robThreeDynamicsExtractParam.m: Extract parameters from the robot three dynamics model.
robThreeDynamicsLatentGradients.m: Gradients of the X vector given the dynamics model.
robThreeDynamicsLogLikeGradients.m: Gradients of the robot three dynamics wrt parameters.
robThreeDynamicsLogLikelihood.m: Give the log likelihood of the robot three dynamics part.
robThreeDynamicsSetLatentValues.m: Set the latent values inside the model.
robThreeSetLatentValues.m: Set the latent values inside the model.
robTwoDynamicsCreate.m: Create the dynamics model. 
robTwoDynamicsDisplay.m: Display the robot dynamics model. 
robTwoDynamicsExpandParam.m: Place the parameters vector into the model for first robot dynamics.
robTwoDynamicsExtractParam.m: Extract parameters from the robot two dynamics model.
robTwoDynamicsLatentGradients.m: Gradients of the X vector given the dynamics model.
robTwoDynamicsLogLikeGradients.m: Gradients of the robot two dynamics wrt parameters.
robTwoDynamicsLogLikelihood.m: Give the log likelihood of the robot one dynamics part.
robTwoDynamicsSetLatentValues.m: Set the latent values inside the model.
