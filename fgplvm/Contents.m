% FGPLVM toolbox
% Version 0.132		Friday 24 Feb 2006 at 19:26
% Copyright (c) 2006 Neil D. Lawrence
% 
% CMDSROADDATA This script uses classical MDS to visualise some road distance data.
% DEMOIL1 Oil data with fully independent training conditional.
% DEMOIL2 Oil data with fully independent training conditional, and MLP back constraints.
% DEMOIL3 Oil data with deterministic training conditional.
% DEMOIL4 Oil data with deterministic training conditional, and MLP back constraints.
% DEMOIL5 Oil data with partially independent training conditional.
% DEMOIL6 Oil data with partially independent training conditional, and MLP back constraints.
% DEMROBOTWIRELESSNAVIGATE Take some test data for the robot and navigate with it.
% DEMROBOTTRACES1 Wireless Robot data from University of Washington, with tailored dynamics.
% DEMROBOTWIRELESS1 Wireless Robot data from University of Washington, without dynamics and without back constraints.
% DEMROBOTWIRELESS2 Wireless Robot data from University of Washington, without dynamics and without back constraints.
% DEMROBOTWIRELESS4 Wireless Robot data from University of Washington with dynamics and no back constraints.
% DEMROBOTWIRELESS4 Wireless Robot data from University of Washington with dynamics and back constraints.
% DEMROBOTWIRELESSNAVIGATE Take some test data for the robot and navigate with it.
% DEMSPGP1D1 Do a simple 1-D regression after Snelson & Ghahramani's example.
% DEMSPGP1D2 Do a simple 1-D regression after Snelson & Ghahramani's example.
% DEMSPGP1D3 Do a simple 1-D regression after Snelson & Ghahramani's example.
% DEMSPGP1D4 Do a simple 1-D regression after Snelson & Ghahramani's example.
% DEMSPGP1DPLOT Plot results from 1-D sparse GP.
% DEMSTICK1 Model the stick man using an RBF kernel.
% DEMSTICK2 Model the stick man using an RBF kernel and dynamics.
% DEMSTICK3 Model the stick man using an RBF kernel and mlp back constraints.
% DEMSTICK4 Model the stick man using an RBF kernel and 3-D latent space.
% DEMVOWELS1 Model the vowels data with a 2-D FGPLVM using RBF kernel.
% DEMVOWELS2 Model the vowels data with a 2-D FGPLVM using RBF kernel and back constraints.
% DEMVOWELS3 Model the vowels data with a 2-D FGPLVM using RBF kernel and back constraints, but without PCA initialisation.
% DEMVOWELSISOMAP Model the vowels data with a 2-D FGPLVM using RBF kernel.
% DEMVOWELSLLE Model the vowels data with a 2-D FGPLVM using RBF kernel.
% DEMWALKSITJOG1 Model the stick man using an RBF kernel.
% DEMWALKSITJOG2 Model the stick man using an RBF kernel.
% DEMWALKSITJOG3 Model the stick man using an RBF kernel.
% DEMWALKSITJOG4 Model the stick man using an RBF kernel.
% DEMWALKSITJOGDYNAMICSLEARN Learn the stick man dynamics.
% DYNAMICSTEST Run some tests on the specified dynamics model.
% FGPLVMADDDYNAMICS Add a dynamics kernel to the model.
% FGPLVMCLASSVISUALISE Callback function for visualising data in 2-D.
% FGPLVMCOVGRADSTEST Test the gradients of the covariance.
% FGPLVMCREATE Create a GPLVM model with inducing varibles.
% FGPLVMDISPLAY Display an FGPLVM model.
% FGPLVMDYNAMICSFIELDPLOT 2-D field plot of the dynamics.
% FGPLVMDYNAMICSPLOT 2-D scatter plot of the latent points.
% FGPLVMDYNAMICSPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% FGPLVMDYNAMICSRUN Visualise the manifold.
% FGPLVMDYNAMICSSAMPLE Sample a field from the GP.
% FGPLVMEXPANDPARAM Expand a parameter vector into a GP-LVM model.
% FGPLVMEXTRACTPARAM Extract a parameter vector from a GP-LVM model.
% FGPLVMFIELDPLOT 2-D field plot of the dynamics.
% FGPLVMGRADIENT GP-LVM gradient wrapper.
% FGPLVMKERNDYNAMICSSAMPLE Sample a field from a given kernel.
% FGPLVMLOADRESULT Load a previously saved result.
% FGPLVMLOGLIKEGRADIENTS Compute the gradients of the EZFT sparse covariance.
% FGPLVMLOGLIKELIHOOD Log-likelihood for a GP-LVM.
% FGPLVMNEARESTNEIGHBOUR Give the number of errors in latent space for 1 nearest neighbour.
% FGPLVMOBJECTIVE Wrapper function for GPLVM objective.
% FGPLVMOPTIMISE Optimise the inducing variable based kernel.
% FGPLVMOPTIMISEPOINT Optimise the postion of a point.
% FGPLVMOPTIONS Return default options for FGPLVM model.
% FGPLVMPOINTGRADIENT Wrapper function for gradient of a single point.
% FGPLVMPOINTLOGLIKEGRADIENT Log-likelihood gradient for of a point of the GP-LVM.
% FGPLVMPOINTLOGLIKELIHOOD Log-likelihood of a point for the GP-LVM.
% FGPLVMPOINTOBJECTIVE Wrapper function for objective of a single point.
% FGPLVMPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% FGPLVMPRINTPLOT Print latent space for learnt model.
% FGPLVMREADFROMFID Load from a FID produced by the C++ implementation.
% FGPLVMREADFROMFILE Load a file produced by the c++ implementation.
% FGPLVMRESULTSDYNAMIC Load a results file and visualise them.
% FGPLVMTEST Test the gradients of the gpCovGrads function and the fgplvm models.
% FGPLVMVISUALISE Visualise the manifold.
% GPBLOCKINDICES Return indices of given block.
% GPCOMPUTEALPHA Update the vector `alpha' for computing posterior mean quickly.
% GPCOMPUTEM Compute the matrix m given the model.
% GPCOVGRADS Sparse objective function gradients wrt Covariance functions for inducing variables.
% GPCOVGRADSTEST Test the gradients of the covariance.
% GPCREATE Create a GP model with inducing varibles/pseudo-inputs.
% GPDATAINDICES Return indices of present data.
% GPDISPLAY Display a Gaussian process model.
% GPDYNAMICSCREATE Create the dynamics model. 
% GPDYNAMICSDISPLAY Display a GP dynamics model.
% GPDYNAMICSEXPANDPARAM Place the parameters vector into the model for GP dynamics.
% GPDYNAMICSEXTRACTPARAM Extract parameters from the GP dynamics model.
% GPDYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.
% GPDYNAMICSLOGLIKEGRADIENTS Gradients of the GP dynamics wrt parameters.
% GPDYNAMICSLOGLIKELIHOOD Give the log likelihood of the dynamics part.
% GPDYNAMICSSAMP Sample from the dynamics for a given input.
% GPDYNAMICSSETLATENTVALUES Set the latent values inside the model.
% GPEXPANDPARAM Expand a parameter vector into a GP model.
% GPEXTRACTPARAM Extract a parameter vector from a GP model.
% GPGRADIENT Gradient wrapper for a GP model.
% GPLOGLIKEGRADIENTS Compute the gradients for the parameters and X.
% GPLOGLIKELIHOOD Compute the log likelihood of a GP.
% GPOBJECTIVE Wrapper function for GP objective.
% GPOPTIMISE Optimise the inducing variable based kernel.
% GPOPTIONS Return default options for GP model.
% GPOUT Evaluate the output of an Gaussian process model.
% GPPOSTERIORGRADMEANVAR Gadient of the mean and variances of the posterior at points given by X.
% GPPOSTERIORMEANCOVAR Mean and covariances of the posterior at points given by X.
% GPPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% GPREVERSIBLEDYNAMICSCREATE Create the dynamics model. 
% GPREVERSIBLEDYNAMICSDISPLAY Display a GP dynamics model.
% GPREVERSIBLEDYNAMICSEXPANDPARAM Place the parameters vector into the model for GP dynamics.
% GPREVERSIBLEDYNAMICSEXTRACTPARAM Extract parameters from the GP reversible dynamics model.
% GPREVERSIBLEDYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.
% GPREVERSIBLEDYNAMICSLOGLIKEGRADIENTS Gradients of the GP reversible dynamics wrt parameters.
% GPREVERSIBLEDYNAMICSLOGLIKELIHOOD Give the log likelihood of the dynamics part.
% GPREVERSIBLEOPTIONS Return default options for GP reversible dynamics model.
% GPREVERSIBLEDYNAMICSSAMP Sample from the dynamics for a given input.
% GPREVERSIBLEDYNAMICSSETLATENTVALUES Set the latent values inside the model.
% GPSCALEBIASGRADIENT Compute the gradient of the scale and bias.
% GPUPDATEKERNELS Update the kernels that are needed.
% MODELLATENTGRADIENTS Gradients of the latent variables for dynamics models in the GPLVM.
% MODELSETLATENTVALUES Set the latent variables for dynamics models in the GPLVM.
% ROBONEDYNAMICSCREATE Create the dynamics model. 
% ROBONEDYNAMICSDISPLAY Display the robot dynamics model. 
% ROBONEDYNAMICSEXPANDPARAM Place the parameters vector into the model for first robot dynamics.
% ROBONEDYNAMICSEXTRACTPARAM Extract parameters from the robot one dynamics model.
% ROBONEDYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.
% ROBONEDYNAMICSLOGLIKEGRADIENTS Gradients of the robot one dynamics wrt parameters.
% ROBONEDYNAMICSLOGLIKELIHOOD Give the log likelihood of the robot one dynamics part.
% ROBONEDYNAMICSSETLATENTVALUES Set the latent values inside the model.
% ROBTHREEDYNAMICSCREATE Create the dynamics model. 
% ROBTHREEDYNAMICSDISPLAY Display the robot dynamics model. 
% ROBTHREEDYNAMICSEXPANDPARAM Place the parameters vector into the model for first robot dynamics.
% ROBTHREEDYNAMICSEXTRACTPARAM Extract parameters from the robot three dynamics model.
% ROBTHREEDYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.
% ROBTHREEDYNAMICSLOGLIKEGRADIENTS Gradients of the robot three dynamics wrt parameters.
% ROBTHREEDYNAMICSLOGLIKELIHOOD Give the log likelihood of the robot three dynamics part.
% ROBTHREEDYNAMICSSETLATENTVALUES Set the latent values inside the model.
% ROBTHREEDYNAMICSSETLATENTVALUES Set the latent values inside the model.
% ROBTWODYNAMICSCREATE Create the dynamics model. 
% ROBTWODYNAMICSDISPLAY Display the robot dynamics model. 
% ROBTWODYNAMICSEXPANDPARAM Place the parameters vector into the model for first robot dynamics.
% ROBTWODYNAMICSEXTRACTPARAM Extract parameters from the robot two dynamics model.
% ROBTWODYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.
% ROBTWODYNAMICSLOGLIKEGRADIENTS Gradients of the robot two dynamics wrt parameters.
% ROBTWODYNAMICSLOGLIKELIHOOD Give the log likelihood of the robot one dynamics part.
% ROBTWODYNAMICSSETLATENTVALUES Set the latent values inside the model.
