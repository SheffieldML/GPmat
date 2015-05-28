function kern = lmcKernParamInit(kern)

% LMCKERNPARAMINIT LMC kernel parameter initialisation.
% The Linear Model of Coregionalization kernel contains the structure to
% compute the kernel matrix of a multiple-output systems that employs the
% Linear Model of Coregionalization to model the covariance function. In
% detail, if the covariance for the outputs is K, then K is expressed as a
% sum of Kronecker products 
%
%   K = \sum_{q=1}^Q B_q \otimes K_q,
%
% where B_q is the coregionalization matrix q, K_q is a basic covariance
% associated to a basic process q and \otimes represents the Kronecker
% product. In this kernel, we restrict Q=1. This is also known as the
% intrinsic coregionalization model. For Q>1, we globally use the multigp
% toolbox. LMC covariance contains two parts: the coregionalization 
% matrix and the covariance of the basic process. For the coregionalization 
% matrix we specify the rank of the matrix (default is 1, which is 
% equivalent to the semiparametric latent factor model) and for the basic 
% covariance process, specify the type of basic covariance (default is 
% 'gaussian').
% FORMAT
% DESC initialises the linear model of coregionalization kernel structure  
% with some default parameters for the initial conditions.
% ARG kern : the kernel structure which requires initialisation.
% ARG options : options for the switching dynamical structure.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, multigpCreate
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if isfield(kern, 'options') && isfield(kern.options, 'basicKernelType')
    kern.basicKernelType = kern.options.basicKernelType;
else
    kern.basicKernelType = 'gaussian';
end

if isfield(kern, 'options') && isfield(kern.options, 'rankCorregMatrix')
   kern.rankCorregMatrix = kern.options.rankCorregMatrix; 
else
   kern.rankCorregMatrix = 1; 
end

fhandle = str2func([kern.basicKernelType 'KernParamInit']);
if isfield(kern, 'options') && isfield(kern.options, 'isArd')
   kern = fhandle(kern, kern.options.isArd); 
else
   kern = fhandle(kern);
end

if isfield(kern, 'options') && isfield(kern.options, 'nout')
    kern.nout = kern.options.nout;
    kern.A = rand(kern.nout, kern.rankCorregMatrix); 
else
    % If the number of outputs is not provided, we assume is one
    kern.nout = 1;
    kern.A = rand(1, kern.rankCorregMatrix);
end
if kern.nout<kern.rankCorregMatrix
    error('The rank of the matrix is greater than the number of outputs')
end
kern.B = kern.A*kern.A';
% Count again the number of parameters
kern.nParamsBK = kern.nParams;
kern.nParams = kern.nParams + kern.rankCorregMatrix*kern.nout;
















