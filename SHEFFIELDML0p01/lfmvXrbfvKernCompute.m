function K = lfmvXrbfvKernCompute(lfmKern, rbfKern, t1, t2)

% LFMVXRBFVKERNCOMPUTE Compute a cross kernel between the LFMV and RBFV kernels.
%
%	Description:
%
%	K = LFMVXRBFVKERNCOMPUTE(LFMKERN, RBFKERN, T) computes cross kernel
%	terms between LFMV and RBFV kernels for the multiple output kernel. 
%	This function is employed in the SDLFM kernel to compute derivatives
%	with respect to the switching points.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  LFMKERN - the kernel structure associated with the LFMV kernel.
%	  RBFKERN - the kernel structure associated with the RBFV kernel.
%	  T - inputs for which kernel is to be computed.
%
%	K = LFMVXRBFVKERNCOMPUTE(LFMKERN, RBFKERN, T1, T2) computes cross
%	kernel terms between LFMV and RBF kernels for the multiple output
%	kernel.  This function is employed in the SDLFM kernel to compute
%	derivatives with respect to the switching points.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  LFMKERN - the kernel structure associated with the LFMV kernel.
%	  RBFKERN - the kernel structure associated with the RBFV kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.


%	Copyright (c) 2010 Mauricio Alvarez


if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 || size(t2, 2) > 1
  error('Input can only have one column');
end


if lfmKern.inverseWidth ~= rbfKern.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end
  
% Get length scale out.
sigma2 = 2/lfmKern.inverseWidth;
sigma = sqrt(sigma2);

% Parameters of the kernel
alpha = lfmKern.damper./(2*lfmKern.mass);
omega = sqrt(lfmKern.spring./lfmKern.mass - alpha.*alpha);

gamma1 = alpha + j*omega;
gamma2 = alpha - j*omega;

sK = lfmvvComputeUpsilonMatrix(gamma2,sigma2,t1, t2, 1) - ...
    lfmvvComputeUpsilonMatrix(gamma1,sigma2,t1, t2, 1);

if lfmKern.isNormalised
   K0 = lfmKern.sensitivity/(j*4*sqrt(2)*lfmKern.mass*omega);
else
   K0 = sqrt(pi)*sigma*lfmKern.sensitivity/(j*4*lfmKern.mass*omega);
end
    
K = K0*sK;
