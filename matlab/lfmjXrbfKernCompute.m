function K = lfmjXrbfKernCompute(lfmKern, rbfKern, t1, t2)

% LFMJXRBFKERNCOMPUTE Compute cross kernel between the LFMJ and RBF kernels.
% FORMAT
% DESC computes cross kernel terms between LFMJ and RBF kernels for
% the multiple output kernel. 
% ARG lfmKern : the kernel structure associated with the LFMJ
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between LFMJ and RBF kernels for
% the multiple output kernel. 
% ARG lfmKern : the kernel structure associated with the LFMJ
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% COPYRIGHT : Mauricio Alvarez, 2010

% KERN

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

sK = lfmjpComputeUpsilonMatrix(gamma2,sigma2,t1, t2, 0) - ...
    lfmjpComputeUpsilonMatrix(gamma1,sigma2,t1, t2, 0);

if lfmKern.isNormalised
   K0 = lfmKern.sensitivity/(j*4*sqrt(2)*lfmKern.mass*omega);
else
   K0 = sqrt(pi)*sigma*lfmKern.sensitivity/(j*4*lfmKern.mass*omega);
end

K = K0*sK;
