function K = lfmXrbfKernCompute(lfmKern, rbfKern, t1, t2)

% LFMXRBFKERNCOMPUTE Compute a cross kernel between the LFM and RBF kernels.
% FORMAT
% DESC computes cross kernel terms between LFM and RBF kernels for
% the multiple output kernel. 
% ARG lfmKern : the kernel structure associated with the LFM
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between LFM and RBF kernels for
% the multiple output kernel. 
% ARG lfmKern : the kernel structure associated with the LFM
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, lfmKernParamInit, rbfKernParamInit
%
% COPYRIGHT : David Luengo, 2007, 2008, Mauricio Alvarez, 2008
%
% MODIFICATIONS : Neil D. Lawrence, 2007


% KERN

if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if lfmKern.inverseWidth ~= rbfKern.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end
  
if lfmKern.isNormalised ~= rbfKern.isNormalised
  error('Kernels cannot be cross combined if they have different normalization settings.')
end

% Get length scale out.
sigma2 = 2/lfmKern.inverseWidth;
sigma = sqrt(sigma2);

% Parameters of the kernel
alpha = lfmKern.damper./(2*lfmKern.mass);
omega = sqrt(lfmKern.spring./lfmKern.mass - alpha.*alpha);

% Kernel evaluation
if isreal(omega)
    gamma = alpha + j*omega;
    sK = imag(lfmComputeUpsilonMatrix(gamma,sigma2, t1,t2));
    if lfmKern.isNormalised 
        K0 = (lfmKern.sensitivity/(2*sqrt(2)*lfmKern.mass*omega));
    else
        K0 = (sqrt(pi)*sigma*lfmKern.sensitivity/(2*lfmKern.mass*omega));
    end    
    K = -K0*sK;    
else
    gamma1 = alpha + j*omega;
    gamma2 = alpha - j*omega;
    sK = lfmComputeUpsilonMatrix(gamma2,sigma2,t1, t2) - ...
          lfmComputeUpsilonMatrix(gamma1,sigma2,t1, t2);
    if lfmKern.isNormalised
        K0 = lfmKern.sensitivity/(j*4*sqrt(2)*lfmKern.mass*omega);
    else
        K0 = sqrt(pi)*sigma*lfmKern.sensitivity/(j*4*lfmKern.mass*omega);
    end
    K = K0*sK;
end


