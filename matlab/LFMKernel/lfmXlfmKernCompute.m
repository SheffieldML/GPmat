function K = lfmXlfmKernCompute(lfmKern1, lfmKern2, t, t2);
  
% LFMXLFMKERNCOMPUTE Compute a cross kernel between two LFM kernels.
% FORMAT
% DESC computes cross kernel terms between two LFM kernels for
% the multiple output kernel. 
% ARG lfmKern1 : the kernel structure associated with the first LFM
% kernel.
% ARG lfmKern2 : the kernel structure associated with the second LFM
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two LFM kernels for
% the multiple output kernel. 
% ARG lfmKern1 : the kernel structure associated with the first LFM
% kernel.
% ARG lfmKern2 : the kernel structure associated with the second LFM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : lfmKernParamInit, lfmKernCompute, lfmKernParamInit
%
% COPYRIGHT : David Luengo, 2007
%  
% MODIFICATIONS : Neil D. Lawrence, 2007; David Luengo, 2008

% LFM
  
if nargin < 4
  t2 = t;
end
if size(t, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if lfmKern1.inverseWidth ~= lfmKern2.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end

% Get length scale out.
sigma2 = 2/lfmKern1.inverseWidth;
sigma = sqrt(sigma2);

% Parameters of the kernel
alpha(1) = lfmKern1.damper./(2*lfmKern1.mass);
alpha(2) = lfmKern2.damper./(2*lfmKern2.mass);
omega(1) = sqrt(lfmKern1.spring./lfmKern1.mass - alpha(1)*alpha(1));
omega(2) = sqrt(lfmKern2.spring./lfmKern2.mass - alpha(2)*alpha(2));

% NEED TO DEAL WITH CRITICALLY DAMPED OMEGA = 0!!! 
% Kernel evaluation
if isreal(omega)
  gamma1 = alpha(1) + j*omega(1);
  gamma2 = alpha(2) + j*omega(2);
  
  K = sigma*lfmKern1.sensitivity*lfmKern2.sensitivity* ...
        real(lfmComputeH(conj(gamma2),gamma1,sigma2,t,t2) + ...
             lfmComputeH(gamma1,conj(gamma2),sigma2,t2,t) - ...
             lfmComputeH(gamma2,gamma1,sigma2,t,t2) - ...
             lfmComputeH(gamma1,gamma2,sigma2,t2,t));
  K = K*sqrt(pi)/(4*lfmKern1.mass*lfmKern2.mass*prod(omega));
else
  gamma1_p = alpha(1) + j*omega(1);
  gamma1_m = alpha(1) - j*omega(1);
  gamma2_p = alpha(2) + j*omega(2);
  gamma2_m = alpha(2) - j*omega(2);

  K = sigma*lfmKern1.sensitivity*lfmKern2.sensitivity* ...
        (lfmComputeH(gamma2_m,gamma1_p,sigma2,t,t2) + ...
         lfmComputeH(gamma1_p,gamma2_m,sigma2,t2,t) + ...
         lfmComputeH(gamma2_p,gamma1_m,sigma2,t,t2) + ...
         lfmComputeH(gamma1_m,gamma2_p,sigma2,t2,t) - ...
         lfmComputeH(gamma2_m,gamma1_m,sigma2,t,t2) - ...
         lfmComputeH(gamma1_m,gamma2_m,sigma2,t2,t) - ...
         lfmComputeH(gamma2_p,gamma1_p,sigma2,t,t2) - ...
         lfmComputeH(gamma1_p,gamma2_p,sigma2,t2,t));
  K = K*sqrt(pi)/(8*lfmKern1.mass*lfmKern2.mass*prod(omega));
end

K = real(K); % introduced Mauricio Alvarez 2008

