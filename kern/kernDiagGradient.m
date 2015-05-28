function g = kernDiagGradient(kern, x, covDiag)

% KERNDIAGGRADIENT Compute the gradient of the kernel's parameters for the diagonal.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% kernel matrix with respect to the parameters of the kernel. The 
% parameters' gradients are returned in the order given by the 
% kernExtractParam command.
% ARG kern : the kernel structure for which the gradients are computed.
% ARG x : the input data for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest
% with respect to the diagonal elements of the kernel matrix.
% RETURN g : gradients of the relevant function with respect to each of the parameters. Ordering should match the ordering given in kernExtractParam.
%
% SEEALSO : kernDiagGradient, kernExtractParam, kernGradient

% KERN


fileName = [kern.type 'KernDiagGradient'];
if exist(fileName) == 2
  fhandle = str2func(fileName);
  g = fhandle(kern, x, covDiag);
else
  fhandle = str2func([kern.type 'KernGradient']);
  g = zeros(1, kern.nParams);
  for i = 1:size(x, 1)
    g = g ...
        + fhandle(kern, x(i, :), covDiag(i));
  end
end
% Check if parameters are being optimised in a transformed space.
factors = kernFactors(kern, 'gradfact');
g(factors.index) = g(factors.index).*factors.val;

  
