function g = cmpndKernDiagGradient(kern, x, covDiag)


% CMPNDKERNDIAGGRADIENT Compute the gradient of the CMPND kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% compound kernel matrix with respect to the parameters of the kernel. The
% parameters' gradients are returned in the order given by the
% cmpndKernExtractParam command.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG x : the input data for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% cmpndKernExtractParam.
%
% SEEALSO : cmpndKernParamInit, kernDiagGradient, cmpndKernExtractParam, cmpndKernGradient
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


g = zeros(1, kern.nParams);
startVal = 1;
endVal = 0;

for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved in the kernel.
    g(1, startVal:endVal) = kernDiagGradient(kern.comp{i}, ...
                                             x(:, kern.comp{i}.index), ...
                                             covDiag);
  else
    % all the data is involved with the kernel.
    g(1, startVal:endVal) = kernDiagGradient(kern.comp{i}, x, covDiag);
  end
  startVal = endVal + 1;
end
g = g*kern.paramGroups;
