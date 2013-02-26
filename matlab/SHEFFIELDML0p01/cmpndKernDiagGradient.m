function g = cmpndKernDiagGradient(kern, x, covDiag)

% CMPNDKERNDIAGGRADIENT Compute the gradient of the CMPND kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = CMPNDKERNDIAGGRADIENT(KERN, X, FACTORS) computes the gradient of
%	functions of the diagonal of the compound kernel matrix with respect
%	to the parameters of the kernel. The parameters' gradients are
%	returned in the order given by the cmpndKernExtractParam command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   cmpndKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	CMPNDKERNPARAMINIT, KERNDIAGGRADIENT, CMPNDKERNEXTRACTPARAM, CMPNDKERNGRADIENT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



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
