function g = tensorKernGradient(kern, x, varargin)

% TENSORKERNGRADIENT Gradient of compound kernel's parameters.

% KERN


% Last of varargin is covGrad.
g = zeros(1, kern.nParams);
startVal = 1;
endVal = 0;
twoXin = 0;
if length(varargin) > 1
  twoXin = 1;
  x2 = varargin{1};
  covGrad = varargin{end};
else
  covGrad = varargin{end};
end

for i = 1:length(kern.comp)
  tempKern = tensorKernSlash(kern, i);
  if twoXin
    tempCovGrad = covGrad.*kernCompute(tempKern, x, x2);
  else
    tempCovGrad = covGrad.*kernCompute(tempKern, x);
  end
  endVal = endVal + kern.comp{i}.nParams;
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved in the kernel.
    if ~twoXin
      g(1, startVal:endVal)  = kernGradient(kern.comp{i}, ...
                                            x(:, kern.comp{i}.index), ...
                                            tempCovGrad);
    else
      g(1, startVal:endVal) = kernGradient(kern.comp{i}, ...
                                           x(:, kern.comp{i}.index), ...
                                           x2(:, kern.comp{i}.index), ...
                                           tempCovGrad);
    end
  else
    if ~twoXin
      % all the data is involved with the kernel.
      g(1, startVal:endVal)  = kernGradient(kern.comp{i}, x, ...
                                            tempCovGrad);
    else
      g(1, startVal:endVal) = kernGradient(kern.comp{i}, x, x2, ...
                                           tempCovGrad);
    end
  end
  startVal = endVal + 1;
end
g = g*kern.paramGroups;