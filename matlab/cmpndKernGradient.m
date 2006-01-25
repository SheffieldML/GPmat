function g = cmpndKernGradient(kern, x, varargin)

% CMPNDKERNGRADIENT Gradient of compound kernel's parameters.

% KERN

g = zeros(1, kern.nParams);
startVal = 1;
endVal = 0;

for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved in the kernel.
    if nargin < 4
      g(1, startVal:endVal)  = kernGradient(kern.comp{i}, ...
                                            x(:, kern.comp{i}.index), ...
                                            varargin{:});
    else
      g(1, startVal:endVal) = kernGradient(kern.comp{i}, ...
                                            x(:, kern.comp{i}.index), ...
                                           varargin{1}(:, kern.comp{i}.index), ...
                                           varargin{2:end});
    end
  else
    % all the data is involved with the kernel.
    g(1, startVal:endVal)  = kernGradient(kern.comp{i}, x, varargin{:});
  end
  startVal = endVal + 1;
end
g = g*kern.paramGroups;