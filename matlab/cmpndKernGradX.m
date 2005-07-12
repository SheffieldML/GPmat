function gX = cmpndKernGradX(kern, X, X2)

% CMPNDKERNGRADX Gradient of compound kernel with respect to a point X.

% KERN

i = 1;
fhandle = str2func([kern.comp{i}.type 'KernGradX']);

if ~isempty(kern.comp{i}.index)
  % only part of the data is involved with the kernel.

  gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
  gX(:, kern.comp{i}.index, :) = fhandle(kern.comp{i}, ...
                                      X(:, kern.comp{i}.index), ...
                                      X2(:, kern.comp{i}.index));
else
  % all the data is involved with the kernel.
  gX = fhandle(kern.comp{i}, X, X2);
end
for i = 2:length(kern.comp)
  fhandle = str2func([kern.comp{i}.type 'KernGradX']);
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved with the kernel.
    gX(:, kern.comp{i}.index, :) = ...
        gX(:, kern.comp{i}.index, :) + ...
        fhandle(kern.comp{i}, ...
                X(:, kern.comp{i}.index), ...
                X2(:, kern.comp{i}.index));
  else
    % all the data is involved with the kernel.
    gX = gX + fhandle(kern.comp{i}, X, X2);
  end
end
