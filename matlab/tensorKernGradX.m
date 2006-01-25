function gX = tensorKernGradX(kern, X, X2)

% CMPNDKERNGRADX Gradient of compound kernel with respect to a point X.

% KERN

if nargin > 2
  twoXin = 1;
else
  twoXin = 0;
end
i = 1;

% Create a kernel with component i missing.
tempKern = tensorKernSlash(kern, i);
% Compute kernel matrix from that kernel.
Kslash = kernCompute(tempKern, X, X2);


if ~isempty(kern.comp{i}.index)
  % only part of the data is involved with the kernel.

  gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
  gX(:, kern.comp{i}.index, :) = kernGradX(kern.comp{i}, ...
                                      X(:, kern.comp{i}.index), ...
                                      X2(:, kern.comp{i}.index));
  gX(:, kern.comp{i}.index, :) = ...
      gX(:, kern.comp{i}.index, :) ...
      .*repmat(reshape(Kslash, [size(X2, 1) 1 size(X, 1)]), ...
               [1 length(kern.comp{i}.index) 1]);
else
  % all the data is involved with the kernel.
  gX = kernGradX(kern.comp{i}, X, X2);

  gX = gX.*repmat(reshape(Kslash, ...
                          [size(X2, 1) 1 size(X, 1)]), ...
                  [1 size(X2, 2) 1]);
  
end
for i = 2:length(kern.comp)
  % Create a kernel with component i missing.
  tempKern = tensorKernSlash(kern, i);
  % Compute kernel matrix from that kernel.
  Kslash = kernCompute(tempKern, X, X2);
  
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved with the kernel.
    gX(:, kern.comp{i}.index, :) = ...
        gX(:, kern.comp{i}.index, :) + ...
        kernGradX(kern.comp{i}, ...
                  X(:, kern.comp{i}.index), ...
                  X2(:, kern.comp{i}.index))   ...
        .*repmat(reshape(Kslash, [size(X2, 1) 1 size(X, 1)]), ...
                 [1 length(kern.comp{i}.index) 1]);
  else
    % all the data is involved with the kernel.
    gX = gX + kernGradX(kern.comp{i}, X, X2).* ...
         repmat(reshape(Kslash, ...
                        [size(X2, 1) 1 size(X, 1)]), ...
                [1 size(X2, 2) 1]);
  end
end
