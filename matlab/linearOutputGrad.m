function g = linearOutputGrad(model, X)

% LINEAROUTPUTGRAD Evaluate derivatives of linear model outputs with respect to parameters.

% MLTOOLS

numData = size(X, 1);
g(:, :, 1) = [X zeros(size(X)) ones(numData, 1) ...
                    zeros(numData, 1)];
g(:, :, 2) = [zeros(size(X)) X zeros(numData, 1) ...
                    ones(numData, 1)];
