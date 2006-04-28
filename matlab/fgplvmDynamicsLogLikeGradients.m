function gX = fgplvmDynamicsLogLikeGradients(model)

% FGPLVMDYNAMICSLOGLIKEGRADIENTS Gradients of the dynamics portion.


switch model.approx
 case {'ftc', 'dtc', 'fitc', 'pitc','nftc'}
  gKX = kernGradX(model.dynamics.kern, model.X(1:end-1, :), ...
                  model.X(1:end-1, :));
  gKX = gKX*2;
  dgKX = kernDiagGradX(model.kern, model.X(1:end-1, :));
  for i = 1:model.N-1
    gKX(i, :, i) = dgKX(i, :);
  end
  gX = zeros(model.N, model.q);
  for k = 1:model.q
    gK = fgplvmDynamicsCovarianceGradients(model, k);
    for i = 1:model.N-1
      for j = 1:model.q
       gX(i, j) = gX(i, j) + gKX(:, j, i)'*gK(:, i);
      end
    end
  end
  gX(2:end, :) = gX(2:end, :) - model.dynamics.invK*model.X(2:end, :);

end  
  
function gK = fgplvmDynamicsCovarianceGradients(model, dimension)

% FGPLVMDYNAMICSCOVARIANCEGRADIENTS

invKy = model.dynamics.invK*model.X(2:end, dimension);
gK = -model.dynamics.invK + invKy*invKy';
gK = gK*.5;

  