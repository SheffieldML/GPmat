function G = cmpndKernGradientK(kern, x)

% CMPNDKERNGRADIENTK Gradient of compound kernel wrt its parameters.

G = cell(1, kern.nParams);
startVal = 1;
endVal = 0;
for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved in the kernel.
    G(1, startVal:endVal)  = kernGradientK(kern.comp{i}, ...
                                           x(:, kern.comp{i}.index));
    % See note below on these 4 lines.
    factors = kernFactors(kern.comp{i}, 'gradfact');
    for j = 1:length(factors.index);
      G{startVal -1 + factors.index(j)} = G{startVal - 1+factors.index(j)}*factors.val(j);
    end

  else
    % all the data is involved with the kernel.
    G(1, startVal:endVal)  = kernGradientK(kern.comp{i}, x);
    % See note below on these 4 lines.
    factors = kernFactors(kern.comp{i}, 'gradfact');
    for j = 1:length(factors.index);
      G{startVal -1 + factors.index(j)} = G{startVal - 1+factors.index(j)}*factors.val(j);
    end
  end
  startVal = endVal + 1;
end

% Note on Gradfactor lines.
% When parameters are being optimised in a transformed space, we need to
% account for this with their gradients. Normally this is done by
% checking the transform field of the kernel. However the compound kernel
% is implemented without a transform field because these factors are
% accounted for when computing the component kernels. TO account for
% gradient factors we need to include the following.
% Check if parameters are being optimised in a transformed space.
