function [g, gdata] = gpsimMapFunctionalLogLikeGradients(model)

% GPSIMMAPFUNCTIONALLOGLIKEGRADIENTS Compute the functional gradient for GPSIMMAP.
% FORMAT
% DESC computes the functional gradient of the log likelihood for
% use in the MAP approximation to the GPSIM posterior solution.
% ARG model : the model for which the gradient is to be computed.
% RETURN g : the gradient of the log likelihood with respect to the
% points of the function.
% 
% SEEALSO : gpsimMapCreate, gpsimMapUpdateYpred
%
% COPYRIGHT : Magnus Rattray and Neil D. Lawrence, 2006

gdata = zeros(1, model.numMapPts);
numData = length(model.t);
for k=model.times_index(1)+1:(model.numMapPts)
  temp=0;
  for i=1:numData
    arg = model.t(i)-model.mapt(k);
    if arg >= 0 
      for j=1:model.numGenes
        ind = (i + (j-1)*numData);
        beta_ij=1/model.yvar(ind);
        factor=(model.ypred(model.times_index(i), j)-model.y(ind))*beta_ij;
        temp=temp+factor*model.g_grad(k)*exp(-model.D(j)*arg)*model.S(j);
        %disp(factor)
      end
    end
  end
  gdata(k) = -temp*model.step;
end

% Add term from prior.
g = gdata - (model.invK*model.f)';