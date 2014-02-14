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
%  
% MODIFIED : Pei Gao, 2008

% SHEFFIELDML

gdata = zeros(1, model.numMapPts);
numData = length(model.t);

if isfield(model,'includeNoise') && model.includeNoise
  noiseMat = ones(numData, 1)*model.noiseVar;
  yvar = model.yvar + noiseMat;
else
  yvar = model.yvar;
end

for k=model.times_index(1)+1:(model.numMapPts)
  temp=0;
  for i=1:numData
    arg = model.t(i)-model.mapt(k);
    if arg >= 0 
      for j=1:model.numGenes
        if model.ngParam
          gInd = j;
        else
          gInd = 1;
        end
        ind = (i + (j-1)*numData);
        beta_ij=1/yvar(ind);
        factor=(model.ypred(model.times_index(i), j)-model.y(ind))*beta_ij;
        temp=temp+factor*model.g_grad(k,gInd)*exp(-model.D(j)*arg)*model.S(j);
        %disp(factor)
      end
    end
  end
  gdata(k) = -temp*model.step;
end

% Add term from prior.
g = gdata - (model.invK*model.f)';

% Add constraints
gCons = zeros(size(g));
if isfield(model, 'priorProtein') && ~isempty(model.priorProtein)
  nCons = length(model.priorProtein);
  for i=1:nCons
    ftimeIndex = find((model.priorProteinTimes(i)-model.mapt)==0);
    g(ftimeIndex) = g(ftimeIndex)-model.consLambda*(model.f(ftimeIndex)- ...
                                          model.priorProtein(i));
  end
end

