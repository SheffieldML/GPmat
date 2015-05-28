function [dfuncGrad] = gpsimMapFunctionalLikeGrad2(model, p)
  
% GPSIMMAPFUNCTIONALLIKEGRAD2 Compute the functional gradient for GPSIMMAP.
% FORMAT
% DESC computes the functional gradient of the log likelihood for
% use in the MAP approximation to the GPSIM posterior solution.
% ARG model : the model for which the gradient is to be computed.
% RETURN p: the gradient of the log likelihood with respect to the
% points of the function.
% 
% SEEALSO : gpsimMapLikeGradientImplicit
%
% COPYRIGHT : Pei Gao, Magnus Rattray and Neil D. Lawrence, 2008

% SHEFFIELDML
  
dfuncGrad = zeros(1, model.numParams-model.kern.nParams);
gB = zeros(model.numGenes, 1);
gD = zeros(model.numGenes, 1);
gS = zeros(model.numGenes, 1);
ngParamk = model.ngParam/model.numGenes;

if isfield(model,'includeNoise') && model.includeNoise
  gn = zeros(model.numGenes, 1);
  noiseMat = ones(length(model.t),1)*model.noiseVar;
  yvar = model.yvar + noiseMat;
else
  yvar = model.yvar;
end

if isfield(model, 'includeRepression') && model.includeRepression
  galpha = zeros(model.numGenes, 1);
  nBaseParams = 4;
else
  nBaseParams = 3;
end
  
if model.ngParam > 0
  ggParam= zeros(model.numGenes, ngParamk);
end  

numData = length(model.t);
for k = 1:model.numGenes
  if model.ngParam
    gInd = k;
  else
    gInd = 1;
  end
  for i=1:numData
    arg = model.t(i)-model.mapt(p);
    if arg >= 0
      ind = i + (k-1)*numData;
      beta_ik=1/yvar(ind);

      [dxdB dxdD dxdS dxdalpha dxdgParam] = gpsimXGradient(model, i, k);
        
      gB(k) = gB(k) - beta_ik*model.g_grad(p,gInd)*dxdB*exp(-model.D(k)*arg ...
                               +log(model.S(k)) +log(model.step));

      factor = model.ypred(model.times_index(i), k) - model.y(ind);
      gD(k) = gD(k) - beta_ik*model.g_grad(p,gInd)* (dxdD-arg*factor) ...
              *exp(-model.D(k)*arg+log(model.S(k)) +log(model.step));
       
      gS(k) = gS(k) - beta_ik*model.g_grad(p,gInd)* (dxdS*model.S(k) + ...
                       factor) *exp(-model.D(k)*arg+log(model.step));
      
      if isfield(model, 'includeRepression') && model.includeRepression      
        galpha(k) = galpha(k) - beta_ik*model.g_grad(p,gInd)*dxdalpha ...
            *exp(-model.D(k)*arg+log(model.S(k)) +log(model.step));
      end
        
      if model.ngParam > 0
        ggParam(k,:) = ggParam(k,:) - beta_ik* (model.g_grad(p,gInd)* ...
            dxdgParam+factor*model.dggrad(p, gInd)) *exp(-model.D(k)* ...
                                arg+log(model.S(k))+log(model.step));
      end
      
      if isfield(model,'includeNoise') && model.includeNoise
        gn(k) = gn(k) + 2*sqrt(model.noiseVar(k))*beta_ik^2* ...
                model.g_grad(p,gInd)*model.S(k)*factor*exp(-model.D(k)*arg+log(model.step));
      end
    end
  end
end
  
geneParam = nBaseParams + ngParamk;  
startPoint = 1;
endPoint = geneParam;
  
for k = 1:model.numGenes
  dfuncGrad(startPoint:(startPoint+2)) = [gB(k) gS(k) gD(k)];
  
  if isfield(model, 'includeRepression') && model.includeRepression
    dfuncGrad(startPoint+3) = galpha(k);
  end
  
  if model.ngParam > 0
    dfuncGrad((startPoint+nBaseParams):endPoint) = ggParam(k,:)';
  end    
  startPoint = startPoint + geneParam;
  endPoint = endPoint + geneParam;
end

if isfield(model,'includeNoise') && model.includeNoise
  dfuncGrad(end-length(gn)+1:end) = gn';
end
  
  
