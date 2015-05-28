function [dWdB, dWdD, dWdS, dWdalpha, dWdgParam, dWdn] = gpsimMapWGradient(model, ...
                                                  k)
% GPSIMMAPWGRADIENT Compute the gradients of W with respect to the parameters of the k-th gene.
% FORMAT
% DESC computes the gradients of W with respect to the parameters of the k-th gene given
% Gaussian process for use in a single input motif protein network.
% ARG model : the model for which the log likelihood is computed.
% ARG k : the k-th gene.  
% RETURN dWdB : the gradients of W w.r.t the model paramter Bk.
% RETURN dWdD : the gradients of W w.r.t the model paramter Dk.
% RETURN dWdS : the gradients of W w.r.t the model paramter Sk.
% RETURN dWdalpha : the gradients of W w.r.t the model paramter alpha_k.
% RETURN dWdgParam : the gradients of W w.r.t the model paramter gamma_k
% RETURN dWdn : the gradients of W w.r.t the general noise variance of the k-gene.
%
% SEEALSO : gpsimMapCreate, gpsimMapLogLikelihood,
% gpsimMapGradient, gpsimMapFunctionalLogLikeGradients
%
% COPYRIGHT : Pei Gao, Magnus Rattray and Neil D. Lawrence, 2008

% GPSIM  
  
intPoints = model.times_index(1)+1:(model.numMapPts);
step2 = model.step*model.step;
S2 = model.S.*model.S;
numData = length(model.t);

[w1, w2] = size(model.W);

dWdB = zeros(w1, w2);
dWdD = zeros(w1, w2);
dWdS = zeros(w1, w2);
dWdalpha = [];
dWdgParam = [];
dWdn = [];
dWdalpha = [];
if model.ngParam > 0
  ngParamk = model.ngParam/model.numGenes;
  dWdgParam= zeros(w1, w2, ngParamk);
  gInd = k;
else
  gInd = 1;
end

% check if it's multiple g(f).
if isfield(model, 'isGroupNonlinearity') && strcmp(model.nonLinearity{k}, ...
                                                   'repression')
  dWdalpha = zeros(w1, w2);
end


if isfield(model, 'includeNoise') && model.includeNoise
  noiseMat = ones(numData, 1)*model.noiseVar;
  yvar = model.yvar + noiseMat;
  dWdn = zeros(w1, w2);
else
  yvar = model.yvar;
end

for p = intPoints
    for i=1:numData
      arg = model.t(i)-model.mapt(p);
      if arg >= 0
        ind = i + (k-1)*numData;
        beta_ik=1/yvar(ind);
        
        [dxdB dxdD dxdS dxdalpha dxdgParam] = gpsimXGradient(model, i, k);

        dWdB(p, p)=dWdB(p, p)+beta_ik*model.g_grad2(p,gInd)*dxdB* ...
                exp(-model.D(k)*arg+log(model.S(k)) +log(model.step));
        
        if isfield(model, 'isGroupNonlinearity') 
          if strcmp(model.nonLinearity{k}, 'repression')
            dWdalpha(p,p) = dWdalpha(p, p)+beta_ik*model.g_grad2(p,gInd)* ...
                dxdalpha*exp(-model.D(k)*arg+log(model.S(k)) +log(model.step));
          end
        end
        
        factor = model.ypred(model.times_index(i), k)-model.y(ind);

        dWdD(p, p) = dWdD(p, p)+model.step*beta_ik*model.g_grad2(p,gInd)* ...
                (dxdD-factor*arg)*exp(-model.D(k)*arg+log(model.S(k))) ;

        dWdS(p, p) = dWdS(p, p)+model.step*beta_ik*exp(-model.D(k)* ...
                           arg)*(dxdS*model.S(k)*model.g_grad2(p,gInd)+ ...
                                 factor*model.g_grad2(p,gInd));

        if model.ngParam > 0
          for gParamInd = 1:ngParamk
            dWdgParam(p, p, gParamInd)= dWdgParam(p, p, gParamInd)+ ...
                model.step*beta_ik*exp(-model.D(k)*arg)* ...
                (dxdgParam(gParamInd)*model.S(k)*model.g_grad2(p, ...
                gInd)+ factor*model.S(k)*model.dggrad2(p,gInd));
          end
        end
        
        if isfield(model,'includeNoise') && model.includeNoise
          dWdn(p, p) = dWdn(p, p)-2*sqrt(model.noiseVar(k))* ...
              beta_ik^2*factor*model.g_grad2(p,gInd)*exp(-model.D(k)*arg+ ...
                                                         log(model.step)+ ...
                                                         log(model.S(k)));
        end
      end
    end
end

for p = intPoints
  for q = intPoints
    for i = 1:numData
      arg1 = model.t(i)-model.mapt(p);
      arg2 = model.t(i)-model.mapt(q);
      if arg1 >= 0 && arg2 >= 0
        ind = i + (k-1)*numData;
        beta_ik = 1/yvar(ind);

        dWdD(p,q) = dWdD(p,q)-beta_ik*model.g_grad(p,gInd)* ...
            model.g_grad(q,gInd)*(arg1+arg2)*exp(-model.D(k)*(arg1+arg2)+ ...
                                                 log(S2(k))+log(step2));

        dWdS(p,q) = dWdS(p,q)+2*beta_ik*model.S(k)* model.g_grad(q,gInd)* ...
            model.g_grad(p,gInd) * exp(-model.D(k)*(arg1+arg2)+log(step2));

        if model.ngParam > 0
          for gParamInd = 1:ngParamk
            dWdgParam(p, q, gParamInd)= dWdgParam(p, q, gParamInd)+beta_ik*(model.dggrad(q,gInd)*model.g_grad(p,gInd)+model.g_grad(q,gInd)*model.dggrad(p,gInd))*exp(-model.D(k)*(arg1+arg2)+log(S2(k))+log(step2));
          end
        end
        
        if isfield(model, 'includeNoise') && model.includeNoise
          dWdn(p, q) = dWdn(p, q)-2*sqrt(model.noiseVar(k))* ...
              beta_ik^2*model.g_grad(p,gInd)*model.g_grad(q,gInd)*exp(- ...
                    model.D(k)*(arg1+arg2)+log(step2)+log(S2(k)));
        end
      end
    end
  end
end




