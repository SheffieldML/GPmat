function dWdf =  gpsimMapFunctionalWGradient(model, l)

% GPSIMMAPFUNCTIONALWGRADIENT computes the implicit part of the gradient

% FORMAT
% DESC updates the portion of the Hessian associated with the
% likelihood of the data.
% ARG model : the model for which the data component of the Hessian is
% to be updated.
% RETURN model : the model with the data component of the Hessian up
% to date.
%
% SEEALSO : gpsimMapCreate, gpsimMapFunctionalGradient, gpsimMapUpdatePosteriorCovariance
%
% COPYRIGHT : Pei Gao, Magnus Rattray and Neil D. Lawrence, 2008
  
% SHEFFIELDML

intPoints = model.times_index(1)+1:(model.numMapPts);
step2 = model.step*model.step;
S2 = model.S.*model.S;
numData = length(model.t);

dWdf = zeros(size(model.W));

if isfield(model,'includeNoise') && model.includeNoise
  noiseMat = ones(numData, 1)*model.noiseVar;
  yvar = model.yvar + noiseMat;
else
  yvar = model.yvar;
end

for p = intPoints
    % part1: p=l or q=l
    for q = intPoints
        for i = 1:numData
            arg1 = model.t(i)-model.mapt(p);
            arg2 = model.t(i)-model.mapt(q);
            if arg1 >= 0 && arg2 >= 0
                if p==l
                    for j = 1:model.numGenes
                      if model.ngParam
                        gInd = j;
                      else
                        gInd = 1;
                      end
                      ind = i + (j-1)*numData;
                      beta_ij = 1/yvar(ind);
                      dWdf(p,q) = dWdf(p,q)+beta_ij*model.g_grad2(p,gInd)*model.g_grad(q,gInd)*exp(-model.D(j)*(arg1+arg2) ...
                            +log(S2(j))+log(step2));
                    end
                end
                if q==l
                    for j = 1:model.numGenes
                      if model.ngParam
                        gInd = j;
                      else
                        gInd = 1;
                      end
                      ind = i + (j-1)*numData;
                      beta_ij = 1/yvar(ind);
                      dWdf(p,q) = dWdf(p,q)+beta_ij*model.g_grad2(q,gInd)*model.g_grad(p,gInd)*exp(-model.D(j)*(arg1+arg2) ...
                            +log(S2(j))+log(step2));
                    end
                end
            end
        end
    end
    
    % part2: p=q
    for i=1:numData
        arg1 = model.t(i)-model.mapt(p);
        arg2 = model.t(i)-model.mapt(l);
        if arg1 >= 0 && arg2 >= 0
            for j=1:model.numGenes
              if model.ngParam
                gInd = j;
              else
                gInd = 1;
              end
              ind = i + (j-1)*numData;
              beta_ij=1/yvar(ind);
              dWdf(p,p)=dWdf(p,p)+beta_ij*model.g_grad(l,gInd)*model.g_grad2(p,gInd)*exp(-model.D(j)*(arg1+arg2)  ...
                    +log(S2(j)) +log(step2));
            end
        end
    end

end

% part3: p=q=l
temp = 0;
for i=1:numData
    arg = model.t(i)-model.mapt(l);
    if arg >= 0  && l > 1 %% need re-visit
        for j=1:model.numGenes
          if model.ngParam
            gInd = j;
          else
            gInd = 1;
          end
          ind = i + (j-1)*numData;
          beta_ij=1/yvar(ind);
          factor=(model.ypred(model.times_index(i), j)-model.y(ind))*beta_ij;
          temp=temp+factor*model.g_grad3(l,gInd)*exp(-model.D(j)*arg + ...
                log(model.S(j)) +log(model.step));
        end
    end
end
dWdf(l, l) =  dWdf(l, l) + temp;

