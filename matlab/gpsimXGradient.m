function [dxdB dxdD dxdS dxdalpha dxdgParam] = gpsimXGradient(model, i, j)

% GPSIMXGRADIENT ...
%
% COPYRIGHT : Pei Gao, 2008
  
% SHEFFIELDML

% i: ti
% j: jth gene

dxdB = 1./model.D(j);

endPoint = model.times_index(i);

dxdD = 0;
dxdS = 0;
dxdalpha = [];
dxdgParam = [];

if model.ngParam > 0
  ngParamk = model.ngParam/model.numGenes;
  dxdgParam = zeros(1, ngParamk);
  gInd = j;
else
  gInd = 1;
end

for m = 1:endPoint  
    arg = model.t(i)-model.mapt(m);
    if arg >= 0
        dxdD = dxdD+model.g(m,gInd)*arg*exp(-model.D(j)*arg+log(model.step)+log(model.S(j)));
        dxdS = dxdS+exp(-model.D(j)*arg+log(model.step))*model.g(m,gInd); % g=Sf
      
        if model.ngParam > 0
          for gParamInd = 1:ngParamk
            dxdgParam(gParamInd) = dxdgParam(gParamInd) + exp(- ...
            model.D(j)*arg+log(model.step)+log(model.S(j)))*model.dg(m,gInd);
          end
        end
    end
end

dxdD = -model.B(j)/(model.D(j)*model.D(j))-dxdD;

% check if g(f) is different for each gene.
if isfield(model, 'isGroupNonlinearity') && strcmp(model.nonLinearity{j}, ...
                                                   'repression')
      dxdalpha = exp(-model.D(j)*model.t(i));
      dxdD = -model.t(i)*model.alpha(j)*exp(-model.D(j)*model.t(i))+dxdD;
end
          
