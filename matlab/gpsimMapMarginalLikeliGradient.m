function g = gpsimMapMarginalLikeliGradient(model)
  
% GPSIMMAPMARGINALLIKELIGRADIENT Compute the gradients of the log marginal likelihood of a GPSIMMAP model with respect to the model parameters.

% FORMAT
% DESC computes the gradients of the log marginal likelihood of the given
% Gaussian process for use in a single input motif protein network.
% ARG model : the model for which the log likelihood is computed.
% RETURN g : the gradients of the parameters of the model.
% 
% SEEALSO : gpsimMapCreate, gpsimMapLogLikelihood,
% gpsimMapGradient, gpsimMapFunctionalLogLikeGradients
%
% COPYRIGHT : Pei Gao, Magnus Rattray and Neil D. Lawrence, 2008

% SHEFFIELDML

numData = length(model.t);
ParamL = length(model.B);

if isfield(model, 'includeNoise') && model.includeNoise
  noiseMat = ones(numData, 1)*model.noiseVar;
  yvar = model.yvar + noiseMat;
  dlogPdn = zeros(ParamL,1);
else
  yvar = model.yvar;
end

if isfield(model, 'includeRepression') && model.includeRepression
    dlogPdalpha = zeros(ParamL,1);
end

dlogPdB = zeros(ParamL,1);
dlogPdD = zeros(ParamL,1);
dlogPdS = zeros(ParamL,1);

if model.ngParam > 0
  ngParamk = model.ngParam/model.numGenes;
  dlogPdgParam= zeros(ParamL, ngParamk);
end

% C = eye(size(model.W)) + model.K*model.W;
% invC = inv(C)*model.K;
invC = model.covf;

for k = 1:model.numGenes
    for i=1:numData
        ind = i + (k-1)*numData;
        beta_ik=1/yvar(ind);
        factor = (model.ypred(model.times_index(i), k)-model.y(ind))*beta_ik;

        [dxdB dxdD dxdS dxdalpha dxdgParam] = gpsimXGradient(model, i, k);

        dlogPdB(k)=dlogPdB(k)-factor*dxdB;

        if isfield(model, 'includeRepression') && model.includeRepression
          dlogPdalpha(k) = dlogPdalpha(k)-factor*dxdalpha;
        end             

        dlogPdD(k) = dlogPdD(k)-factor*dxdD;

        dlogPdS(k) = dlogPdS(k)-factor*dxdS;

        if model.ngParam > 0
            dlogPdgParam(k,:)= dlogPdgParam(k,:)-factor*dxdgParam;
        end

        if isfield(model,'includeNoise') && model.includeNoise
          dlogPdn(k) = dlogPdn(k) + sqrt(model.noiseVar(k))*factor*factor - ...
              sqrt(model.noiseVar(k))*beta_ik;
        end
    end

%    gB(k) = dlogPdB(k);
%    gD(k) = dlogPdD(k);
%    gS(k) = dlogPdS(k);
%    ggParam(:,k) = dlogPdgParam(k,:)';
%    gn(k) = dlogPdn(k);

    [dWdB, dWdD, dWdS, dWdalpha, dWdgParam, dWdn] = gpsimMapWGradient(model, k);

    gB(k) = -0.5*trace(invC*dWdB)+dlogPdB(k);
    gD(k) = -0.5*trace(invC*dWdD)+dlogPdD(k);
    gS(k) = -0.5*trace(invC*dWdS)+dlogPdS(k);
    
    if isfield(model, 'includeRepression') && model.includeRepression
      galpha(k) = -0.5*trace(invC*dWdalpha)+dlogPdalpha(k);
    end      

    if model.ngParam > 0
      for gParam_index = 1:ngParamk
        ggParam(gParam_index,k) = -0.5*trace(invC*dWdgParam(:,:, ...
                         gParam_index))+dlogPdgParam(k, gParam_index);
      end
    end
    
    if isfield(model, 'includeNoise') && model.includeNoise
      gn(k) = -0.5*trace(invC*dWdn(:,:)) + dlogPdn(k);
    end

end

fhandle = str2func([model.Transform 'Transform']);
g=[];

for i = 1:model.numGenes
  if isfield(model,'bTransform') && isempty(model.bTransform)
    g = [g gB(i)];
  else
    g = [g gB(i)*fhandle(model.B(i), 'gradfact')];
  end
  
  g = [g gS(i)*fhandle(model.S(i), 'gradfact') gD(i)*fhandle(model.D(i), 'gradfact')];
    
  if isfield(model, 'includeRepression') && model.includeRepression
    if isfield(model,'alphaTransform') && isempty(model.alphaTransform)
      g = [g galpha(i)];
    else
      g = [g galpha(i)*fhandle(model.alpha(i), 'gradfact')];      
    end    
  end      
   
  if model.ngParam > 0
    g = [g ggParam(:,i).*fhandle(model.gParam(:,i), 'gradfact')];
  end
end

if isfield(model, 'includeNoise') && model.includeNoise
  g = [g gn.*fhandle(sqrt(model.noiseVar), 'gradfact')];
end
