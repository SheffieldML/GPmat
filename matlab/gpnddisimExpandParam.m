function model = gpnddisimExpandParam(model, params, varargin)

% GPDISIMEXPANDPARAM Expand the given parameters into a GPDISIM structure.
% FORMAT
% DESC takes the given vector of parameters and places them in the
% model structure, it then updates any stored representations that
% are dependent on those parameters, for example kernel matrices
% etc..
% ARG model : the model structure for which parameters are to be
% updated.
% ARG params : a vector of parameters for placing in the model
% structure.
% RETURN model : a returned model structure containing the updated
% parameters.
% 
% SEEALSO : gpsimCreate, gpsimExtractParam, modelExtractParam, gpsimUpdateKernels
%
% COPYRIGHT : Neil D. Lawrence, 2006

% GPSIM


if length(varargin)==2,
  update_kernel=varargin{1};
  update_mean=varargin{2};
else
  update_kernel=1;
  update_mean=1;
end;


%fprintf('Params received by nddisimExpandParam:\n');
%params

%params
%pause

params = real(params);
if isfield(model, 'fix')
  for i = 1:length(model.fix)
    params(model.fix(i).index) = model.fix(i).value;
  end
end

%fprintf('Params received by nddisimExpandParam, after fixing values:\n');
%params


if length(params) ~= model.numParams
  error(sprintf('Parameter vector is incorrect length %d, expected %d',...
                length(params),model.numParams));
end
startVal = 1;
endVal = model.kern.nParams;


%fprintf(1,'Calling kernExpandParam with these parameters:\n');
%params(startVal:endVal)

model.kern = kernExpandParam(model.kern, params(startVal:endVal));

if model.numGenes>0,
  for i = 1:model.numGenes,  
    model.B(i) = doTransform(params(endVal+i), 'atox', model.bTransform(i));
  end;
  endVal=endVal+model.numGenes;
  
  if model.use_disimstartmean==1,
    for i = 1:model.numGenes,  
      model.disimStartMean(i) = doTransform(params(endVal+i), 'atox', model.disimStartMeanTransform(i));
    end;    
    endVal=endVal+model.numGenes;    
  end;
end;

model.simMean = doTransform(params(endVal+1), 'atox', model.simMeanTransform);



% The decays and sensitivities and delays are actually stored in the kernel.
% We'll put them here as well for convenience.
if model.includeNoise,
  simMultiKern = model.kern.comp{1};
else
  simMultiKern = model.kern;
end

if model.numGenes>0,
  model.sigma = sqrt(simMultiKern.comp{2}.di_variance);
  for i = 2:simMultiKern.numBlocks
    model.D(i-1) = simMultiKern.comp{i}.decay;
    model.S(i-1) = sqrt(simMultiKern.comp{i}.variance);
    model.delay(i-1) = simMultiKern.comp{i}.delay;
    %model.delay(i-1)
    %pause
  end
end;


if (~isempty(update_kernel)) && (update_kernel==1),
  model = gpnddisimUpdateKernels(model);
end;


if (~isempty(update_mean)) && (update_mean==1),

  model.mu = zeros(size(model.y));
  
  % update POL2 mean
  if iscell(model.t)==0, tempt=model.t; else tempt=model.t{1}; end;
  lengthObs = size(tempt, 1);
  ind = 1:lengthObs;
  model.mu(ind) = model.simMean;
  %model.m(ind) = model.y(ind)-model.mu(ind);
  
  % update RNA means
  if iscell(model.t)==0, tempt=model.t; else tempt=model.t{1}; end;
  timeshift = min(tempt);
  nt=size(tempt,1);
  tempind1=nt+1;
  for k=1:model.numGenes,
    if iscell(model.t)==0, tempt=model.t; else tempt=model.t{k+1}; end;
    nt=size(tempt,1);

    if (model.use_disimstartmean==1),
      delayedt=tempt-model.delay(k);
      I=find(delayedt<0);
      delayedt(I)=0;
      model.mu(tempind1:tempind1+nt-1)=...
	  model.disimStartMean(k)*exp(model.D(k)*(-tempt + timeshift)) ...
          +(model.B(k)/model.D(k))*(1-exp(-model.D(k)*tempt)) ...
          +(model.simMean*model.S(k)/model.D(k))*(1-exp(-model.D(k)*delayedt));
      %model.mu
    else
      delayedt=model.t-model.delay(k);
      I=find(delayedt<0);
      delayedt(I)=0;
      model.mu(tempind1:tempind1+nt-1)=...
	  ((model.B(k)+model.simMean*model.S(k))/model.D(k))*exp(model.D(k)*(-tempt))...
          +((model.B(k)+model.simMean*model.S(k))/model.D(k))*(1-exp(-model.D(k)*delayedt));
    end;

    tempind1=tempind1+nt;
  end;
       
  model.m = model.y - model.mu;
end;
