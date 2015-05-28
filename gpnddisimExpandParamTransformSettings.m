function model = gpnddisimExpandParamTransformSettings(model, paramtransformsettings)

% GPNDDISIMEXPANDPARAM Expand the given parameters' transform settings into a GPNDDISIM structure.
% FORMAT
% DESC takes the given vector of parameters and places them in the
% model structure, it then updates any stored representations that
% are dependent on those parameters, for example kernel matrices
% etc..
% ARG model : the model structure for which parameter transform
% settings are to be updated.
% ARG paramtransformsettings : a cell array (nParams x 1) of
% transform settings for each parameter, to be placed in the model
% structure.
% RETURN model : a returned model structure containing the updated
% parameter transform settings.
% 
% SEEALSO : gpsimCreate, gpsimExtractParam, modelExtractParam, gpsimUpdateKernels
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Jaakko Peltonen, 2011

% GPSIM

% fprintf(1,'gpnddisimExpandParamTransformSettings step1\n');


if length(paramtransformsettings) ~= model.numParams
  error(sprintf('Parameter transform settings cell array is incorrect length %d, expected %d',...
                length(paramtransformsettings),model.numParams));
end
startVal = 1;
endVal = model.kern.nParams;
model.kern = kernExpandParamTransformSettings(model.kern, paramtransformsettings(startVal:endVal));

% fprintf(1,'gpnddisimExpandParamTransformSettings step2\n');


if model.numGenes>0,
  for k=1:model.numGenes,
    model.bTransform(k).transformsettings = paramtransformsettings{endVal+k};
  end
  if isfield(model, 'bprior'),
    model.bprior = priorSetBounds(model.bprior, ...
				  paramtransformsettings(endVal+1:endVal+model.numGenes));
  end
  endVal=endVal+model.numGenes;
  
  if model.use_disimstartmean==1,
    for k=1:model.numGenes,
      model.disimStartMeanTransform(k).transformsettings = ...
	  paramtransformsettings{endVal+k};
    end
    if isfield(model, 'disimStartMeanPrior'),
      model.disimStartMeanPrior = priorSetBounds(model.disimStartMeanPrior, ...
						 paramtransformsettings(endVal+1:endVal+model.numGenes));
    end
    endVal=endVal+model.numGenes;  
  end;  
end;

% fprintf(1,'gpnddisimExpandParamTransformSettings step3\n');


if isfield(model,'disimdecayindices') && ~isempty(model.disimdecayindices),
  % find and store transformation settings related to DISIM
  % decay. Assumes model.disimdecayindices has already been created.
  for disimindex=1:model.numGenes,
    %fprintf('taking disimdecay transformsettings from index %d\n',model.disimdecayindices(disimindex));
    model.disimdecaytransformation(disimindex).transformsettings=...
	paramtransformsettings{model.disimdecayindices(disimindex)};
  end;
end;

if isfield(model,'disimvarianceindices') && ~isempty(model.disimvarianceindices),
  % find and store transformation settings related to DISIM
  % variance. Assumes model.disimvarianceindices has already been created.
  disimvariancetransformationsettings=cell(model.numGenes,1);
  for disimindex=1:model.numGenes,
    model.disimvariancetransformation(disimindex).transformsettings=...
	paramtransformsettings{model.disimvarianceindices(disimindex)};
  end;
end;

if isfield(model,'disimdelayindices') && ~isempty(model.disimdelayindices),
  % find and store transformation settings related to DISIM
  % delay. Assumes model.disimdelayindices has already been created.
  for disimindex=1:model.numGenes,
    model.disimdelaytransformation(disimindex).transformsettings=...
	paramtransformsettings{model.disimdelayindices(disimindex)};
  end;
end;

model.simMeanTransform.transformsettings = paramtransformsettings{endVal+1};
if isfield(model, 'simMeanPrior'),
  model.simMeanPrior = priorSetBounds(model.simMeanPrior, ...
				      paramtransformsettings{endVal+1});
end

% fprintf(1,'gpnddisimExpandParamTransformSettings done\n');
