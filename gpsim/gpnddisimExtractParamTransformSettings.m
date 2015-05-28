function [paramtransformsettings, names] = gpnddisimExtractParamTransformSettings(model)

% NDDISIMEXTRACTPARAMTRANSFORMSETTINGS Extract the parameter transform settings of a GPNDDISIM model.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a Gaussian process for single input motif
% modelling.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
%
% SEEALSO : gpdisimCreate, gpdisimExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007
%
% COPYRIGHT : Jaakko Peltonen, 2011

% GPSIM

if nargout>1
  [paramtransformsettings, names] = kernExtractParamTransformSettings(model.kern);

  for i=1:model.numGenes,
    names{end+1}=['Basal transcription ' num2str(i)];
  end

  if model.numGenes>0,
    if (model.use_disimstartmean==1),
      for i=1:model.numGenes,
	names{end+1}=['NDDISIM startmean ' num2str(i)];
      end;
    end;
  end;
  
  names{end+1}=['NDSIM mean ' num2str(i)];
else
  paramtransformsettings = kernExtractParamTransformSettings(model.kern);
end


% Transformation settings for basal transcription rates
if model.numGenes>0,
  for k=1:model.numGenes,
    paramtransformsettings = {paramtransformsettings{:}, model.bTransform(k).transformsettings};
  end;
end;

% Transformation settings for DISIM start means
if model.numGenes>0,
  if (model.use_disimstartmean==1),
    for i=1:model.numGenes,
      paramtransformsettings = {paramtransformsettings{:}, model.disimStartMeanTransform(k).transformsettings};
    end;
  end;
end;

% Transformation settings for SIM mean
paramtransformsettings = {paramtransformsettings{:}, model.simMeanTransform.transformsettings};



