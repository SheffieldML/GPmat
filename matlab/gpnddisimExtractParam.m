function [param, names] = gpnddisimExtractParam(model)

% NDDISIMEXTRACTPARAM Extract the parameters of a NDDISIM model.
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
  [param, names] = kernExtractParam(model.kern);

%  fprintf(1,'gpdisimExtractParam: parameters from kernels\n');
%  param
%  names

  %  model.mu
  %  length(model.mu)
  %for i=1:length(model.mu),
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
  param = kernExtractParam(model.kern);
end




%length(param)
%model.numGenes
%length(model.B)
if model.numGenes>0,
  for k=1:model.numGenes,
    fhandle = str2func([model.bTransform 'Transform']);
    param = [param fhandle(model.B(k), 'xtoa',model.bTransformSettings{k})];
  end;
end;
%length(param)

if model.numGenes>0,
  if (model.use_disimstartmean==1),
    for i=1:model.numGenes,
      fhandle = str2func([model.disimStartMeanTransform 'Transform']);
      param = [param fhandle(model.disimStartMean(k),'xtoa',model.disimStartMeanTransformSettings{k})];
    end;
  end;
end;


param = [param fhandle(model.simMean,'xtoa',model.simMeanTransformSettings)];



if isfield(model, 'fix')
%model.fix
  for i = 1:length(model.fix)
    param(model.fix(i).index) = model.fix(i).value;
  end
end
param = real(param);
