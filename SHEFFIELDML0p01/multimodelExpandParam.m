function model = multimodelExpandParam(model, params)

% MULTIMODELEXPANDPARAM Create model structure from MULTIMODEL model's parameters.
%
%	Description:
%
%	MODEL = MULTIMODELEXPANDPARAM(MODEL, PARAM) returns a multi-task
%	learning wrapper model structure filled with the parameters in the
%	given vector. This is used as a helper function to enable parameters
%	to be optimised in, for example, the NETLAB optimisation functions.
%	 Returns:
%	  MODEL - model structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  MODEL - the model structure in which the parameters are to be
%	   placed.
%	  PARAM - vector of parameters which are to be placed in the model
%	   structure.
%	
%	
%
%	See also
%	MULTIMODELCREATE, MULTIMODELEXTRACTPARAM, MODELEXPANDPARAM


%	Copyright (c) 2007, 2008 Neil D. Lawrence


%	With modifications by Mauricio Alvarez 2008, 2009


if numel(model.outputDim) == 1
    endVal = model.numParams - model.numModels*model.numSep;
    baseParams = params(1:endVal);
    passParams(model.sharedIndices) = baseParams;
    if ~isempty(model.separateIndices)
        for i = 1:length(model.comp)
            startVal = endVal + 1;
            endVal = endVal + model.numSep;
            passParams(model.separateIndices) = params(startVal:endVal);
            model.comp{i} = modelExpandParam(model.comp{i}, passParams);
        end
    else
        for i = 1:length(model.comp)
            model.comp{i} = modelExpandParam(model.comp{i}, passParams);
        end
    end
else
    startVal = 1;
    endVal = 0;
    for i = 1:length(model.comp)
        endVal = endVal + model.comp{i}.nParams;        
        model.comp{i} = modelExpandParam(model.comp{i}, params(startVal:endVal));
        startVal = endVal + 1;
    end
end