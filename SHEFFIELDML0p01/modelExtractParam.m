function [params, names] = modelExtractParam(model, dim)

% MODELEXTRACTPARAM Extract the parameters of a model.
%
%	Description:
%
%	PARAM = MODELEXTRACTPARAM(MODEL) Extract parameters from the model
%	into a vector of parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the model.
%	 Arguments:
%	  MODEL - the model structure containing the parameters to be
%	   extracted.
%	
%	
%	
%
%	See also
%	MODELEXPANDPARAM, SCG, CONJGRAD


%	Copyright (c) 2005, 2006 Neil D. Lawrence


%	With modifications by Cark Henrik Ek 2007


%	With modifications by Mauricio Alvarez 2008


fhandle = str2func([model.type 'ExtractParam']);
if nargout < 2
    if(nargin<2)
        params = fhandle(model);
    else
        params = fhandle(model,dim);
    end
else
    [params, namesTemp] = fhandle(model);
end
if isfield(model, 'paramGroups')
    paramGroups = model.paramGroups;
    for i = 1:size(paramGroups, 2)
        ind = find(paramGroups(:, i));
        if nargout > 1
            names{i} = namesTemp{ind(1)};
            if length(ind) > 1
                for j = 2:length(ind)
                    names{i} = [names{i} ', ' namesTemp{ind(j)}];
                end
            end
        end
        paramGroups(ind(2:end), i) = 0;
    end
    params = params*paramGroups;
else
    if nargout>1
        names = namesTemp;
    end
end
end
