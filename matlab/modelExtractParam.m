function [params, names] = modelExtractParam(model, dim)

% MODELEXTRACTPARAM Extract the parameters of a model.
% FORMAT
% DESC Extract parameters from the model into a vector of
% parameters for optimisation.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model.
%
% SEEALSO : modelExpandParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% MODIFICATIONS : Cark Henrik Ek, 2007
%
% MODIFICATIONS : Mauricio Alvarez, 2008

% MLTOOLS

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
