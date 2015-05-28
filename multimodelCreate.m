function model = multimodelCreate(inputDim, outputDim, varargin)

% MULTIMODELCREATE Create a MULTIMODEL model.
% The MULTIMODEL is a way of performing multi-task learning by sharing
% model parameters across a range of models. The default (simple)
% assumption is that the data is conditionally independent given the
% parameters, i.e. the log likelihood is the sum of the log likelihood of
% the models.
%
% SEEALSO : modelCreate
%
% FORMAT
% DESC creates a multi-task learning wrapper
% model structure given an options structure. 
% ARG inputDim : the input dimension of the model.
% ARG outputDim : the output dimension of the model.
% ARG options : an options structure that determines the form of the model.
% RETURN model : the model structure with the default parameters placed in.
%
% SEEALSO : multimodelOptions, multimodelParamInit, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2007
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009

% MLTOOLS

options = varargin{end};
model.numModels = options.numModels;
model.type = 'multimodel';
model.compType = options.type;
model.inputDim = inputDim;
model.outputDim = outputDim;

if numel(outputDim) == 1
    % Indices of parameters to be trained separately for each model.
    model.separateIndices = options.separate;
    model.numSep = length(model.separateIndices);
    for i = 1:model.numModels
        varargput = cell(1, length(varargin)-1);
        for j = 1:length(varargput)
            varargput{j} = varargin{j}{i};
        end
        % MAURICIO : temporarily changed to allow different options for each model
        if ~iscell(options.compOptions)
            model.comp{i} = modelCreate(model.compType, inputDim, outputDim, ...
                varargput{:}, options.compOptions);
        else
            model.comp{i} = modelCreate(model.compType, inputDim, outputDim, ...
                varargput{:}, options.compOptions{i});
        end
    end
    if isfield(model.comp{i}, 'numParams');
        model.numParams = model.comp{1}.numParams;
    else
        model.numParams = length(modelExtractParam(model.comp{1}));
    end
    model.sharedIndices = 1:model.numParams;
    model.sharedIndices(model.separateIndices) = [];    
    model.numParams = model.numParams + (model.numModels-1)*model.numSep;    
else    
    model.separateIndices = options.separate;
    model.numSep = length(model.separateIndices);  
    for i = 1:model.numModels
        varargput = cell(1, length(varargin)-1);
        for j = 1:length(varargput)
            varargput{j} = varargin{j}{i};
        end
        fprintf('Creating model number: %d\n', i)
        model.comp{i} = modelCreate(model.compType, inputDim, outputDim(i), ...
            varargput{:}, options.compOptions{i}); 
        
    end
    if isfield(model.comp{1}, 'nParams');
        model.numParams = 0;
        for i =1:model.numModels,
            model.numParams = model.numParams + model.comp{i}.nParams;
        end
    else
        model.numParams = length(modelExtractParam(model.comp{1}));
    end
    model.numParams = model.numParams + (model.numModels-1)*model.numSep;    
end

if isfield(options, 'optimiser') && ~isempty(options.optimiser)
    model.optimiser = options.optimiser;    
end

