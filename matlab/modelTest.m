function modelRet = modelTest(modelType, numIn, varargin)

% MODELTEST Run some tests on the specified model.
% FORMAT
% DESC runs some tests on the specified model to ensure it is
% correctly implemented.
% ARG modelType : type of model to test. For example, 'linear' or
% 'mlp'.
% RETURN model : the model that was generated for the tests.
% 
% SEEALSO : modelCreate
%
% COPYRIGHT : Neil D. Lawrence,  2006

% MLTOOLS

%if exist([modelType 'Test']) == 2
%  feval([modelType 'Test'])
%end
if ~isstruct(modelType)
  if nargin < 2
    numIn = 4;
  end
  numData = 20;
  numOut = 3;
  
  % Generate some x positions.
  x = randn(numData, numIn);
  x2 = randn(numData, numIn);
  options = modelOptions(modelType);
  model = modelCreate(modelType, numIn, numOut, options);
  model = modelParamInit(model);
  model.X = x;
  model.y = modelOut(model, x2);
  
  % Set the parameters randomly.
  params = modelExtractParam(model);
  params = randn(size(params))./sqrt(randn(size(params)).^2);
  model = modelExpandParam(model, params);
else
  model = modelType;
end

% Check model parameter gradient.
modelGradientCheck(model, varargin{:});

if exist([model.type 'OutputGrad'])==2;
  epsilon = 1e-6;
  [params, names] = modelExtractParam(model);
  if length(names) == 0
    for i = 1:model.numParams
      names{i} = ['Param ' num2str(i)];
    end
  end
  origParams = params;
  Lplus = zeros(size(model.X, 1), model.numParams, model.outputDim);
  Lminus = zeros(size(model.X, 1), model.numParams, model.outputDim);
  for i = 1:length(params);
    params = origParams;
    params(i) = origParams(i) + epsilon;
    model = modelExpandParam(model, params);
    Lplus(:, i, :) = reshape(modelOut(model, model.X), ...
                             size(model.X, 1), 1, model.outputDim);
    params(i) = origParams(i) - epsilon;
    model = modelExpandParam(model, params);
    Lminus(:, i, :) = reshape(modelOut(model, model.X), ...
                              size(model.X, 1), 1, model.outputDim);
  end
  params = origParams;
  model = modelExpandParam(model, params);
  gLDiff = .5*(Lplus - Lminus)/epsilon;
  g = modelOutputGrad(model, model.X);
  
  
  outParamMaxDiff = max(max(max(abs(gLDiff-g))));
  if outParamMaxDiff > 2*epsilon
    l = 0;
    for i = 1:length(names)
      if l < length(names{i})
        l = length(names{i});
      end
    end
    
    fprintf([char(repmat(32, 1, l)) '\tanalytic   diffs     delta\n']);
    for i = 1:length(names)
      spaceLen = l - length(names{i});
      space = char(repmat(32, 1, spaceLen));
      fprintf([space names{i} ':\t%4.6f\t%4.6f\t%4.6f\n'], ...
              g(i), gLDiff(i), gLDiff(i) - g(i));
    end
  end
  fprintf('Output param max diff: %2.6f.\n', outParamMaxDiff);
else
  fprintf('No grad of output with respect to param implemented.\n');
end


if exist([model.type 'OutputGradX'])==2;
  epsilon = 1e-6;
  X = model.X;
  origX = X;
  Lplus = zeros(size(model.X, 1), size(model.X, 2), model.outputDim);
  Lminus = zeros(size(model.X, 1), size(model.X, 2), model.outputDim);
  for i = 1:size(model.X, 2)
    X = origX;
    X(:, i) = origX(:, i) + epsilon;
    Lplus(:, i, :) = reshape(modelOut(model, X), ...
                             size(model.X, 1), 1, model.outputDim);
    X(:, i) = origX(:, i) - epsilon;
    Lminus(:, i, :) = reshape(modelOut(model, X), ...
                              size(model.X, 1), 1, model.outputDim);
  end
  X = origX;
  gLDiff = .5*(Lplus - Lminus)/epsilon;
  g = modelOutputGradX(model, X);
  
  
  outputXMaxDiff = max(max(max(abs(gLDiff-g))));
  if outputXMaxDiff > 2*epsilon
    fprintf('gX\n')
    disp(g)
    fprintf('gXDiff\n')
    disp(gLDiff)
  end
  fprintf('X max diff: %2.6f.\n', outputXMaxDiff);
else
  fprintf('No grad of output with respect to X implemented.\n');
end
if nargout > 0
  modelRet = model;
else
  modelDisplay(model);
end

