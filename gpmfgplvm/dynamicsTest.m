function dynamics = dynamicsTest(type);

% DYNAMICSTEST Run some tests on the specified dynamics model.

% FGPLVM

type = [type 'Dynamics'];
X = randn(10, 2);
model = modelCreate(type, 2, 2, X);

% Set the parameters randomly.
params = modelExtractParam(model);
params = randn(size(params))./sqrt(randn(size(params)).^2);
model = modelExpandParam(model, params);

epsilon = 1e-6;
params = modelExtractParam(model);
if ~isempty(params)
  origParams = params;
  for i = 1:length(params);
    params = origParams;
    params(i) = origParams(i) + epsilon;
    model = modelExpandParam(model, params);
    Lplus(i) = modelLogLikelihood(model);
    params(i) = origParams(i) - epsilon;
    model = modelExpandParam(model, params);
    Lminus(i) = modelLogLikelihood(model);
  end
  params = origParams;
  model = modelExpandParam(model, params);
  
  gLDiff = .5*(Lplus - Lminus)/epsilon;
  g = modelGradientParam(model);
  nameBase = 'Parameter ';
  for i = 1:length(params)
    names{i} = [nameBase num2str(i)];
  end
  paramMaxDiff = max(max(abs(gLDiff-g)));
  if paramMaxDiff > 2*epsilon
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
end
Lplus = zeros(size(X));
Lminus = zeros(size(X));
origX = X;
epsilon = 0.001;
for i = 1:size(X, 1)
  for j = 1:size(X, 2)
    X = origX;
    X(i, j) = origX(i, j) + epsilon;
    model = modelSetLatentValues(model, X);
    params = modelExtractParam(model);
    model = modelExpandParam(model, params);
    Lplus(i, j) = modelLogLikelihood(model);
    X(i, j) = origX(i, j) - epsilon;
    model = modelSetLatentValues(model, X);
    params = modelExtractParam(model);
    model = modelExpandParam(model, params);
    Lminus(i, j) = modelLogLikelihood(model);
  end
end
X = origX;
model = modelSetLatentValues(model, X);
gXDiff = .5*(Lplus - Lminus)/epsilon;
g = modelLatentGradients(model);
XMaxDiff = max(max(abs(g-gXDiff)));
fprintf('X max diff: %2.6f.\n', XMaxDiff)
if XMaxDiff > 1e-6
  disp(g)
  disp(gXDiff)
  disp(g - gXDiff);
end
