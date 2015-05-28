function prior = priorTest(priorType);

% PRIORTEST Run some tests on the specified prior.
% FORMAT
% DESC runs some checks on the specified prior (gradients etc.).
% ARG priorType : the prior type to check.
% RETURN prior : the prior structure that was created to test.
%
% SEEALSO : priorParamInit, priorExpandParam, priorExtractParam,
% priorLogProb, priorGradient
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004
%
% MODIFICATIONS : Antti Honkela, 2012, 2014

% PRIOR

prior.type = priorType;
prior = priorParamInit(prior);

% Set the parameters randomly.
params = priorExtractParam(prior);
params = randn(size(params))./sqrt(randn(size(params)).^2);
prior = priorExpandParam(prior, params);
if isfield(prior, 'isBounded') && prior.isBounded,
  prior.b = 1+9*rand(1);
  x = rand(1, 10) * (prior.b - prior.a) + prior.a;
else
  x = randn(1, 10);
end
epsilon = 1e-6;

if exist([priorType 'PriorGradientParams'])
  params = priorExtractParam(prior);
  origParams = params;
  for i = 1:length(params);
    params = origParams;
    params(i) = origParams(i) + epsilon;
    prior = priorExpandParam(prior, params);
    Lplus(i) = priorLogProb(prior, x);
    params(i) = origParams(i) - epsilon;
    prior = priorExpandParam(prior, params);
    Lminus(i) = priorLogProb(prior, x);
  end
  params = origParams;
  prior = priorExpandParam(prior, params);
  [void, names] = priorExtractParam(prior);
  gLDiff = .5*(Lplus - Lminus)/epsilon;
  g = priorGradient(prior, x);
  
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
if exist([priorType 'PriorGradient'])
  origX = x;
  for i = 1:length(x);
    x = origX;
    x(i) = origX(i) + epsilon;
    Lplus(i) = priorLogProb(prior, x);
    x(i) = origX(i) - epsilon;
    Lminus(i) = priorLogProb(prior, x);
  end
  x = origX;
  gLDiff = .5*(Lplus - Lminus)/epsilon;
  g = priorGradient(prior, x);
  
  paramMaxDiff = max(max(abs(gLDiff-g)));
%  if paramMaxDiff > 2*epsilon
    
    fprintf([ '\tanalytic   diffs     delta\n']);
    for i = 1:length(x)
      fprintf([num2str(i) ':\t%4.6f\t%4.6f\t%4.6f\n'], ...
              g(i), gLDiff(i), gLDiff(i) - g(i));
    end
%  end
end
