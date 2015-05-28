function g = gpdisimLogLikeGradients(model)

% GPDISIMLOGLIKEGRADIENTS Compute the gradients of the log likelihood of a GPDISIM model.
% FORMAT
% DESC computes the gradients of the log likelihood of the given
% Gaussian process for use in a single input motif protein network.
% ARG model : the model for which the log likelihood is computed.
% RETURN g : the gradients of the parameters of the model.
% 
% SEEALSO : gpsimCreate, gpsimLogLikelihood, gpsimGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007
  
% SHEFFIELDML

covGrad = -model.invK + model.invK*model.m*model.m'*model.invK;
covGrad = 0.5*covGrad;
if isfield(model, 'proteinPrior') && ~isempty(model.proteinPrior)
  g = kernGradient(model.kern, model.timesCell, covGrad);
else
  g = kernGradient(model.kern, model.t, covGrad);
end

% In case we need priors in.
% Add contribution of any priors 
if isfield(model, 'bprior'),
  g = g + kernPriorGradient(model.kern);
end


gmuFull = model.m'*model.invK;

% if isfield(model, 'proteinPrior') && ~isempty(model.proteinPrior)
%   if model.includeNoise
%     ind = model.kern.comp{1}.diagBlockDim{1} + (1:model.kern.comp{1}.diagBlockDim{2});
%     gmu = zeros(size(1, model.numGenes));

%     for i = 1:model.numGenes
%       gmu(i) = sum(gmuFull(ind));
%       ind = ind + model.kern.comp{1}.diagBlockDim{i+1};
%     end
%   else
%     ind = model.kern.diagBlockDim{1} + (1:model.kern.diagBlockDim{2});
%     gmu = zeros(size(1, model.numGenes));

%     for i = 1:model.numGenes
%       gmu(i) = sum(gmuFull(ind));
%       ind = ind + model.kern.diagBlockDim{i+1};
%     end
%   end
  
% else
  numData = size(model.t, 1);
  ind = 1:numData;
  ind = ind + numData;
  gmu = zeros(size(1, model.numGenes));
  for i = 1:model.numGenes
    gmu(i) = sum(gmuFull(ind));
    ind = ind + numData;
  end
%end

gb = gmu./model.D;
fhandle = str2func([model.bTransform 'Transform']);
% In case we need priors in.
% Add prior on B if it exists.
if isfield(model, 'bprior');
  gb = gb + priorGradient(model.bprior, model.B);
end
 
gb = gb.*fhandle(model.B, 'gradfact');

% This is a nasty hack to add the influence of the D in the mean to
% the gradient already computed for the kernel. This is all very
% clunky and sensitive to changes that take place elsewhere in the
% code ...
gd = -gmu.*model.B./(model.D.*model.D);
decayIndices = [5];
for i = 3:model.kern.numBlocks
  decayIndices(end+1) = decayIndices(end) + 2;
end    

% Account for decay in mean.
g(decayIndices) = g(decayIndices) ...
    + gd.*expTransform(model.D, 'gradfact');

g = [g gb];

if isfield(model, 'fix')
  for i = 1:length(model.fix)
    g(model.fix(i).index) = 0;
  end
end
