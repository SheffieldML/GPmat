function model = gpsimAddCandidate(model, numCandidates, times, geneVals, ...
                                   geneVars, options)

% GPSIMADDCANDIDATE Add candidate genes to a GPSIM model.
% FORMAT
% DESC adds candidate genes to the GPSIM model for optimisation and log
% likelihood evaluation.
% ARG model : the model to add dynamics to.
% ARG numCandidates : number of candiate genes to add in.
% ARG times : times of the candidate observations.
% ARG geneVals : the values of the target gene expressions.
% ARG geneVars : the variances of the target gene expressions.
% ARG options : optins structure as obtained from gpsimCandidateOptions
% RETURN model : model with the candidate field added.
%
% SEEALSO : gpsimCreate, gpsimCandidateOptions, gpsimCandidateOptimise, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2007

% SHEFFIELDML

model.candidate.y = geneVals(:);
model.candidate.t = times(:);
model.candidate.yvar = geneVars(:);
kern = model.kern;
kern.numBlocks = model.kern.numBlocks + numCandidates;


for i = model.kern.numBlocks + 1:kern.numBlocks
  kern.comp{i} = 'sim';
end

numNewParams = 0;
counter = model.kern.numBlocks;
for i = 1:numCandidates
  counter = counter + 1;
  if ~kern.comp{i}.isStationary
    kern.isStationary = false;
  end
  kern.comp{counter} = kernCreate(model.candidate.t, kern.comp{counter});
  kern.comp{counter} = kernParamInit(kern.comp{counter});
  numNewParams = numNewParams + kern.comp{counter}.nParams;
  kern.comp{counter}.index = [];
  for j = 1:counter-1
    func = [kern.comp{counter}.type 'X' kern.comp{j}.type];
    if exist([func 'KernCompute']) == 2
      kern.block{counter}.cross{j} = func;
      kern.block{counter}.transpose(j) = false;
    else
      func = [kern.comp{j}.type 'X' kern.comp{counter}.type];
      if exist([func 'KernCompute']) == 2
        kern.block{counter}.cross{j} = func;
        kern.block{counter}.transpose(j) = true;
      else
        warning(['No cross covariance found between ' kern.comp{counter}.type ...
                 ' and ' kern.comp{j}.type ' assuming independence.'])
        kern.block{counter}.cross{j} = [];
        kern.block{counter}.transpose(j) = 0;
      end
    end
  end
end

kern.nParams = kern.nParams + numNewParams;

%
kern.paramGroups = [kern.paramGroups zeros(size(kern.paramGroups, 1), numNewParams);  
                    zeros(numNewParams, size(kern.paramGroups, 2)) speye(numNewParams)];

% Tie RBF kernel parameters.
tieParam = [2 model.kern.nParams+2:3:kern.nParams];
kern = modelTieParam(kern, {tieParam});
param = kernExtractParam(kern);
kern = kernExpandParam(kern, param);

%kern = modelTieParam(kern, 
model.candidate.kern = kern;

counter = 0;
for i = model.kern.numBlocks+1:model.candidate.kern.numBlocks
  counter = counter + 1;
  model.candidate.D(counter) = model.candidate.kern.comp{i}.decay;
  model.candidate.S(counter) = sqrt(model.candidate.kern.comp{i}.variance);
end

model.candidate.numParams = numCandidates + model.candidate.kern.nParams-model.kern.nParams;
model.candidate.numGenes = numCandidates;
model.candidate.mu = mean(geneVals);
model.candidate.B = model.candidate.D.*model.candidate.mu;
model.candidate.m = model.candidate.y;
model.candidate.t = times;

model.candidate.optimiser = options.optimiser;

if isfield(options, 'fix')
  model.candidate.fix = options.fix;
end

% The basal transcriptions rates must be postitive.
model.candidate.bTransform = optimiDefaultConstraint('positive');

% This forces kernel compute.
params = gpsimCandidateExtractParam(model);
model = gpsimCandidateExpandParam(model, params);


