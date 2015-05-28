function model = lfmCreate(numDisplacements, numForces, times, displacementVals,  options)

% LFMCREATE Create a LFM model.
% The LFM model is a model for estimating the force
% driving a small system where several displacements are
% governed by one force. The model is based on Gaussian processes
% and simple linear differential equations of the form
%
% m d2x(t)/dt2 + Cdx(t)/dt + Dx(t)= Sf(t) 
%
% where x(t) is a given displacement and f(t) is the governing force.
%
% FORMAT
% DESC creates a model for single input motifs with Gaussian
% processes.
% ARG numDisplacements : number of displacements to be modelled in the system.
% ARG numForces : number of forces to be modelled in the
% system.
% ARG times : the time points where the data is to be modelled.
% ARG displacementVals : the values of each displacement at the different time points.
% ARG options : options structure, the default options can be
% displacementrated using lfmOptions.
% RETURN model : model structure containing default
% parameterisation.
%
% SEEALSO : modelCreate, lfmOptions
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN

if(numDisplacements ~= size(displacementVals, 2))
  error('The number of displacements given does not match the dimension of the displacement values given.')
end

if(size(times, 1) ~= size(displacementVals, 1))
  error('The number of time points given does not match the number of displacement values given')
end

model.type = 'lfm';

kernType1{1} = 'multi';
kernType2{1} = 'multi';
tieParam = [4]; % These are the indices of the inverse widths which
                % need to be constrained to be equal.
for i = 1:numDisplacements
  kernType1{i+1} = 'lfm';
  if i>1
    tieParam = [tieParam tieParam(end)+5];
  end
end

model.y = displacementVals(:);
model.kern = kernCreate(times, kernType1);
%/~ This is if we need to place priors on parameters ...
% for i = 1:length(model.kern.numBlocks)
%   % Priors on the lfm kernels.
%   model.kern.comp{i}.priors = priorCreate('gamma');
%   model.kern.comp{i}.priors.a = 1;
%   model.kern.comp{i}.priors.b = 1;
%   if i == 1
%     % For first kernel place prior on inverse width.
%     model.kern.comp{i}.priors.index = [1 2 3];
%   else
%     % For other kernels don't place prior on inverse width --- as
%     % they are all tied together and it will be counted multiple
%     % times.
%     model.kern.comp{i}.priors.index = [1 3];
%   end
% end

% Prior on the b values.
% model.bprior = priorCreate('gamma');
% model.bprior.a = 1;
% model.bprior.b = 1;
%~/
model.kern = modelTieParam(model.kern, {tieParam});

% The decays and sensitivities are actually stored in the kernel.
% We'll put them here as well for convenience.
for i = 1:model.kern.numBlocks
  model.spring(i) = model.kern.comp{i}.spring;
  model.sensitivity(i) = model.kern.comp{i}.sensitivity;
  model.mass(i) = model.kern.comp{i}.mass;
  model.damper(i) = model.kern.comp{i}.damper;
end

model.numParams = model.kern.nParams;
model.numDisplacements = numDisplacements;
model.m = model.y;
model.t = times;

model.optimiser = options.optimiser;

if isfield(options, 'fix')
  model.fix = options.fix;
end

% This forces kernel compute.
params = lfmExtractParam(model);
model = lfmExpandParam(model, params);

