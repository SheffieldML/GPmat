function model = gpsimMapCreate(numGenes, numProteins, times, geneVals, ...
                             geneVars, options)

% GPSIMMAPCREATE Create a GPSIMMAP model.
% The GPSIMMAP model is a model for estimating the protein
% concentration in a small gene network where several genes are
% governed by one protein. The model is based on Gaussian processes
% and simple linear differential equations of the form
%
% dx(t)/dt = B + Cf(t) - Dx(t)
%
% where x(t) is a given genes concentration and f(t) is the protein
% concentration. 
%
% FORMAT
% DESC creates a model for single input motifs with Gaussian
% processes.
% ARG numGenes : number of genes to be modelled in the system.
% ARG numProteins : number of proteins to be modelled in the
% system.
% ARG times : the time points where the data is to be modelled.
% ARG geneVals : the values of each gene at the different time points.
% ARG geneVars : the varuabces of each gene at the different time points.
% ARG options : options structure, the default options can be
% generated using gpsimMapOptions.
% RETURN model : model structure containing default
% parameterisation.
%
% SEEALSO : modelCreate, gpsimMapOptions
%
% COPYRIGHT : Neil D. Lawrence, 2006

% GPSIM

if any(size(geneVars)~=size(geneVals))
  error('The gene variances have a different size matrix to the gene values.');
end

if(numGenes ~= size(geneVals, 2))
  error('The number of genes given does not match the dimension of the gene values given.')
end

if(size(times, 1) ~= size(geneVals, 1))
  error('The number of time points given does not match the number of gene values given')
end

% Initial parameters for the differential equation.
if length(options.B)~=numGenes
  error('Incorrect length of options.B');
end
if length(options.D)~=numGenes
  error('Incorrect length of options.D');
end
if length(options.S)~=numGenes
  error('Incorrect length of options.S');
end


%data collection times
model.yvar = geneVars;
model.y = geneVals;
model.t = times;
model.type = 'gpsimMap';

% Assign options
model.nonLinearity = options.nonLinearity;
model.B = options.B;
model.D = options.D;
model.S = options.S;
model.optimiser = options.optimiser;

%spacing used for path integration.
startPoint = min(times);
endPoint = max(times);
span = endPoint - startPoint;
startPoint = startPoint - span/6;
endPoint = endPoint + span/6;
fullSpan = endPoint - startPoint;
if isfield(options, 'intPoints')
  model.step = fullSpan/(options.intPoints-1);
else 
  model.step = options.step;
end
model.mapt=[];
model.mapt=[model.mapt,startPoint:model.step:times(1)];
model.times_index = [length(model.mapt)];
for i=1:length(times)-1
  model.mapt=[model.mapt,times(i)+model.step:model.step:times(i+1)];
  model.times_index = [model.times_index,length(model.mapt)];
end
model.mapt=[model.mapt,times(length(times))+model.step:model.step:endPoint]'; %time vector

model.numMapPts = length(model.mapt);
model.numGenes = numGenes;


% Initialise the kernel.
if isstruct(options.kern) 
  model.kern = options.kern;
else
  model.kern = kernCreate(model.mapt, options.kern);
end
model.numParams = model.kern.nParams;

% Initialise posterior at zero
model.f = zeros(length(model.mapt), 1);
model.varf = ones(size(model.f));   %Variances

% Initialise data portion of Hessian to zero.
model.W = zeros(length(model.mapt));
% Update the data Hessian.
model.updateW = true;

% Update some internal representations.
model = gpsimMapUpdateG(model);
model = gpsimMapUpdateYpred(model);

% Switch off update of f to update kernels & W.
model.updateF = false;
model = gpsimMapUpdateKernels(model);
model = gpsimMapFunctionalUpdateW(model);
model.updateF = true;