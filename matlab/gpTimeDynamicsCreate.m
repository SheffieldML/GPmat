function model = gpTimeDynamicsCreate(q, d, latentVals, options, varargin)

% GPTIMEDYNAMICSCREATE Create the time dynamics model. 
% FORMAT
% DESC creates a Gaussian process model for dealing with dynamics
% in the latent space of a GP-LVM. The input to the dynamics model
% is time, rather than the previous observation as is the case for GPDYNAMICS.
% ARG q : the latent space dimension.
% ARG q : the latent space dimension.
% ARG X : the latent variables.
% ARG options : options structure as defined by gpOptions.m.
% ARG t : the input time values.
% ARG diff : Whether or not to use differences between points in
% the latent space as the targets for the GP or absolute location of
% points (default 1).
% ARG learn : Whether or not to learn the parameters of the
% dynamics model (default 0).
% ARG seq : array containing the indices of the last frame of each
% sequence in the data. The array should be the same length as the
% number of sequences in the data (default N, where N is total
% number of frames in the model).
% the latent space as the inputs to the GP or absolute location of points.
% RETURN model : model structure containing the Gaussian process.
%
% SEEALSO : gpCreate, gpDynamicsCreate, gpTimeDynamicsLatentGradients, gpTimeDynamicsSetLatentValues, gpTimeDynamicsLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Carl Henrik Ek, 2008
%

% FGPLVM

if nargin>4
  t = varargin{1};
else
  t = [1:size(latentVals, 1)]';
end
if nargin>5
  diff = varargin{2};
else 
  diff = 0;
end
if nargin > 6
  learn = varargin{3};
else
  learn = 0;
end
if nargin > 7
  seq = varargin{4};
else
  seq = size(latentVals, 1);
end
  
if(iscell(options))
  varargin = options;
  options = varargin{1};
  t = varargin{2};
  diff = varargin{3};
  learn = varargin{4};
  seq = varargin{5};
end

% set source and target
startPoints = [];t2 = []; y = [];
startVal = 1;
for i = 1:length(seq)
  endVal = seq(i);
  startPoints = [startPoints; latentVals(startVal:endVal-1, :)];
  t2 = [t2; t(startVal+1:endVal)];
  y = [y; latentVals(startVal+1:endVal, :)];
  startVal = endVal + 1;
end

if diff
  y = y - startPoints;
end
model = gpCreate(size(t2, 2), d, t2, y, options);

model.diff = diff;
model.learn = learn;
model.type = 'gpTimeDynamics';
model.seq = seq;
model.dynamicsType = 'regressive';
