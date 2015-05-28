function kern = lfmwhiteKernParamInit(kern)

% LFMWHITEKERNPARAMINIT LFM-WHITE kernel parameter initialisation.
% FORMAT
% DESC initialises the LFM-White (Latent Force Model - White) kernel
% structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : David Luengo, 2009
%
% COPYRIGHT : Neil D. Lawrence, 2009

  
% KERN


if kern.inputDimension > 1
  error('LFM-WHITE kernel only valid for one-D input.')
end

kern.nParams = 5;
kern.mass = 1;
kern.spring = 1;
kern.damper = 1;
kern.variance = 1;
kern.sensitivity = 1;

kern.delay = 0;
kern.initVal = 1;

kern.transforms.index = [1 2 3 4];
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = false;

% Serial number used to distinguish LFM kernels
maxSerial = double(intmax('uint64'));
kern.serialNumber = uint64(1+rand(1)*maxSerial);

% Force any precomputation contained in lfmKernExpandParam
params = lfmwhiteKernExtractParam(kern);
kern = lfmwhiteKernExpandParam(kern, params);
kern.positiveTime = true;
