function  kern = sdrbfKernParamInit(kern)

% SDRBFKERNPARAMINIT SDRBF kernel initialization
% FORMAT
% DESC
% Initializes the switching dynamical latent force model structure with
% some initial parameters for the RBF underlying kernels. The initial 
% parameters are passed through an option in kern.
% RETURN kern : the kernel structure with the default parameters placed in.
% ARG kern : the kernel structure which requires initialisation.
%
% SEEALSO : kernCreate, kernParamInit, rbfKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

% Create the basic structure based on the lfm kernel.
kern = rbfKernParamInit(kern);


if isfield(kern, 'options') && isfield(kern.options, 'isNormalised') && ...
        kern.options.isNormalised,
    kern.isNormalised = 1;
end

if isfield(kern, 'options') && isfield(kern.options, 'nIntervals')
    kern.nIntervals = kern.options.nIntervals;
    if isfield(kern.options, 'nlfPerInt')
        kern.nlfPerInt = kern.options.nlfPerInt;
    else
        kern.nlfPerInt = 1;
    end
    kern.inverseWidth = ones(kern.nlfPerInt, kern.nIntervals); % Total number of inverse widths.
    % An option for the initialization of the switching times
    if isfield(kern.options, 'switchingTimes')
        if kern.nIntervals == length(kern.options.switchingTimes)
            kern.switchingTimes = kern.options.switchingTimes;
        else
            error('The number of intervals does not match the information of the swicthing time points')
        end
    else
        partition = linspace(-0.1,1, kern.nIntervals + 1);
        kern.switchingTimes = partition(1:end-1);
    end
else
    kern.nIntervals = 1;
    kern.switchingTimes = -0.1;
    warning('LFM:Instead:SDLFM', 'Use the LFM kernel instead.')
end

% Number of parameters computed as: inverseWidth (= nlf*nIntervals) and the
% switching points (= nIntervals). The variance, which is a parameter in 
% the original rbf kernel, is not counted as a parameter. It's left here, 
% so that it can be used when making sdrbf kernel computations and gradients.

kern.nParams = kern.nlfPerInt*kern.nIntervals + kern.nIntervals;
dimParam = [kern.nlfPerInt*kern.nIntervals kern.nIntervals];
kern.inverseWidthIndx = 1:dimParam(1);
kern.switchingTimesIndx = dimParam(1)+1:sum(dimParam(1:2));
%kern.transforms.index = 1:(kern.nlfPerInt*kern.nIntervals);
kern.transforms.index = [1:(kern.nlfPerInt*kern.nIntervals) (kern.nlfPerInt*kern.nIntervals)+2:kern.nParams];
kern.transforms.type = optimiDefaultConstraint('positive');

