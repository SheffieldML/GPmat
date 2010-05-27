function  kern = sdlfmKernParamInit(kern)

% SDLFMKERNPARAMINIT SDLFM kernel initialization
% FORMAT
% DESC
% Initializes the switching dynamical latent force model kernel structure 
% with some initial parameters. The initial parameters are passed through 
% an option in kern.
% RETURN kern : the kernel structure with the default parameters placed in.
% ARG kern : the kernel structure which requires initialisation.
%
% SEEALSO : kernCreate, kernParamInit, lfmKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

% Create the basic structure based on the lfm kernel.
kern = lfmKernParamInit(kern);

% Remove unnecessary fields from the LFM structure
kern = rmfield(kern, {'delay','initVal','variance', 'serialNumber'});

if isfield(kern, 'options') && isfield(kern.options, 'nIntervals')
    kern.nIntervals = kern.options.nIntervals;
    if isfield(kern.options, 'nlfPerInt')
        kern.nlfPerInt = kern.options.nlfPerInt;
    else
        kern.nlfPerInt = 1;
    end
    kern.inverseWidth = ones(kern.nlfPerInt, kern.nIntervals); % Total number of inverse widths.
    kern.sensitivity = ones(kern.nlfPerInt, kern.nIntervals); % Total number of sensitivities.
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

% Number of parameters computed as:
% mass + damper + spring + inverseWidth (= nlf*intervals)
% + switching times + sensitivities( = nlf*intervals) 
dimParam = [3 kern.nlfPerInt*kern.nIntervals kern.nIntervals ...
    kern.nlfPerInt*kern.nIntervals];
kern.nParams = sum(dimParam);
kern.outputIndx = [1 2 3]; % Indexex of the mass, damper and decay
kern.inverseWidthIndx = dimParam(1)+1:sum(dimParam(1:2));
kern.switchingTimesIndx = sum(dimParam(1:2))+1:sum(dimParam(1:3));
kern.sensitivityIndx = sum(dimParam(1:3))+1:sum(dimParam);
kern.transforms.index = [1:sum(dimParam(1:2)) (sum(dimParam(1:2))+2):sum(dimParam(1:3))] ;
%kern.transforms.index = 1:sum(dimParam(1:3));
kern.transforms.type = optimiDefaultConstraint('positive');



