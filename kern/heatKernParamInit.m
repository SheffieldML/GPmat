function kern = heatKernParamInit(kern)

% HEATKERNPARAMINIT HEAT kernel parameter initialisation.
% The heat kernel corresponds to the covariance function for the output
% process of a partial differential equation in two variables, one spatial 
% dimension and time:
%
% \frac{\partial y(x,t)}{\partial t} = Sf(t) - \lambda y(x,t) + D
% \frac{\partial^2 y(x,t)}{\partial x^2}, 
%
% where S is the sensitivity coefficient, \lambda is the deca and D is the
% difussion term. The kernel also contains the inverse widths associated to
% RBF covariance of the Gaussian process prior imposed over f(x,t), one for
% the time variable and one for the spatial variable. By default, the decay 
% and the difussion are restricted to be positive. By default, the
% sensitivity is negative. The particular kernel derived from this equation 
% assumes that the spatial interval is fixed and it can be given as an 
% option. Otherwise, it is assumed to be 0.5mm by default.
%
% Also, the solution to this partial differential equation is given in the
% form of a series. We also specify the number of terms in the series. By
% default the number of terms is 5.
%
% FORMAT
% DESC initialises the heat kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if kern.inputDimension ~= 2
  error('HEAT kernel only valid for 2-D inputs.')
end

if isfield(kern, 'options') && isfield(kern.options, 'nTerms')
    kern.nTerms = kern.options.nTerms;
else
    kern.nTerms = 20;
end

if isfield(kern, 'options') && isfield(kern.options, 'lengthX')
    kern.lengthX = kern.options.lengthX;
else
    kern.lengthX = 0.5;
end

kern.decay = 1;
kern.diffusion = 1e-6;
kern.inverseWidthTime = 100;
kern.inverseWidthSpace = 1000;
kern.sensitivity = 1;
% Create templates for th SIM kernel (used to compute the time part) and
% the SHEAT kernel (used to compute the spatial part).
kern.sim.inputDimension = 1;
if isfield(kern, 'options')
    kern.sim.options = kern.options;
end
kern.sim.options.isVarS = true;
kern.sim = simKernParamInit(kern.sim);
kern.nParams = 5;
kern.transforms.type = optimiDefaultConstraint('positive');
kern.transforms.index = 1:(kern.nParams-1); % The sensitivity could be negative

if isfield(kern, 'options') && isfield(kern.options, 'includeIC') ...
        && kern.options.includeIC
    kern.includeIC = true;
    kern.inverseWidthSpaceIC = 1000;
    kern.sensitivityIC = 1;
    kern.nParams = 7;
    kern.transforms.type = optimiDefaultConstraint('positive');
    kern.transforms.index = [1 2 3 4 6]; % The sensitivities could be negative
else
    kern.includeIC = false;
end

% if isfield(kern, 'options') && isfield(kern.options, 'includeIndSens') ...
%         && kern.options.includeIndSens    
%     kern.includeIndSens = true;
%     kern.sensitivitySpace = 1;
%     if kern.includeIC
%         kern.nParams = 8;
%     else
%         kern.nParams = 6; % The last two are sensitivities for [time space]        
%     end
% else
%     kern.includeIndSens = false;
% end


if isfield(kern, 'options') ...
        && isfield(kern.options, 'isStationary') ...
        && kern.options.isStationary,
   kern.isStationary = true;
else
   kern.isStationary = false;
end

% Type of boundary conditions 

if isfield(kern, 'options') ...
        && isfield(kern.options, 'pde')
   kern.pde = kern.options.pde;
else
   kern.pde = 'sin'; 
end




