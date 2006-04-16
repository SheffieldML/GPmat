function kern = kernExpandParam(kern, params)

% KERNEXPANDPARAM Expand parameters to form a kernel structure.
% FORMAT
% DESC returns a kernel structure filled with the parameters in the
% given vector. This is used as a helper function to enable
% parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% As well as extracting the parameters, some transformation of
% parameters is also undertaken in this file. If the field
% transforms is not empty, it dictactes how the kernel parameters
% are to be transformed (for example by a exponential to keep them
% positive).
%
% SEEALSO : kernExtractParam, scg, conjgrad

% KERN

% Check if parameters are being optimised in a transformed space.
if ~isempty(kern.transforms)
  for i = 1:length(kern.transforms)
    index = kern.transforms(i).index;
    fhandle = str2func([kern.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'atox');
  end
end
fhandle = str2func([kern.type 'KernExpandParam']);
kern = fhandle(kern, params);

