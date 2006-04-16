function kern = cmpndKernParamInit(kern)

% CMPNDKERNPARAMINIT CMPND kernel parameter initialisation.
% The compound (CMPND) kernel is a container kernel for allowing
% several different kernels to be added together. It is created by
% using the kernCreate command with the kernel type given as a
% cell. For example, to create a compound kernel that is composed
% of an RBF kernel, a LIN kernel and a WHITE kernel you call
%
% kern = kernCreate(X, {'rbf', 'lin', 'white'});
%
% Each individual kernel is then stored within the returned kernel
% structure. The kernels are stored in order in a field called
% 'comp'. So display obtain the 'rbf' kernel you write:
%
% kernDisplay(kern.comp{1})
%
% Optionally the first argument of the cell can be 'cmpnd'. This is
% to differentiate the call from other possible container calls
% such as 'tensor' and 'multi'.
%
% Note that there is the option to 'tie' CMPND kernel parameters
% together, so that they are optimised as one parameter. Which kernel
% parameters are tied together are identified in the 'paramGroups'
% field.
%
% Finally, note that there is a field 'whiteVariance' which
% summarises the white noise terms from all the kernels that make
% up the compound kernel. This is useful in some code (for example
% the IVM) where a quick assesment of the noise is sometimes
% required.
%
% SEEALSO : tensorKernParamInit, multiKernParamInit
%
% FORMAT
% DESC initialises the compound
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


kern.nParams = 0;
kern.transforms = [];
if ~isfield(kern, 'comp')
  kern.comp=cell(0);
end
for i = 1:length(kern.comp)
  kern.comp{i} = kernParamInit(kern.comp{i});
  kern.nParams = kern.nParams + kern.comp{i}.nParams;
  kern.comp{i}.index = [];
end
kern.paramGroups = speye(kern.nParams);

% Summarise the total white variance in the field whiteVariance.
kern.whiteVariance = 0;
for i = 1:length(kern.comp)
  if strcmp(kern.comp{i}.type, 'white')
    kern.whiteVariance = kern.whiteVariance + kern.comp{i}.variance;
  else
    if(isfield(kern.comp{i}, 'whiteVariance'))
      kern.whiteVariance = kern.whiteVariance + ...
          kern.comp{i}.whiteVariance;
    end
  end
end