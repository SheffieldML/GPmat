function kern = invcmpndKernParamInit(kern)

% INVCMPNDKERNPARAMINIT INV.PRECISION-CMPND kernel parameter initialisation.
% The inverse precision compound (INVCMPND) kernel (call it K) is a container kernel for
% allowing several different kernels (K1, K2, ...) to be considered together through the
% formula inv( inv(K1) + inv(K2) + ... ), i.e. the inverse of the sum of
% the precision matrices of the individual kernels. That's the result one
% gets as a cov. matrix by taking the product of Gaussians.
% It is created by using the kernCreate command with the kernel type given as a
% cell. The kernel type is usually itself a (normal) COMPOUND kernel, since
% each Ki being a kernel {'...', 'white', 'bias'} guarantees that it is
% numerically stable and in practise positive definite, so that K will also
% be positive definite. Also, each Ki is allowed to take as input different
% partitions of X. For this reason, for creating the kernel X must be
% passed as a cell (one cell element for each Ki), and then X will be
% concatenated into a matrix and the relevant indices of each Ki will be
% found in Ki.index.

% For example, to create a compound kernel K that is composed
% of K1 = an RBF, white and bias kernel with input X1, and
% K2 = a matern32, white and bias kernel with input X2, you call:
%
% kern = kernCreate({X1, X2}, {'invcmpnd',{'rbf', 'white', 'bias'}, {'matern32', 'white', 'bias'}});
%
% If each Ki takes the same input X, instead of {X1, X2} you can simply pass X.
% Each individual kernel is then stored within the returned kernel
% structure. The kernels are stored in order in a field called
% 'comp'. So display obtain the {'rbf', 'white', 'bias'} kernel you write:
%
% kernDisplay(kern.comp{1})
%
% To see which input indices Ki is using, you can check
% kern.comp{1}.index
%
% Necessarily, the first argument of the cell can be 'invcmpnd'. This is
% to differentiate the call from other possible container calls
% such as 'cmpnd', 'tensor' and 'multi'.
%
% Note that there is the option to 'tie' CMPND kernel parameters
% together, so that they are optimised as one parameter. Which kernel
% parameters are tied together are identified in the 'paramGroups'
% field. 
% !!!!! This has not yet be tested well in the invcmpnd kernel.
%
%
% SEEALSO : cmpndKernParamInit, kernCreate, tensorKernParamInit, multiKernParamInit
%
% FORMAT
% DESC initialises the inverse precision compound
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Andreas C. Damianou, 2012

% KERN



% We don't use the following, because it sets the indexes to [].
% kern = cmpndKernParamInit(kern);

kern.jitter = 1e-5;

kern.nParams = 0;
kern.transforms = [];
if ~isfield(kern, 'comp')
  kern.comp=cell(0);
end
for i = 1:length(kern.comp)
  
  % Set up the component kernels.
  %kern.comp{i} = kernParamInit(kern.comp{i});
  kern.nParams = kern.nParams + kern.comp{i}.nParams;
  %kern.comp{i}.index = [];
  
  % Check whether the kernel is a multi-output block kernel.
  if isfield(kern.comp{i}, 'numBlocks')
    if i == 1
      kern.numBlocks = kern.comp{i}.numBlocks;
    else
      if ~isfield(kern, 'numBlocks') | kern.numBlocks ~= kern.comp{i}.numBlocks
        error('Compound of multi kerns with different numbers of blocks')
      end
    end
  else
    if isfield(kern, 'numBlocks')
      error('Attempt to combine multi-kernel with non multi kernel.')
    end
  end
end
kern.paramGroups = speye(kern.nParams);

% Summarise the total white variance in the field whiteVariance and
% find out whether the kernel is stationary.
kern.whiteVariance = 0;
kern.isStationary = true;
for i = 1:length(kern.comp)
  if ~kern.comp{i}.isStationary
    kern.isStationary = false;
  end
  if strcmp(kern.comp{i}.type, 'white')
    kern.whiteVariance = kern.whiteVariance + kern.comp{i}.variance;
  else
    if(isfield(kern.comp{i}, 'whiteVariance'))
      kern.whiteVariance = kern.whiteVariance + ...
          kern.comp{i}.whiteVariance;
    end
  end
end




