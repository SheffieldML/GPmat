function K = multiKernCacheBlock(kern, fhandle, arg, i, j, X, X2)

% MULTIKERNCACHEBLOCK
% FORMAT
% DESC Maintains a cache of values of multiKern blocks with constant
% parameters and returns values
% ARG kern : the structure containing the kernel.
% ARG fhandle : function handle for computing the kernel
% ARG arg : cell array of arguments to fhandle
% ARG i : the row of the block of the kernel to be computed.
% ARG j : the column of the block of the kernel to be computed.
% ARG X1 : first set of kernel inputs.
% ARG X2 : second set of kernel inputs.
% RETURN K : the kernel matrix for the given inputs.
%
% FORMAT
% DESC Clears the cache.
%
% SEEALSO : multiKernCreate, multiKernCompute, multiKernComputeBlock,
% multiKernFixBlocks
%
% COPYRIGHT : Antti Honkela, 2008

% KERN

persistent cache

if nargin == 1,
  cache.(kern.uuid) = cell(kern.numBlocks);
  K = [];
  return;
end

mycache = cache.(kern.uuid){i,j};

key = [X(:); X2(:)];

for k=1:length(mycache),
  if length(key) == length(mycache{k}{1}) && all(key == mycache{k}{1})
    K = mycache{k}{2};
    return;
  end
end

% No match if we get all the way here
if isempty(X2)
  K = fhandle(arg{:}, X);
else
  K = fhandle(arg{:}, X, X2);
end
cache.(kern.uuid){i, j}{end+1} = {key, K};
