function K = multiKernComputeBlock(kern, X, X2, i, j)

% MULTIKERNCOMPUTEBLOCK
% FORMAT
% DESC computes a block of a multi-output kernel given two matrices
% of input. 
% ARG kern : the structure containing the kernel.
% ARG X1 : first set of kernel inputs.
% ARG X2 : second set of kernel inputs.
% ARG i : the row of the block of the kernel to be computed.
% ARG j : the column of the block of the kernel to be computed.
% RETURN K : the kernel matrix for the given inputs.
%
% FORMAT
% DESC computes a block of a multi-output kernel given a single matrix
% of input.
% ARG kern : the structure containing the kernel.
% ARG X : first set of kernel inputs.
% ARG i : the row of the block of the kernel to be computed.
% ARG j : the column of the block of the kernel to be computed.
% RETURN K : the kernel matrix for the given inputs.
%
% SEEALSO : multiKernCreate, multiKernCompute, multiKernGradientBlock
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS: Antti Honkela, 2008

% KERN


%fprintf(1,'multiKernComputeBlock step1\n');
%X
%X2
%i
%j
%nargin


if nargin < 5
  j = i;
  i = X2;
  X2 = [];
end  
if i == j
  fhandle = [kern.comp{i}.type 'KernCompute'];
  transpose = 0;
  arg{1} = kern.comp{i};
else
  if j<i
    fhandle = [kern.block{i}.cross{j} 'KernCompute'];
    transpose = kern.block{i}.transpose(j);
  else
    fhandle = [kern.block{j}.cross{i} 'KernCompute'];
    transpose = ~kern.block{j}.transpose(i);
  end
  if transpose
    arg{1} = kern.comp{j};
    arg{2} = kern.comp{i};
  else
    arg{1} = kern.comp{i};
    arg{2} = kern.comp{j};
  end
end
fhandle = str2func(fhandle);

if isfield(kern, 'fixedBlocks') && kern.fixedBlocks(i) && ...
      kern.fixedBlocks(j),
  K = multiKernCacheBlock(kern, fhandle, arg, i, j, X, X2);
else
  % PROBLEM: inconsistent behavior!

  if isempty(X2)
    K = fhandle(arg{:}, X);
  else
    K = fhandle(arg{:}, X, X2);
  end
end
%if transpose
%  K = K';
%end
