function k = multiKernDiagCompute(kern, x)

% MULTIKERNDIAGCOMPUTE Compute diagonal of MULTI kernel.
%
%	Description:
%
%	K = MULTIKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the multiple output block kernel given a design
%	matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%	
%
%	See also
%	MULTIKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, MULTIKERNCOMPUTE


%	Copyright (c) 2006 Neil D. Lawrence
%	Copyright (c) 2007 Pei Gao

% MODIFICATIONS : Mauricio Alvarez, 2008


if iscell(x)
  dim = 0;
  for i = 1:length(x)
    dim = dim + size(x{i},1);
  end
  k = zeros(dim, 1);
  startVal = 1;
  endVal = size(x{1},1);
  for i = 1:length(kern.comp)
      if ~isempty(x{i})
          k(startVal:endVal) = kernDiagCompute(kern.comp{i}, x{i});
      end
      startVal = endVal + 1;
      if i+1 <= length(kern.comp)
          endVal = endVal + size(x{i+1},1);
      end
  end
else
  k = zeros(size(x, 1)*kern.numBlocks, 1);
  endVal = size(x, 1);
  startVal = 1;
  for i = 1:length(kern.comp)
    k(startVal:endVal, 1)  = kernDiagCompute(kern.comp{i}, x);
    startVal = endVal + 1;
    endVal = endVal + size(x, 1);
  end
end

