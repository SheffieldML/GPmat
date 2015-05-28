function k = cmpndKernDiagCompute(kern, x)

% CMPNDKERNDIAGCOMPUTE Compute diagonal of CMPND kernel.
%
%	Description:
%
%	K = CMPNDKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the compound kernel given a design matrix of
%	inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	CMPNDKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, CMPNDKERNCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



i = 1;
if ~isempty(kern.comp{i}.index)
  % only part of the data is involved with the kernel.
  k  = kernDiagCompute(kern.comp{i}, x(:, kern.comp{i}.index));
else
  % all the data is involved with the kernel.
  k  = kernDiagCompute(kern.comp{i}, x);
end
for i = 2:length(kern.comp)
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved with the kernel.
    k  = k + kernDiagCompute(kern.comp{i}, x(:, kern.comp{i}.index));
  else
    % all the data is involved with the kernel.
    k  = k + kernDiagCompute(kern.comp{i}, x);
  end
end
