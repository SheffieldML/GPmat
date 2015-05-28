function k = ndsimKernDiagCompute(kern, t)

% NDSIMKERNDIAGCOMPUTE Compute diagonal of NDSIM kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the single input motif kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : simKernParamInit, kernDiagCompute, kernCreate, simKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Antti Honkela, 2008
%
% MODIFICATIONS : David Luengo, 2009
%
% MODIFICATIONS : Mauricio Alvarez, 2009
%
% COPYRIGHT : Jaakko Peltonen, 2011

% KERN

if size(t, 2) > 1 
  error('Input can only have one column');
end


sigma = sqrt(2/kern.inverseWidth);
if isfield(kern, 'isNegativeS') && (kern.isNegativeS == true)
    variancemultiplier = (kern.sensitivity*kern.sensitivity);
else
    variancemultiplier = kern.variance;
end

% dim1 = size(t1, 1);
% dim2 = size(t2, 1);
% t1Mat = t1(:, ones(1, dim2));
% t2Mat = t2(:, ones(1, dim1))';
% diffT = (t1Mat - t2Mat);

k1 = 2*t.*erf(t/sigma);
k2 = k1;
k3 = 0;
k4 = +2*sigma/sqrt(pi)*(exp(-(t/sigma).^2));
k5 = k4;
k6 = +2*sigma/sqrt(pi)*(-1-1);

k = (k1+k2+k3+k4+k5+k6)*sqrt(pi)*sigma/4*variancemultiplier;


