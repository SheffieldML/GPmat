function k = ndsimKernCompute(kern, t1, t2)

% NDSIMKERNCOMPUTE Compute the NDSIM kernel with no decay given the parameters and X.
%
% FORMAT
% DESC computes the kernel parameters for the single input motif
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t1 : the input matrix associated with the rows of the kernel.
% ARG t2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
% RETURN sk : unscaled kernel matrix (i.e. only 0.5 times the sum of h's
% part).
%
% FORMAT
% DESC computes the kernel matrix for the single input motif
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
% RETURN sk : unscaled kernel matrix (i.e. only 0.5 times the sum of h's
% part).
%
% SEEALSO : simKernParamInit, kernCompute, kernCreate, simKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : David Luengo, 2009
%
% MODIFICATIONS : Mauricio Alvarez, 2009
%
% MODIFICATIONS : Antti Honkela, 2009
%
% COPYRIGHT : Jaakko Peltonen, 2011

% KERN

if nargin < 3
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

sigma = sqrt(2/kern.inverseWidth);
if isfield(kern, 'isNegativeS') && (kern.isNegativeS == true)
    variancemultiplier = (kern.sensitivity*kern.sensitivity);
else
    variancemultiplier = kern.variance;
end

dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);

k1 = 2*t1Mat.*erf(t1Mat/sigma);
k2 = 2*t2Mat.*erf(t2Mat/sigma);
k3 = - 2*diffT.*erf(diffT/sigma);
k4 = +2*sigma/sqrt(pi)*(exp(-(t1Mat/sigma).^2));
k5 = +2*sigma/sqrt(pi)*(exp(-(t2Mat/sigma).^2));
k6 = +2*sigma/sqrt(pi)*(-exp(-(diffT/sigma).^2)-1);

%t1Mat
%t2Mat
%diffT
%pause
%k1
%k2
%k3
%pause
%k4
%k5
%k6
%pause

k = (k1+k2+k3+k4+k5+k6)*sqrt(pi)*sigma/4*variancemultiplier;


%if isfield(kern, 'isNormalised') && (kern.isNormalised == true)
%TODO    k = sk;
%end;

% gaussianInitial currently unsupported
