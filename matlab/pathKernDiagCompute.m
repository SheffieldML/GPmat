function [k, gk, kern] = pathKernDiagCompute(kern, x)

% PATHKERNDIAGCOMPUTE Compute diagonal of PATH kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the path kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% RETURN gk : evaluation of ground kernel
% RETURN kern: updated kernel structure
% computed at the given points.
%
% SEEALSO : pathKernParamInit, kernDiagCompute, kernCreate, pathKernCompute
%
% COPYRIGHT : Andrea Baisero, Carl Henrik Ek, 2013

% SHEFFIELDML


maxl=max(cellfun(@(x)size(x,1),x));
kern=pathKernUpdateWMat(kern,maxl);

num=length(x);

gk=cell(1,num);
k=zeros(1,num);
for i=1:num
    li=size(x{i},1);
    w=kern.wmat(1:li,1:li);
    w=.5*(w+w(end:-1:1,end:-1:1));

    gk{i}=kernCompute(kern.gkern,x{i});
    k(i)=sum(sum(w.*gk{i}));
end

