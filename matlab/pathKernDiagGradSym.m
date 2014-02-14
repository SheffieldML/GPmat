function [gX, kern] = pathKernDiagGradSym(kern, X)

% PATHKERNDIAGGRADSYM Compute the gradient of the PATH kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% path kernel matrix with respect to the parameters of the
% symmetric kernel.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG x : the input data for which the gradient is being computed.
% RETURN gX : gradients of the relevant function with respect to each
% of the parameters. 
% RETURN kern : the updated kernel structure
%
% SEEALSO : pathKernParamInit, pathDiagGradient, pathKernExtractParam, pathKernGradient
%
% COPYRIGHT : Andrea Baisero, Carl Henrik Ek, 2013

% SHEFFIELDML


maxl=max(cellfun(@(x)size(x,1),X));
kern=pathKernUpdateWMat(kern,maxl);

num=length(X);

gX=cell(1,num);
for i=1:num
    gX{i}=pathKernDiagGradSymSeq(kern,X{i});
end

function gX = pathKernDiagGradSymSeq(kern, x)

l=size(x,1);
dim=size(x,2);

t=kernGradX(kern.gkern,x,x);
kw=kern.wmat(1:l,1:l);

gX=zeros(size(x));
for i=1:l
    gX(i,:)=sum((kw(i,:)'*ones(1,dim)).*t(:,:,i));
end
