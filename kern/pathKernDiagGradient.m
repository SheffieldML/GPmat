function [g, kern] = pathKernDiagGradient(kern, x, covDiag)

% PATHKERNDIAGGRADIENT Compute the gradient of the PATH kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% path kernel matrix with respect to the parameters of the kernel. The
% parameters' gradients are returned in the order given by the
% pathKernExtractParam command.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG x : the input data for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% pathKernExtractParam.
% RETURN kern : the updated kernel structure
%
% SEEALSO : pathKernParamInit, pathDiagGradient, pathKernExtractParam, pathKernGradient
%
% COPYRIGHT : Andrea Baisero, Carl Henrik Ek, 2013

% SHEFFIELDML


maxl=max(cellfun(@(x)size(x,1),x));
kern=pathKernUpdateWMat(kern,maxl);

num=length(x);

gwd=zeros(1,num);
gwhv=zeros(1,num);
ggk=zeros(num,kern.gkern.nParams);
for i=1:num
    li=size(x{i},1);
    w=kern.wmat(1:li,1:li);
    w=.5*(w+w(end:-1:1,end:-1:1));
    dwd=kern.dwdmat(1:li,1:li);
    dwd=.5*(dwd+dwd(end:-1:1,end:-1:1));
    dwhv=kern.dwhvmat(1:li,1:li);
    dwhv=.5*(dwhv+dwhv(end:-1:1,end:-1:1));

    gk=kernCompute(kern.gkern,x{i});
    gwd(i)=sum(sum(dwd.*gk));
    gwhv(i)=sum(sum(dwhv.*gk));
    ggk(i,:)=kernGradient(kern.gkern,x{i},w);
end

g(1,1)=sum(covDiag.*gwd);
g(1,2)=sum(covDiag.*gwhv);
for i=1:kern.gkern.nParams
    g(1,i+2)=sum(covDiag.*ggk(:,i));
end
