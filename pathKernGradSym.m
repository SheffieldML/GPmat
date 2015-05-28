function [gX, kern] = pathKernGradSym(kern, X, X2)

% PATHKERNGRADSYM Gradient of the symetric PATH kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% path
% kernel's parameters. As well as the kernel structure and the
% input positions, the user provides a matrix PARTIAL which gives
% the partial derivatives of the function with respect to the
% relevant elements of the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x : the input locations for which the gradients are being
% computed. 
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in X.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters. The ordering of the vector should match
% that provided by the function kernExtractParam.
% RETURN kern : the updated kernel structure
%
% FORMAT
% DESC computes the derivatives as above, but input locations are
% now provided in two matrices associated with rows and columns of
% the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x1 : the input locations associated with the rows of the
% kernel matrix.
% ARG x2 : the input locations associated with the columns of the
% kernel matrix.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as X1 and the same number of columns
% as X2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
% RETURN kern : the updated kernel structure
%
% SEEALSO pathKernParamInit, kernGradient, pathKernDiagGradient
%
% COPYRIGHT : Andrea Baisero, Carl Henrik Ek, 2013

% SHEFFIELDML


maxl=max(cellfun(@(x)size(x,1),[X,X2]));
kern=pathKernUpdateWMat(kern,maxl);

num=length(X);
num2=length(X2);

gX=cell(num,num2);
for i=1:num
    for j=1:num2
        gX{i,j}=pathKernGradSymSeq(kern,X{i},X2{j});
    end
end

function gX = pathKernGradSymSeq(kern, x, x2)

l=size(x,1);
l2=size(x2,1);
dim=size(x,2);

t=kernGradX(kern.gkern,x,x2);
kw=kern.wmat(1:l,1:l2);

gX=zeros(size(x));
for i=1:l
    gX(i,:)=sum((kw(i,:)'*ones(1,dim)).*t(:,:,i));
end
