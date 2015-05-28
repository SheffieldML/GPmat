function [g, kern] = pathKernGradient(kern, x, varargin)

% PATHKERNGRADIENT Gradient of PATH kernel's parameters.
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


if nargin<4
    maxl=max(cellfun(@(x)size(x,1),x));
else
    maxl=max(cellfun(@(x)size(x,1),[x,varargin{1}]));
end
kern=pathKernUpdateWMat(kern,maxl);

if nargin<4
    num=length(x);

    gwd=zeros(num);
    gwhv=zeros(num);
    ggk=zeros(num,num,kern.gkern.nParams);
    for i=1:num
        li=size(x{i},1);
        w=kern.wmat(1:li,1:li);
        w=.5*(w+w(end:-1:1,end:-1:1));
        dwd=kern.dwdmat(1:li,1:li);
        dwd=.5*(dwd+dwd(end:-1:1,end:-1:1));
        dwhv=kern.dwhvmat(1:li,1:li);
        dwhv=.5*(dwhv+dwhv(end:-1:1,end:-1:1));

        gk=kernCompute(kern.gkern,x{i});
        gwd(i,i)=sum(sum(dwd.*gk));
        gwhv(i,i)=sum(sum(dwhv.*gk));
        ggk(i,i,:)=kernGradient(kern.gkern,x{i},w);

        for j=i+1:num
            lj=size(x{j},1);
            w=kern.wmat(1:li,1:lj);
            w=.5*(w+w(end:-1:1,end:-1:1));
            dwd=kern.dwdmat(1:li,1:lj);
            dwd=.5*(dwd+dwd(end:-1:1,end:-1:1));
            dwhv=kern.dwhvmat(1:li,1:lj);
            dwhv=.5*(dwhv+dwhv(end:-1:1,end:-1:1));

            gk=kernCompute(kern.gkern,x{i},x{j});
            gwd(i,j)=sum(sum(dwd.*gk));
            gwd(j,i)=gwd(i,j);
            gwhv(i,j)=sum(sum(dwhv.*gk));
            gwhv(j,i)=gwhv(i,j);
            ggk(i,j,:)=kernGradient(kern.gkern,x{i},x{j},w);
            ggk(j,i,:)=ggk(i,j,:);
        end
    end
else
    x2=varargin{1};
    
    num=length(x);
    num2=length(x2);

    gwd=zeros(num,num2);
    gwhv=zeros(num,num2);
    ggk=zeros(num,num2,kern.gkern.nParams);
    for i=1:num
        li=size(x{i},1);
        for j=1:num2
            lj=size(x2{j},1);
            w=kern.wmat(1:li,1:lj);
            w=.5*(w+w(end:-1:1,end:-1:1));
            dwd=kern.dwdmat(1:li,1:lj);
            dwd=.5*(dwd+dwd(end:-1:1,end:-1:1));
            dwhv=kern.dwhvmat(1:li,1:lj);
            dwhv=.5*(dwhv+dwhv(end:-1:1,end:-1:1));

            gk=kernCompute(kern.gkern,x{i},x2{j});
            gwd(i,j)=sum(sum(dwd.*gk));
            gwhv(i,j)=sum(sum(dwhv.*gk));
            ggk(i,j,:)=kernGradient(kern.gkern,x{i},x2{j},w);
        end
    end
end

g(1,1)=sum(sum(varargin{end}.*gwd));
g(1,2)=sum(sum(varargin{end}.*gwhv));
for i=1:kern.gkern.nParams
    g(1,i+2)=sum(sum(varargin{end}.*ggk(:,:,i)));
end
