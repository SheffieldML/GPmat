function [k, gk, kern] = pathKernCompute(kern, x, x2)

% PATHKERNCOMPUTE Compute the Path kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the path
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
% RETURN gk : the evaluation of the ground kernel
% RETURN kern : the updated kernel
%
% FORMAT
% DESC computes the kernel parameters for the path
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix in the form of a design matrix
% RETURN k : the kernel matrix computed at the given points.
% RETURN gk : the evaluation of the ground kernel
% RETURN kern : the updated kernel
%
% SEEALSO : pathKernParamInit, kernCompute, kernCreate, pathKernCreate
%
% COPYRIGHT : Andrea Baisero, Carl Henrik Ek, 2013
%
% 

% SHEFFIELDML


if nargin<3
    maxl=max(cellfun(@(x)size(x,1),x));
    kern=pathKernUpdateWMat(kern,maxl);

    num=length(x);

    gk=cell(num);
    k=zeros(num);
    for i=1:num
        li=size(x{i},1);
        w=kern.wmat(1:li,1:li);
        w=.5*(w+w(end:-1:1,end:-1:1));

        gk{i,i}=kernCompute(kern.gkern,x{i},x{i});
        k(i,i)=sum(sum(w.*gk{i,i}));

        for j=i+1:num
            lj=size(x{j},1);
            w=kern.wmat(1:li,1:lj);
            w=.5*(w+w(end:-1:1,end:-1:1));

            gk{i,j}=kernCompute(kern.gkern,x{i},x{j});
            gk{j,i}=gk{i,j};
            k(i,j)=sum(sum(w.*gk{i,j}));
            k(j,i)=k(i,j);
        end
    end
else
    maxl=max(max(cellfun(@(x)size(x,1),x)),max(cellfun(@(x2)size(x2,1),x2)));
    kern=pathKernUpdateWMat(kern,maxl);

    num=length(x);
    num2=length(x2);

    gk=cell(num,num2);
    k=zeros(num,num2);
    for i=1:num
        li=size(x{i},1);
        for j=1:num2
            lj=size(x2{j},1);
            w=kern.wmat(1:li,1:lj);
            w=.5*(w+w(end:-1:1,end:-1:1));

            gk{i,j}=kernCompute(kern.gkern,x{i},x2{j});
            k(i,j)=sum(sum(w.*gk{i,j}));
        end
    end
end

