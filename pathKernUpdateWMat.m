function kern = pathKernUpdateWMat(kern,dim,force_compute)

% PATHKERNUPDATEWMAT Update the weight matrix for the PATH kernel
% FORMAT
% DESC updates the weight matrix for the path kernel 
% ARG kern : the kernel structure
% RETURN kern : the kernel structure with the updated matrix
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Andrea Baisero, Carl Henrik Ek, 2013

% SHEFFIELDML


if(nargin<3)
    force_compute = false;
    if(nargin<2)
        dim = [];
    end
end

if(force_compute)
    if(~isempty(dim))
        kern.wmatDim = dim;
    else
        if(kern.wmatDim<=0)
            kern.wmatDim = 1;
        end
    end
    [w,dwd,dwhv]=pathWMatPos(kern.wmatDim,kern.cd,kern.chv);
    kern.wmat=w;
    kern.dwdmat=dwd;
    kern.dwhvmat=dwhv;
else
    if(~isempty(dim))
        if(kern.wmatDim<dim)
            kern.wmatDim=dim;
            [w,dwd,dwhv]=pathWMatPos(dim,kern.cd,kern.chv);
            kern.wmat=w;
            kern.dwdmat=dwd;
            kern.dwhvmat=dwhv;
        end
    end
end
