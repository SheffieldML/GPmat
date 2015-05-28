function [K, P, Pinv, Lrinv, Lqrinv, Kbase, factorVar1, factorNoise, n2] = ...
    ggwhiteXgaussianwhiteKernCompute(ggwhiteKern, gaussianwhiteKern, x, x2)

% GGWHITEXGAUSSIANWHITEKERNCOMPUTE Compute a cross kernel between the GG white and GAUSSIAN white kernels.
% FORMAT
% DESC computes cross kernel
%	terms between GG white and GAUSSIAN white kernels for the multiple output kernel.
% RETURN k :  block of values from kernel matrix.
% ARG ggKern : the kernel structure associated with the GG kernel.
% ARG gaussianKern : the kernel structure associated with the GAUSSIAN kernel.
% ARG x :  inputs for which kernel is to be computed.
%
% FORMAT
% DESC computes cross
%	kernel terms between GG white and GAUSSIAN white kernels for the multiple output
%	kernel.
% RETURN K : block of values from kernel matrix.
% ARG ggwhiteKern :  the kernel structure associated with the GG kernel.
% ARG gaussianwhiteKern the kernel structure associated with the GAUSSIAN kernel.
% ARG x : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
%
% SEEALSO : multiKernParamInit, multiKernCompute, ggwhiteKernParamInit,
%           gaussianwhiteKernParamInit
%
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009.

% KERN

if ggwhiteKern.isArd ~= gaussianwhiteKern.isArd
    error(['For current implementation of the code, both output kernel' ...
        ' and inducing kernel should be ARD or not']);
end

if nargin < 4
    x2 = x;
end

if ggwhiteKern.isArd
    if gaussianwhiteKern.nIndFunct == 1,
        Lqr = ggwhiteKern.precisionG;
        Lr  = gaussianwhiteKern.precisionT;
        Lqrinv = 1./Lqr;
        Lrinv = 1./Lr;
        Pinv = Lqrinv + Lrinv;
        P = 1./Pinv;
        factorNum = 2^(ggwhiteKern.inputDimension/2)*...
            prod(Lqr)^(-1/4)*prod(Lr)^(-1/4);
        factorDen = prod(Pinv)^(1/2);
        sqrtP = sqrt(P);
        sqrtPx = x*sparseDiag(sqrtP);
        sqrtPx2 = x2*sparseDiag(sqrtP);
        n2 = dist2(sqrtPx, sqrtPx2);
        factor = ggwhiteKern.variance*gaussianwhiteKern.sigma2Noise*...
            factorNum/factorDen;
        Kbase = exp(-0.5*n2);
        K = factor.*Kbase;
        if nargout > 4
            factorVar1 = gaussianwhiteKern.sigma2Noise*factorNum./factorDen;
            factorNoise = ggwhiteKern.variance*factorNum./factorDen;
        end
    else
        kBase = zeros(size(x,1), size(x2,1), ggwhiteKern.inputDimension);
        factorBase = zeros(size(x,1), size(x2,1), ggwhiteKern.inputDimension);
        P = zeros(size(x,1), size(x2,1), ggwhiteKern.inputDimension);
        Pinv = zeros(size(x,1), size(x2,1), ggwhiteKern.inputDimension);        
        n2 = zeros(size(x,1), size(x2,1), ggwhiteKern.inputDimension);        
        Lrinv = zeros(size(x,1), size(x2,1), ggwhiteKern.inputDimension);        
        Lqrinv = zeros(size(ggwhiteKern.precisionG));
        for i=1:ggwhiteKern.inputDimension,
            n2(:,:,i) = dist2(x(:,i), x2(:,i));
            Lqr = ggwhiteKern.precisionG(i);
            if gaussianwhiteKern.nIndFunct~=size(x2,1)
                error(['The number of inducing functions must be equal the' ...
                    'number of inducing points']);
            end
            Lr = repmat(gaussianwhiteKern.precisionT(i,:), size(x,1), 1);
            %/~MAURICIO : This code is necessary if we are doing kernTest of this multioutput Kernel
            %             if nargin < 4
            %                if gaussianwhiteKern.nIndFunct~=size(x2,1)
            %                    error(['The number of inducing functions must be equal the' ...
            %                        'number of inducing points']);
            %                end
            %                Lr = repmat(gaussianwhiteKern.precisionT(i,:), size(x,1), 1);
            %             else
            %                if size(x,1) >= size(x2,1)
            %                    Lr = repmat(gaussianwhiteKern.precisionT(i,:)', 1, size(x2,1));
            %                else
            %                    Lr = repmat(gaussianwhiteKern.precisionT(i,:), size(x,1), 1);
            %                end
            %             end 
            %~/
            Lqrinv(i) = 1/Lqr;
            Lrinv(:,:,i) = 1./Lr;
            Pinv(:,:,i) = Lqrinv(i) + Lrinv(:,:,i);
            P(:,:,i) = 1./Pinv(:,:,i);
            factorNum = 2^(1/2)*(Lqr^(-1/4))*(Lr.^(-1/4));
            factorDen = Pinv(:,:,i).^(1/2);
            factorBase(:,:,i) = factorNum./factorDen;
            kBase(:,:,i) = exp(-0.5.*P(:,:,i).*n2(:,:,i));
        end
        Kbase = kBase(:,:,1);
        factor = factorBase(:,:,1);
        for i = 2:ggwhiteKern.inputDimension,
           Kbase = Kbase.*kBase(:,:,i);
           factor = factor.*factorBase(:,:,i); 
        end
        K = (gaussianwhiteKern.sigma2Noise*ggwhiteKern.variance)*factor.*Kbase;
        if nargout > 4
            factorVar1 = gaussianwhiteKern.sigma2Noise*factor;
            factorNoise = ggwhiteKern.variance*factor;
        end
    end
else
    if gaussianwhiteKern.nIndFunct == 1,
        n2 = dist2(x, x2);
        Lqr = ggwhiteKern.precisionG;
        Lr = gaussianwhiteKern.precisionT;
        Lqrinv = 1/Lqr;
        Lrinv = 1/Lr;
        Pinv = Lqrinv + Lrinv;
        P = 1/Pinv;
        factorNum = 2^(ggwhiteKern.inputDimension/2)*...
            (Lqr^(-ggwhiteKern.inputDimension/4))*...
            (Lr^(-ggwhiteKern.inputDimension/4));
        factorDen = Pinv^(ggwhiteKern.inputDimension/2);
        factor = (gaussianwhiteKern.sigma2Noise*...
            ggwhiteKern.variance)*factorNum/factorDen;
        Kbase = exp(-0.5*P*n2);
        K = factor*Kbase;
        if nargout > 4
            factorVar1 = gaussianwhiteKern.sigma2Noise*factorNum/factorDen;
            factorNoise = ggwhiteKern.variance*factorNum/factorDen;
        end
    else
        n2 = dist2(x, x2);
        Lqr = ggwhiteKern.precisionG;
        if gaussianwhiteKern.nIndFunct~=size(x2,1)
            error(['The number of inducing functions must be equal the' ...
                'number of inducing points']);
        end
        Lr = repmat(gaussianwhiteKern.precisionT, size(x,1), 1);
        %/~MAURICIO : This code is necessary if we are doing kernTest of this multioutput Kernel
        %if nargin < 4
        %    if gaussianwhiteKern.nIndFunct~=size(x2,1)
        %        error(['The number of inducing functions must be equal the' ...
        %            'number of inducing points']);
        %    end
        %    Lr = repmat(gaussianwhiteKern.precisionT, size(x,1), 1);
        %else
        %    if size(x,1) >= size(x2,1)
        %        Lr = repmat(gaussianwhiteKern.precisionT', 1, size(x2,1));
        %    else
        %        Lr = repmat(gaussianwhiteKern.precisionT, size(x,1), 1);
        %    end
        %end 
        %~/
        Lqrinv = 1/Lqr;
        Lrinv = 1./Lr;
        Pinv = Lqrinv + Lrinv;
        P = 1./Pinv;
        factorNum = 2^(ggwhiteKern.inputDimension/2)*...
            (Lqr^(-ggwhiteKern.inputDimension/4))*...
            (Lr.^(-ggwhiteKern.inputDimension/4));
        factorDen = Pinv.^(ggwhiteKern.inputDimension/2);
        factor = (gaussianwhiteKern.sigma2Noise*...
            ggwhiteKern.variance)*factorNum./factorDen;
        Kbase = exp(-0.5.*P.*n2);
        K = factor.*Kbase;
        if nargout > 4
            factorVar1 = gaussianwhiteKern.sigma2Noise*factorNum./factorDen;
            factorNoise = ggwhiteKern.variance*factorNum./factorDen;
        end
    end
end

