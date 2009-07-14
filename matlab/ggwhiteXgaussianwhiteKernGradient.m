function [g1, g2] = ggwhiteXgaussianwhiteKernGradient(ggwhiteKern, gaussianwhiteKern,...
    x, x2, covGrad)
% GGWHITEXGAUSSIANWHITEKERNGRADIENT Compute gradient between the GG white
%                                   and GAUSSIAN white kernels.
% FORMAT
% DESC computes the
%	gradient of an objective function with respect to cross kernel terms
%	between GG white and GAUSSIAN white kernels for the multiple output kernel.
% RETURN g1 : gradient of objective function with respect to kernel
%	   parameters of GG white kernel.
% RETURN g2 : gradient of objective function with respect to kernel
%	   parameters of GAUSSIAN white kernel.
% ARG ggwhitekern : the kernel structure associated with the GG white kernel.
% ARG gaussianwhiteKern :  the kernel structure associated with the GAUSSIAN white kernel.
% ARG x : inputs for which kernel is to be computed.
%
% FORMAT
% DESC  computes
%	the gradient of an objective function with respect to cross kernel
%	terms between GG white and GAUSSIAN white kernels for the multiple output kernel.
% RETURN g1 : gradient of objective function with respect to kernel
%	   parameters of GG white kernel.
% RETURN g2 : gradient of objective function with respect to kernel
%	   parameters of GAUSSIAN white kernel.
% ARG ggwhiteKern : the kernel structure associated with the GG white kernel.
% ARG gaussianwhiteKern : the kernel structure associated with the GAUSSIAN white kernel.
% ARG x1 : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
%
% SEEALSO : multiKernParamInit, multiKernCompute, ggwhiteKernParamInit,
% gaussianwhiteKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009.


% KERN

if ggwhiteKern.isArd ~= gaussianwhiteKern.isArd
    error(['For current implementation of the code, both output kernel' ...
        ' and inducing kernel should be ARD or not']);
end

if nargin < 5
    covGrad = x2;
    x2 = x;
end

[K, P, Pinv, Lrinv, Lqrinv, Kbase, factorVar1, factorNoise, dist] = ...
    ggwhiteXgaussianwhiteKernCompute(ggwhiteKern, gaussianwhiteKern, x, x2);

if ggwhiteKern.isArd
    if gaussianwhiteKern.nIndFunct == 1
        matGradLr = zeros(ggwhiteKern.inputDimension,1);
        matGradLqr = zeros(ggwhiteKern.inputDimension,1);
        for i=1: ggwhiteKern.inputDimension,
            X = repmat(x(:,i),1, size(x2,1));
            X2 = repmat(x2(:,i)',size(x,1),1);
            X_X2 = (X - X2).*(X - X2);
            matGradLr(i) = sum(sum(0.5*covGrad.*K.*...
                (-0.5*Lrinv(i) + Lrinv(i)*P(i)*Lrinv(i) - Lrinv(i)*P(i)*X_X2*P(i)*Lrinv(i))));
            matGradLqr(i) = sum(sum(0.5*covGrad.*K.*...
                (-0.5*Lqrinv(i) + Lqrinv(i)*P(i)*Lqrinv(i) - Lqrinv(i)*P(i)*X_X2*P(i)*Lqrinv(i))));
        end
        grad_sigma2Noise = factorNoise*sum(sum(covGrad.*Kbase));
        grad_variance    = factorVar1*sum(sum(covGrad.*Kbase));
    else
        matGradLr  = zeros(size(gaussianwhiteKern.precisionT));
        matGradLqr = zeros(size(ggwhiteKern.precisionG));
        for i=1:ggwhiteKern.inputDimension
           temp1 = ((P(:,:,i).*Lrinv(:,:,i)).^2).*(Pinv(:,:,i) - dist(:,:,i));
           matGradLr(i, :) = sum(0.5*covGrad.*K.*(-0.5*Lrinv(:,:,i) +temp1), 1)';
            %/~ MAURICIO : This is only important if we do kernTest for this kernel
            %             if size(x,1) >= size(x2,1)
            %                 matGradLr(i,:) = sum(0.5*covGrad.*K.*(-0.5*Lrinv(:,:,i) +temp1), 2)';
            %             else
            %                 matGradLr(i,:) = sum(0.5*covGrad.*K.*(-0.5*Lrinv(:,:,i) +temp1), 1)';
            %             end
            %~/
            temp2 = (Lqrinv(i)^2)*(P(:,:,i).^2).*(Pinv(:,:,i) - dist(:,:,i));
            matGradLqr(i) = sum(sum(0.5*covGrad.*K.*(-0.5*Lqrinv(i) +temp2), 2));
        end
        grad_sigma2Noise = sum(sum(factorNoise.*covGrad.*Kbase));
        grad_variance    = sum(sum(factorVar1.*covGrad.*Kbase));        
    end
else
    if gaussianwhiteKern.nIndFunct == 1        
        temp1 = ((P*Lrinv)^2).*(ggwhiteKern.inputDimension*Pinv - dist);
        matGradLr = sum(sum(0.5*covGrad.*K.*(-0.5*ggwhiteKern.inputDimension*Lrinv + temp1), 1));
        temp2 = (Lqrinv^2)*(P^2)*(ggwhiteKern.inputDimension*Pinv - dist);
        matGradLqr = sum(sum(0.5*covGrad.*K.*(-0.5*ggwhiteKern.inputDimension*Lqrinv +temp2), 2));
        grad_sigma2Noise = sum(sum(factorNoise*covGrad.*Kbase));
        grad_variance    = sum(sum(factorVar1*covGrad.*Kbase));        
    else
        temp1 = ((P.*Lrinv).^2).*(ggwhiteKern.inputDimension*Pinv - dist);
        matGradLr = sum(0.5*covGrad.*K.*(-0.5*ggwhiteKern.inputDimension*Lrinv +temp1), 1)';
        %/~ MAURICIO : This is only important if we do kernTest for this kernel
        %    if size(x,1) >= size(x2,1)
        %        matGradLr = sum(0.5*covGrad.*K.*(-0.5*ggwhiteKern.inputDimension*Lrinv +temp1), 2)';
        %    else
        %        matGradLr = sum(0.5*covGrad.*K.*(-0.5*ggwhiteKern.inputDimension*Lrinv +temp1), 1)';
        %    end
        %~/
        temp2 = (Lqrinv^2)*(P.^2).*(ggwhiteKern.inputDimension*Pinv - dist);
        matGradLqr = sum(sum(0.5*covGrad.*K.*(-0.5*ggwhiteKern.inputDimension*Lqrinv +temp2), 2));
        grad_sigma2Noise = sum(sum(factorNoise.*covGrad.*Kbase));
        grad_variance    = sum(sum(factorVar1.*covGrad.*Kbase));
    end
end
% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g1 = [matGradLqr(:)' 0 grad_variance];
g2 = [matGradLr(:)' grad_sigma2Noise];


