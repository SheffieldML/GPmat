function g = simKernDiagGradient(kern, t, covDiag)

% SIMKERNDIAGGRADIENT Compute the gradient of the SIM kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% single input motif kernel matrix with respect to the parameters of the kernel. The
% parameters' gradients are returned in the order given by the
% simKernExtractParam command.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG t : the input times for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% simKernExtractParam.
%
% SEEALSO : simKernParamInit, kernDiagGradient, simKernExtractParam, simKernGradient
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if size(t, 2) > 1
    error('Input can only have one column');
end

sigma = sqrt(2/kern.inverseWidth);

if ~kern.isStationary
    [h, dh_dD, dh_dsigma] = simComputeDiagH(t, kern.decay, kern.delay, sigma);
else
    [h, dh_dD, dh_dsigma] = simComputeDiagStatH(t, kern.decay, sigma);
end

dK_dsigma = 2*dh_dsigma;

if isfield(kern, 'isVarS') && (kern.isVarS)
    K = h;
    if ~kern.isNormalised
        K = sqrt(pi) * K;
        dk_dD = sum(covDiag.*dh_dD)*sqrt(pi)*sigma;
        dk_dsigma = sum(sum(covDiag.*(dK_dsigma*0.5*sqrt(pi)*sigma + K)));
    else
        dk_dD = sum(covDiag.*dh_dD);
        dk_dsigma = 0.5 * sum(sum(covDiag.*dK_dsigma));
    end
    dk_dinvWidth = -0.5*sqrt(2)/(kern.inverseWidth* ...
        sqrt(kern.inverseWidth))*dk_dsigma;
    % only pass the gradient with respect to the inverse width to one
    % of the gradient vectors ... otherwise it is counted twice.
    g = real([dk_dD dk_dinvWidth]);
else
    if isfield(kern, 'isNegativeS') && (kern.isNegativeS)
        C = kern.sensitivity;
    else
        C = sqrt(kern.variance);
    end
    K = h;    
    var2 = C^2;
    if ~kern.isNormalised
        K = sqrt(pi) * K;
        dk_dD = sum(covDiag.*dh_dD)*sqrt(pi)*sigma*var2;
        dk_dsigma = sum(covDiag.*(dK_dsigma*0.5*sqrt(pi)*sigma + K))*var2;
        dk_dC = sigma * C * sum(covDiag.*K);
    else
        dk_dD  = sum(covDiag.*dh_dD)*var2;
        dk_dsigma = 0.5 * var2 * sum(covDiag.*dK_dsigma);
        dk_dC = C * sum(covDiag.*K);
    end
    if isfield(kern, 'isNegativeS') && kern.isNegativeS
        dk_dSimVariance = dk_dC;
    else
        dk_dSimVariance = dk_dC*0.5/C;
    end
    
    dk_dinvWidth = -0.5*sqrt(2)/(kern.inverseWidth* ...
        sqrt(kern.inverseWidth))*dk_dsigma;
    % only pass the gradient with respect to the inverse width to one
    % of the gradient vectors ... otherwise it is counted twice.
    g = real([dk_dD dk_dinvWidth 2*dk_dSimVariance]);
end
end

function [h, dh_dD, dh_dsigma] = simComputeDiagH(t, D, delta, sigma)

t = t - delta;
halfSigmaD = 0.5*sigma*D;
lnPart1 = lnDiffErfs(halfSigmaD + t/sigma, halfSigmaD);
lnPart2 = lnDiffErfs(halfSigmaD, halfSigmaD - t/sigma);
h = exp(halfSigmaD*halfSigmaD   + lnPart1 - log(2*D)) ...
    - exp(halfSigmaD*halfSigmaD - 2*D*t + lnPart2 - log(2*D));
sigma2 = sigma*sigma;
if nargout > 1
    dh_dD = (0.5*D*sigma2*(2*D)-1)*h ...
        + t.*exp(halfSigmaD*halfSigmaD - 2*D*t + lnPart2) ...
        + sigma/sqrt(pi)*(-1 + 2*exp(-t.*t/sigma2-D*t) - exp(-(2*D*t))); 
    dh_dD = real(dh_dD/(2*D));    
    dh_dD_j = t.*exp(halfSigmaD*halfSigmaD-(2*D*t)+lnPart2)-h;
    dh_dD_j = real(dh_dD_j/(2*D));
    dh_dD = dh_dD + dh_dD_j;  
    if nargout > 2
        dh_dsigma = 0.5*D^2*sigma*h + 2/(sqrt(pi)*(2*D))*(-D/2 ...
            + D*exp(-t.*t/sigma2 -D*t) - D/2*exp(-(2*D*t)));
    end
end
end

function [h, dh_dD, dh_dsigma] = simComputeDiagStatH(t, D, sigma)


halfSigmaD = 0.5*sigma*D;
lnPart1 = lnDiffErfs(inf, halfSigmaD);

h = exp(halfSigmaD*halfSigmaD + lnPart1 - log(2*D));

h = h*ones(size(t));

if nargout > 1
    dh_dD = real((-1/(2*D) + sigma*halfSigmaD)*h - sigma/((2*D)*sqrt(pi)));
    dh_dD_j = real(-1/(2*D)*h);
    dh_dD = dh_dD + dh_dD_j;
    if nargout > 2
        dh_dsigma = real(halfSigmaD*D*h - 1/(2*sqrt(pi)));
    end
end

end

% sigma = sqrt(2/kern.inverseWidth);
%
% if (kern.isStationary == false)
%     [h, dh_dD_p, dh_dD_q, dh_dsigma] = ...
%         simComputeH(t, t, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
% else
%     [h, dh_dD_p, dh_dD_q, dh_dsigma] = ...
%         simComputeHStat(t, t, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
% end
%
% if ~isfield(kern, 'isNormalised') || (kern.isNormalised == false)
%     g(1) = sqrt(pi) * sigma * dh_dDp;
%     g(2) = sqrt(pi) * sigma * (dh_dsigma - 1/(2*kern.inverseWidth)*h);
%     g(3) = sqrt(pi) * sigma * h;
% else
%     g(1) = dh_dDp;
%     g(2) = dh_dsigma;
%     g(3) = h;
% end
%
% if isfield(kern, 'isVarS') && (kern.isVarS)
%     g = g(1:2);
% else
%     if isfield(kern, 'isNegativeS') && kern.isNegativeS
%         g(1:2) = g(1:2)*kern.sensitivity*kern.sensitivity;
%         g(3) = g(3)*2*kern.sensitivity;
%     else
%         g(1:2) = g(1:2)*kern.variance;
%     end
%
%     if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
%         g(4) = sum(exp(-2*kern.decay*t) .* covDiag);
%     end
% end
