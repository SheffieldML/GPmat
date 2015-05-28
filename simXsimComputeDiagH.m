function [h, dh_dD_i, dh_dD_j, dh_dsigma] = simXsimComputeDiagH(t, D_i, D_j, delta_i, delta_j, sigma)

% SIMXSIMCOMPUTEDIAGH Helper function for comptuing part of the SIM kernel.
% FORMAT
% DESC computes a portion of the SIM kernel.
% ARG t :  time input (number of time points x 1).
% ARG decay1 : Decay rate for first system.
% ARG decay2 : Decay rate for second system.
% ARG delay1 : Delay for first system.
% ARG delay2 : Delay for second system.
% ARG l : length scale of latent process.
% RETURN h : result of this subcomponent of the kernel for the given values.
%
% DESC computes a portion of the SIM kernel and gradients with
% respect to various parameters.
% ARG t : time input (number of time points x 1).
% ARG decay1 : Decay rate for first system.
% ARG decay2 : Decay rate for second system.
% ARG delay1 : Delay for first system.
% ARG delay2 : Delay for second system.
% ARG l : length scale of latent process.
% RETURN h : result of this subcomponent of the kernel for the given values.
% RETURN grad_D_decay1 : gradient of H with respect to DECAY1.
% RETURN grad_D_decay2 : gradient of H with respect to DECAY1.
% RETURN grad_L : gradient of H with respect to length scale of
% latent process.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if delta_i~=delta_j
    [h, dh_dD_i, dh_dD_j, dh_dsigma] = diag(simComputeH(t, t, D_i, D_j, delta_i, delta_j, sigma));
else
    t = t - delta_i;

    halfSigmaD_i = 0.5*sigma*D_i;
    [lnPart1, sign1] = lnDiffErfs(halfSigmaD_i + t/sigma, halfSigmaD_i);
    
    [lnPart2, sign2] = lnDiffErfs(halfSigmaD_i, halfSigmaD_i - t/sigma);

    h = sign1.* exp(halfSigmaD_i*halfSigmaD_i  + lnPart1 - log(D_i + D_j)) ...
        - sign2.*exp(halfSigmaD_i*halfSigmaD_i - (D_i + D_j)*t + lnPart2 ...
        - log(D_i + D_j));
    
    sigma2 = sigma*sigma;

    if nargout > 1
        
        dh_dD_i = (0.5*D_i*sigma2*(D_i + D_j)-1)*h ... 
            + t.*sign2.*exp(halfSigmaD_i*halfSigmaD_i-(D_i+D_j)*t + lnPart2) ...
            + sigma/sqrt(pi)*(-1 + exp(-t.*t/sigma2-D_i*t) ...
            + exp(-t.*t/sigma2-D_j*t) - exp(-(D_i + D_j)*t));
        
        dh_dD_i = real(dh_dD_i/(D_i+D_j));
        
        
        if nargout > 2
        
            dh_dD_j = t.*sign2.*exp(halfSigmaD_i*halfSigmaD_i-(D_i + D_j)*t+lnPart2)-h;
            dh_dD_j = real(dh_dD_j/(D_i + D_j));

            if nargout > 3
                dh_dsigma = 0.5*D_i*D_i*sigma*h + 2/(sqrt(pi)*(D_i+D_j))*((-D_i/2) ...
                + (-t/sigma2+D_i/2).*exp(-t.*t/sigma2 - D_i*t) ...
                - (-t/sigma2-D_i/2).*exp(-t.*t/sigma2 - D_j*t) ...
                - D_i/2*exp(-(D_i+D_j)*t));
            end
        end
    end
end
