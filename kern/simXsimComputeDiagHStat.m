function [h, dh_dD_i, dh_dD_j, dh_dsigma] = simXsimComputeDiagHStat(t, D_i, D_j, delta_i, delta_j, sigma)

% SIMXSIMCOMPUTEDIAGHSTAT Helper function for computing part of the stationary version
% of the diagonal of the SIM kernel.
% FORMAT
% DESC computes a portion of the stationary version of the diagonal of he SIM kernel.
% ARG t : first time input (number of time points x 1).
% ARG decay1 : Decay rate for first system.
% ARG decay2 : Decay rate for second system.
% ARG l : length scale of latent process.
% RETURN h : result of this subcomponent of the kernel for the given values.
%
% DESC computes a portion of the SIM kernel and gradients with
% respect to various parameters.
% ARG t : first time input (number of time points x 1).
% ARG decay1 : Decay rate for first system.
% ARG decay2 : Decay rate for second system.
% ARG l : length scale of latent process.
% RETURN h : result of this subcomponent of the kernel for the given values.
% RETURN grad_D_decay1 : gradient of H with respect to DECAY1.
% RETURN grad_D_decay2 : gradient of H with respect to DECAY1.
% RETURN grad_L : gradient of H with respect to length scale of
% latent process.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010
%
% SEEALSO : simXsimComputeDiagH, simKernParamInit, lnDiffErfs

% KERN

if delta_i~=delta_j
    [h, dh_dD_i, dh_dD_j, dh_dsigma] = diag(simComputeHStat(t, t, D_i, D_j, delta_i, delta_j, sigma));
else
    t = t - delta_i;
    halfSigmaD_i = 0.5*sigma*D_i;
    lnPart1 = lnDiffErfs(inf, halfSigmaD_i);
    h = exp(halfSigmaD_i*halfSigmaD_i + lnPart1 - log(D_i + D_j));
    h = h*ones(size(t));
    if nargout > 1
        % Gradient w.r.t. the first decay
        dh_dD_i = real((-1/(D_i + D_j) + sigma*halfSigmaD_i) .* h ...
            - sigma/((D_i+D_j)*sqrt(pi)));
        if nargout > 2
            % Gradient w.r.t. the second decay
            dh_dD_j = real(-1/(D_i + D_j) * h);
            if nargout > 3
                % Gradient w.r.t. the length scale
                dh_dsigma = real(halfSigmaD_i * D_i * h - D_i/(sqrt(pi)*(D_i+D_j)));
            end
        end
    end
end
