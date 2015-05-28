function [h, dh_dD_i, dh_dD_j, dh_dsigma] = simComputeHStat(t1, t2, D_i, D_j, delta_i, delta_j, sigma)

% SIMCOMPUTEHSTAT Helper function for computing part of the stationary version
%
%	Description:
%	of the SIM kernel.
%
%	H = SIMCOMPUTEHSTAT(T1, T2, DECAY1, DECAY2, L) computes a portion of
%	the stationary version of the SIM kernel.
%	 Returns:
%	  H - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1).
%	  DECAY1 - Decay rate for first system.
%	  DECAY2 - Decay rate for second system.
%	  L - length scale of latent process.
%	DESC computes a portion of the SIM kernel and gradients with
%	respect to various parameters.
%	ARG t1 : first time input (number of time points x 1).
%	ARG t2 : second time input (number of time points x 1).
%	ARG decay1 : Decay rate for first system.
%	ARG decay2 : Decay rate for second system.
%	ARG l : length scale of latent process.
%	RETURN h : result of this subcomponent of the kernel for the given values.
%	RETURN grad_D_decay1 : gradient of H with respect to DECAY1.
%	RETURN grad_D_decay2 : gradient of H with respect to DECAY1.
%	RETURN grad_L : gradient of H with respect to length scale of
%	latent process.
%	
%	
%	
%	
%
%	See also
%	SIMCOMPUTEH, SIMKERNPARAMINIT, LNDIFFERFS


%	Copyright (c) 2006 Neil D. Lawrence
%	Copyright (c) 2009 David Luengo


%	With modifications by Antti Honkela 2007, 2008


if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1 = t1 - delta_i;
t2 = t2 - delta_j;
diffT = (t1(:, ones(1, dim2)) - t2(:, ones(1, dim1))');
invSigmaDiffT = 1/sigma*diffT;
halfSigmaD_i = 0.5*sigma*D_i;
sigma2 = sigma*sigma;

lnPart1 = lnDiffErfs(inf, halfSigmaD_i - invSigmaDiffT);
h = exp(halfSigmaD_i*halfSigmaD_i - D_i*diffT + lnPart1 - log(D_i + D_j));

if nargout > 1
  % Gradient w.r.t. the first decay
  dh_dD_i = real((-1/(D_i + D_j) + sigma*halfSigmaD_i - diffT) .* h ...
      - sigma/((D_i+D_j)*sqrt(pi)) * exp(-invSigmaDiffT.*invSigmaDiffT));
  if nargout > 2
    % Gradient w.r.t. the second decay
    dh_dD_j = real(-1/(D_i + D_j) * h);   
    if nargout > 3
      % Gradient w.r.t. the length scale
      dh_dsigma = real(halfSigmaD_i * D_i * h ...
          - 2/(sqrt(pi)*(D_i+D_j)) * (0.5*D_i + invSigmaDiffT/sigma) ...
            .* exp(-invSigmaDiffT.*invSigmaDiffT));
    end
  end
end
