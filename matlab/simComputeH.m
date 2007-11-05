function [h, dh_dD_i, dh_dD_j, dh_dsigma] = simComputeH(t1, t2, D_i, D_j, delta_i, delta_j, sigma)

% SIMCOMPUTEH Helper function for comptuing part of the SIM kernel.
% FORMAT
% DESC computes a portion of the SIM kernel.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG decay1 : Decay rate for first system.
% ARG decay2 : Decay rate for second system.
% ARG l : length scale of latent process.
% RETURN h : result of this subcomponent of the kernel for the given values.
%
% DESC computes a portion of the SIM kernel and gradients with
% respect to various parameters.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG decay1 : Decay rate for first system.
% ARG decay2 : Decay rate for second system.
% ARG l : length scale of latent process.
% RETURN h : result of this subcomponent of the kernel for the given values.
% RETURN grad_D_decay1 : gradient of H with respect to DECAY1.
% RETURN grad_D_decay2 : gradient of H with respect to DECAY1.
% RETURN grad_L : gradient of H with respect to length scale of
% latent process.
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Antti Honkela, 2007
%
% SEEALSO : simKernParamInit

% KERN

if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1 = t1 - delta_i;
t2 = t2 - delta_j;
t1Mat = repmat(t1, [1 dim2]);
t2Mat = repmat(t2', [dim1 1]);
diffT = (t1Mat - t2Mat);
invSigmaDiffT = 1/sigma*diffT;
halfSigmaD_i = 0.5*sigma*D_i;

%/~
%ind = find(t1Mat~=0);
%h = zeros(size(t1Mat));
%lnPart1 = zeros(size(h));
%lnPart2 = zeros(size(h));
%~/

lnPart1 = zeros(size(t1Mat));
warnState = warning('query', 'MATLAB:log:logOfZero');
warning('off', 'MATLAB:log:logOfZero');
I = (real(halfSigmaD_i - invSigmaDiffT) > 0);
J = ~I;
lnPart1(I) = log(  erfc( real(halfSigmaD_i - invSigmaDiffT(I))) ...
		 - erfc( real(halfSigmaD_i + t2Mat(I)/sigma)));
lnPart1(J) = log(- erfc(-real(halfSigmaD_i - invSigmaDiffT(J))) ...
		 + erfc(-real(halfSigmaD_i + t2Mat(J)/sigma)));
if real(halfSigmaD_i) > 0,
  lnPart2 = log(  erfc( real(halfSigmaD_i - t1Mat/sigma)) ...
		- erfc( real(halfSigmaD_i)));
else
  lnPart2 = log(- erfc(-real(halfSigmaD_i - t1Mat/sigma)) ...
		+ erfc(-real(halfSigmaD_i)));
end
warning(warnState.state, 'MATLAB:log:logOfZero');

h = exp(halfSigmaD_i*halfSigmaD_i ...
        -D_i*diffT+lnPart1 ...
        -log(D_i + D_j)) ...
    - exp(halfSigmaD_i*halfSigmaD_i ...
          -D_i*t1Mat-D_j*t2Mat ...
          +lnPart2 ...
          -log(D_i + D_j));


if nargout > 1
  sigma2 = sigma*sigma;
  dh_dD_i = (0.5*D_i*sigma2*(D_i + D_j)-1)*h ...
            + (-diffT.*exp(halfSigmaD_i*halfSigmaD_i-D_i*diffT+lnPart1) ...
               +t1Mat.*exp(halfSigmaD_i*halfSigmaD_i-D_i*t1Mat - D_j*t2Mat+lnPart2)) ...
            +sigma/sqrt(pi)*(-exp(-diffT.*diffT/sigma2)...
                             +exp(-t2Mat.*t2Mat/sigma2-D_i*t1Mat) ...
                             +exp(-t1Mat.*t1Mat/sigma2-D_j*t2Mat) ...
                             -exp(-(D_i*t1Mat + D_j*t2Mat)));
  dh_dD_i = real(dh_dD_i/(D_i+D_j));
  if nargout > 2
    dh_dD_j = t2Mat.*exp(halfSigmaD_i*halfSigmaD_i-(D_i*t1Mat + D_j*t2Mat)+lnPart2)-h;
    dh_dD_j = real(dh_dD_j/(D_i + D_j));
    
    if nargout > 3
      dh_dsigma = 0.5*D_i*D_i*sigma*h ...
          + 2/(sqrt(pi)*(D_i+D_j))*((-diffT/sigma2-D_i/ ...
                                     2).*exp(-diffT.* ...
                                 diffT/sigma2) ...
                                    + (-t2Mat/sigma2+D_i/2) ...
                                    .*exp(-t2Mat.*t2Mat/sigma2 ...
                              -D_i*t1Mat) ...
                                    - (-t1Mat/sigma2-D_i/2) ...
                                    .*exp(-t1Mat.*t1Mat/sigma2-D_j*t2Mat) ...
                                    - D_i/2*exp(-(D_i*t1Mat+D_j*t2Mat)));
    end
  end
end
