function [dlnZ_dmu, dlnZ_dvs] = orderedNoiseGradVals(noise, mu, varsigma, y)

% ORDEREDNOISEGRADVALS Gradient of ORDERED noise log Z with respect to input mean and variance.
%
%	Description:
%
%	[DLNZ_DMU, DLNZ_DVS] = ORDEREDNOISEGRADVALS(NOISE, MU, VARSIGMA, Y)
%	computes the gradient of the ordered categorical noise with respect
%	to the input mean and the input variance.
%	 Returns:
%	  DLNZ_DMU - the gradient of log Z with respect to the input mean.
%	  DLNZ_DVS - the gradient of log Z with respect to the input
%	   variance.
%	 Arguments:
%	  NOISE - noise structure for which gradients are being computed.
%	  MU - mean input locations with respect to which gradients are
%	   being computed.
%	  VARSIGMA - variance input locations with respect to which
%	   gradients are being computed.
%	  Y - noise model output observed values associated with the given
%	   points.
%	
%
%	See also
%	% SEEALSO ORDEREDNOISEPARAMINIT, ORDEREDNOISEGRADIENTPARAM, NOISEGRADVALS, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



D = size(y, 2);
c = 1./sqrt(noise.variance + varsigma);

% Missing values are left untouched at zero.
dlnZ_dmu = zeros(size(c));
dlnZ_dvs = zeros(size(c));
for j = 1:D
  % Do lowest category first
  index = find(y(:, j)==0);
  if ~isempty(index)
    mu(index, j) = mu(index, j) + noise.bias(j) ;
    mu(index, j) = mu(index, j).*c(index, j);
    dlnZ_dmu(index, j) = -c(index, j).*gradLogCumGaussian(-mu(index, j));
    dlnZ_dvs(index, j) = -.5*dlnZ_dmu(index, j).*c(index, j).*mu(index, j);
  end

    
  % Intermediate categories
  index = find(y(:, j)>0 & y(:, j) <noise.C-1);
  if ~isempty(index)
    for i = index'
      mu(i, j) = mu(i, j) + noise.bias(j) - sum(noise.widths(1:y(i, j)-1));
    end
    u = mu(index, j).*c(index, j);
    uprime = (mu(index, j)-noise.widths(y(index, j))).*c(index, j);
    B1 = gaussOverDiffCumGaussian(u, uprime, 1);
    B2 = gaussOverDiffCumGaussian(u, uprime, 2);
    dlnZ_dmu(index, j) = c(index, j).*(B1-B2);
    
    dlnZ_dvs(index, j) = -.5*c(index, j).*c(index, j).*(u.*B1 - uprime.*B2);
  end
  
  % Highest category
  index = find(y(:, j) == noise.C-1);
  if ~isempty(index)
    for i = index'
      mu(i, j) = mu(i, j) + noise.bias(j) - sum(noise.widths(1:y(i, j)-1));
    end
    mu(index, j) = mu(index, j).*c(index, j);
    dlnZ_dmu(index, j) = c(index, j).*gradLogCumGaussian(mu(index, j));
    dlnZ_dvs(index, j) = -.5*dlnZ_dmu(index, j).*c(index, j).*mu(index, j);
    
  end
end
