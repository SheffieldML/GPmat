function [dlnZ_dmu, dlnZ_dvs] = orderedNoiseGradVals(noise, mu, varsigma, y)

% ORDEREDNOISEGRADVALS Gradient wrt mu and varsigma of log-likelihood for ordered categorical noise model.

% IVM1

D = size(y, 2);
c = 1./sqrt(varsigma);
dlnZ_dmu = zeros(size(c));
dlnZ_dvs = zeros(size(c));
for j = 1:D
  % Do lowest category first
  index = find(y(:, j)==0);
  if ~isempty(index)
    mu(index, j) = mu(index, j) + noise.bias(j) ;
    mu(index, j) = mu(index, j).*c(index, j);
    dlnZ_dmu(index, j) = -c(index, j).*ngaussian(mu(index, j))...
        ./(cumGaussian(-mu(index, j))+noise.eta/(1-noise.C*noise.eta));
    %/~
    % This is the stable way of computing Gaussian/erf
    % dlnZ_dmu(index, j) = -c(index, j)./...
    %    (sqrt(2*pi)...
    %      *(exp(1/2*mu(index, j).^2)...
    %       -.5*erfcx(-sqrt(2)/2*mu(index, j))));
    %~/
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
    dlnZ_dmu(index, j) = c(index, j).*(ngaussian(u) - ngaussian(uprime))./(cumGaussian(u) - cumGaussian(uprime)+noise.eta/(1-noise.C*noise.eta));
    dlnZ_dvs(index, j) = -.5*c(index, j).*c(index, j).*(u.*ngaussian(u) ...
                                    - uprime.*ngaussian(uprime))...
        ./(cumGaussian(u)...
           - cumGaussian(uprime)+noise.eta/(1-noise.C*noise.eta));
  end
  
  % Highest category
  index = find(y(:, j) == noise.C-1);
  if ~isempty(index)
    for i = index'
      mu(i, j) = mu(i, j) + noise.bias(j) - sum(noise.widths(1:y(i, j)-1));
    end
    mu(index, j) = mu(index, j).*c(index, j);
    dlnZ_dmu(index, j) = c(index, j).*ngaussian(mu(index, j))...
        ./(cumGaussian(mu(index, j)) + noise.eta/(1-noise.C*noise.eta));
    %/~
    %     dlnZ_dmu(index, j) = c(index, j)./...
    %         (sqrt(2*pi)...
    %          *(exp(1/2*mu(index, j).^2)...
    %            -.5*erfcx(sqrt(2)/2*mu(index, j))));
    %~/
    dlnZ_dvs(index, j) = -.5*dlnZ_dmu(index, j).*c(index, j).*mu(index, j);
    
  end
end
