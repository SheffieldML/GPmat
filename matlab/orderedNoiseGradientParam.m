function g = orderedNoiseGradientParam(noise, mu, varsigma, y)

% ORDEREDNOISEGRADIENTPARAM Gradient of the ordered categorical noise model's parameters.

% IVM

D = size(y, 2);
c = 1./sqrt(noise.variance + varsigma);
gnoise.bias = zeros(1, D);
gnoise.widths = zeros(noise.C-2, 1);
for j = 1:D
  % Do lowest category first
  index = find(y(:, j)==0);
  if ~isempty(index)
    mu(index, j) = mu(index, j) + noise.bias(j) ;
    mu(index, j) = mu(index, j).*c(index, j);
    gnoise.bias(j) = gnoise.bias(j) - sum(c(index, j).*gradLogCumGaussian(-mu(index, j)));
  end

  % Intermediate categories
  index = find(y(:, j)>0 & y(:, j) <noise.C-1);
  if ~isempty(index)
    for i = index'
      mu(i, j) = mu(i, j) + noise.bias(j) - sum(noise.widths(1:y(i, j)-1));
    end
    u = mu(index, j).*c(index, j);
    uprime = (mu(index, j) - noise.widths(y(index, j))).*c(index, j);
    denom = (cumGaussian(u) - cumGaussian(uprime))+eps;
    gnoise.bias(j) = gnoise.bias(j) ...
        + sum(c(index, j) ...
              .*(ngaussian(u) - ngaussian(uprime))...
              ./denom);
    for cat = 1:noise.C-2
      
      subIndex = find(y(index, j) == cat);
      if ~isempty(subIndex)
        addpart = sum(c(index(subIndex), j)...
                      .*gaussOverDiffCumGaussian(u(subIndex), ...
                                                 uprime(subIndex), 2));
        gnoise.widths(1:cat) = gnoise.widths(1:cat) ...
            + repmat(addpart, cat, 1);
        if(cat > 1)
          addpart = sum(c(index(subIndex), j)...
                        .*gaussOverDiffCumGaussian(u(subIndex), ...
                                                   uprime(subIndex), 1));
          gnoise.widths(1:cat-1) = gnoise.widths(1:cat-1) ...
              - repmat(addpart, cat-1, 1);
        end
      end
    end
  end
  
  % Highest category
  index = find(y(:, j) == noise.C-1);
  if ~isempty(index)
    for i = index'
      mu(i, j) = mu(i, j) + noise.bias(j) - sum(noise.widths(1:y(i, j)-1));
    end
    mu(index, j) = mu(index, j).*c(index, j);
    addpart = sum(c(index, j).*gradLogCumGaussian(mu(index, j)));
    gnoise.bias(j) = gnoise.bias(j) + addpart;
    if length(noise.widths > 0)
      gnoise.widths = gnoise.widths ...
          - repmat(addpart, noise.C-2, 1);
    end
  end
end
if length(noise.widths>0)
  gnoise.widths = gnoise.widths.*gradFactLinearBound(noise.widths)';
  g = [gnoise.bias gnoise.widths(:)'];
else
  g = gnoise.bias;
end
