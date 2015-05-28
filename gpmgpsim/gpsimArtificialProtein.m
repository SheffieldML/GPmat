function f = gpsimArtificialProtein(t, alpha, mu, sigma);

% GPSIMARTIFICIALPROTEIN Generate an artifical protein sequence.

% SHEFFIELDML

f = zeros(size(t));
for i = 1:length(alpha)
  difft = (t - mu(i))/sigma(i);
  f = f + alpha(i)*exp(-difft.*difft);
end
