numSamp = 100;
t = linspace(100, 120, numSamp)';
kern = kernCreate(t, {'multi', 'rbf', 'lfm', 'lfm', 'lfm'});

kern.comp{2}.mass = 2;
kern.comp{2}.spring = 2;
kern.comp{2}.damper = 0.5;
kern.comp{3}.damper = 4;
kern.comp{4}.damper = 2.01;

for i = 2:length(kern.comp)
  zeta(i) = kern.comp{i}.damper/(2*sqrt(kern.comp{i}.spring* ...
                                      kern.comp{i}.mass));
  omega0(i) = sqrt(kern.comp{i}.spring/kern.comp{i}.mass);
end
params = kernExtractParam(kern);
kern = kernExpandParam(kern, params);

K = kernCompute(kern, t);

imagesc(real(K))

colorbar
numSamp2 = 4;
y =gsamp(zeros(size(K), 1), 0.5*(real(K) + real(K')), numSamp2)';
t = t - 100;
for i = 1:numSamp2
  figure
  x = reshape(y(:, i), numSamp, length(kern.comp));
  a = plot(t, x(:, 1), 'c-');
  hold on
  a= [a plot(t, x(:, 2), 'r-')];
  a = [a plot(t, x(:, 3), 'g-')];
  a = [a plot(t, x(:, 4), 'b-')];
  set(a, 'linewidth', 2)
  set(gca, 'fontname' 'arial');
  set(gca, 'fontsize', 14);
end
  
