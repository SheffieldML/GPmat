% DEMRIEMAN Demonstrate the numerical approximation to the integral.


thin = 2;
thick = 3;
numVals = 161;
D = 0.1;
S = 1;
randn('seed', 1e5)
rand('seed', 1e5)
t = linspace(-2, 14, numVals)';
intInd = 21:10:numVals;
kern = kernCreate(t, 'rbf');
kern.inverseWidth = 0.5;
K = kernCompute(kern, t);
f = gsamp(zeros(numVals, 1)', K, 1)';
g = exp(f);
e1 = exp(D*t);
ge1 = exp(f + D*t);
lin = plot(t, g, 'g');
set(lin, 'linewidth', thick)
hold on
% Print g
pause

lin = [lin plot(t, e1, 'g')];
set(lin, 'linewidth', thick)
% print g and e1
pause

lin2 = plot(t, ge1);
set(lin2, 'linewidth', thick)
set(lin, 'linewidth', thin)
% print ge1

set(lin, 'visible', 'off')
ylim = get(gca, 'ylim');
xlim = get(gca, 'xlim');

lin3 = line([0 0], ylim);

lin4 = line([t(1) t(1)], [0 ge1(2)]);
lin5 =[];
for i = 2:length(intInd)
  i1 = (intInd(i-1));
  i2 = (intInd(i));
  lin5 = [lin5  line(t(i2), ge1(i2))];
  set(lin5, 'marker', '.');
  set(lin5, 'markersize', 30);
  set(lin5, 'color', [1 0 0]); 
  % print new point
  
  set(lin4, 'linewidth', thin);
  lin4 = line([t(i2) t(i2)], [0 ge1(i2)]);
  lin4 = [lin4 line([t(i2) t(i1)], [ge1(i2) ge1(i2)])];
  lin4 = [lin4 line([t(i1) t(i1)], [ge1(i2) 0])];
  set(lin4, 'linewidth', thick)
  set(lin4, 'color', [1 0 0]);
end
