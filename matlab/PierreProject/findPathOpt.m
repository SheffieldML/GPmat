function Xout = findPathOpt(model, xvec, SoFC, Yfinal ,sizePath, skel, display, iters)

if nargin < 5
  iters = 200;
  if nargin < 4
    display = true;
  end
end

options = optOptions;
if display
  options(1) = 1;
  options(9) = 1;
end
options(14) = iters;
options(9) = 1;

if isfield(model, 'optimiser')
  optim = str2func(model.optimiser);
else
  optim = str2func('scg');
end

Xout = optim('optimSequence', xvec,  options, 'optimGradient', model, SoFC, Yfinal, sizePath, skel);          
 
