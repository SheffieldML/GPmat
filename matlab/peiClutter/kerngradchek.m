load step0.mat

% Likelihood gradient check for all parameters
model = model.comp{1};

eps = 1.0e-4;
x = model.mapt;
f = model.f;

step = 1.0;
arg1 = 1 + eps*step;
arg2 = 1 - eps*step;

[k, n2] = rbfKernCompute(model.kern, x);
K = k + eye(size(k, 1))*1e-6;
invK =  pdinv(K);
% g = rbfKernGradient(model.kern, x);
% covGrad = 0.5*f'*invK*g*invK*f;
covGrad1 = invK*(f*f')*invK;
gr = rbfKernGradient(model.kern, x, covGrad1);
covGrad2 = model.invK*(model.f*model.f')*model.invK;
gkern = kernGradient(model.kern, model.mapt, covGrad2);


model.kern.inverseWidth = arg1;
[k, n2] = rbfKernCompute(model.kern, x);
K = k + eye(size(k, 1))*1e-6;
invK =  pdinv(K);
ll1 = -.5*f'*invK*f;

model.kern.inverseWidth = arg2;
[k, n2] = rbfKernCompute(model.kern, x);
K = k + eye(size(k, 1))*1e-6;
invK =  pdinv(K);
ll2 = -.5*f'*invK*f;

diff = (ll1 - ll2)/(2*eps);
g = 0.5*gr(1);

fprintf('Checking gpsimMap gradient ...\n');
delta = g - diff;
fprintf(1, '   analytic   diffs     delta\n\n');
disp([g', diff', delta'])