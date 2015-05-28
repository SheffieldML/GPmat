function k = nddisimKernCompute(kern, t1, t2)

% DISIMKERNCOMPUTE Compute the DISIM kernel given the parameters
% and X, with SIM-level decay set to zero.
% FORMAT
% DESC computes the kernel parameters for the single input motif
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t1 : the input matrix associated with the rows of the kernel.
% ARG t2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the single input motif
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : disimKernParamInit, kernCompute, kernCreate, disimKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006; Antti Honkela, 2007
%
% COPYRIGHT : Jaakko Peltonen, 2011

% KERN

if nargin < 3
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

if isfield(kern,'delay'),
  t1=t1-kern.delay;
  t2=t2-kern.delay;
  
  % crude way to handle times below zero, just truncate to zero
  % since the kernel value at t=0 is zero, which is the same as
  % kernel values for negative t.
  I=find(t1<0);
  t1(I)=0;
  I=find(t2<0);
  t2(I)=0;
end;






% l = sqrt(2/kern.inverseWidth);
% 
% h1sum=nddisimComputeHSum(t1,t2,kern.decay,kern.decay,l);
% h2sum=nddisimComputeHSumPart2(t1,t2,kern.decay,l);
% k = h1sum - h2sum;

% hp = disimComputeHPrime(t, t2, kern.di_decay, kern.decay, kern.decay, l);
% hp2 = disimComputeHPrime(t2, t, kern.di_decay, kern.decay, kern.decay, l);
% k = h1sum + hp + hp2';


if 1,
  dim1 = size(t1, 1);
  dim2 = size(t2, 1);
  t1 = t1;
  t2 = t2;
  t1Mat = t1(:, ones(1, dim2));
  t2Mat = t2(:, ones(1, dim1))';
  diffT = (t1Mat - t2Mat);
  
  D = kern.decay;
  l = sqrt(2/kern.inverseWidth);
  
  h1a=0.5*sqrt(pi)*l*kern.di_variance*kern.variance/(D^2) * ...
      (  (t1Mat - 1/D + exp(-D*t2Mat)/D).*erf(t1Mat/l) ...
	 +(t2Mat - 1/D + exp(-D*t1Mat)/D).*erf(t2Mat/l) ...
	 -(diffT).*erf(diffT/l) );
  h1b=0.5*(l^2)*kern.di_variance*kern.variance/(D^2) * ...
      ( exp(-(t1Mat/l).^2) + exp(-(t2Mat/l).^2) - exp(-(diffT/l).^2) - 1 );
  
  h2a=sqrt(pi)*l*kern.di_variance*kern.variance * ...
      ( 0.25*exp(-3*log(D) + (D*l/2)^2 + D*(diffT) + lnDiffErfs(D*l/2+t1Mat/l,D*l/2+diffT/l)) ...
      +0.25*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat - D*t1Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) ...
      -0.5*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) );
  h2b=sqrt(pi)*l*kern.di_variance*kern.variance * ...
      ( 0.25*exp(-3*log(D) + (D*l/2)^2 - D*(diffT) + lnDiffErfs(D*l/2+t2Mat/l,D*l/2-diffT/l)) ...
      +0.25*exp(-3*log(D) + (D*l/2)^2 - D*t1Mat - D*t2Mat + lnDiffErfs(D*l/2,D*l/2-t1Mat/l)) ...
      -0.5*exp(-3*log(D) + (D*l/2)^2 - D*t1Mat + lnDiffErfs(D*l/2,D*l/2-t1Mat/l)) );
  
  k=h1a + h1b - h2a - h2b;
end;


% k = 0.5*sqrt(pi)*l*k;
% k = kern.rbf_variance*kern.di_variance*kern.variance*k;
k = real(k);


% gaussianInitial currently unsupported
% if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
%   dim1 = size(t, 1);
%   dim2 = size(t2, 1);
%   t1Mat = t(:, ones(1, dim2));
%   t2Mat = t2(:, ones(1, dim1))';
% 
%   k = k + kern.initialVariance * kern.variance * ...
%       (exp(- kern.di_decay * t1Mat) - exp(- kern.decay * t1Mat)) ./ (kern.decay - kern.di_decay) .* ...
%       (exp(- kern.di_decay * t2Mat) - exp(- kern.decay * t2Mat)) ./ (kern.decay - kern.di_decay);
% end


if 0,
% Cantor/Sage version
  h1a(t1,t2,D,l)=1/(D^2) * (  (t1Mat - 1/D + exp(-D*t2Mat)/D).*erf(t1Mat/l) ...
		   +(t2Mat - 1/D + exp(-D*t1Mat)/D).*erf(t2Mat/l) ...
		   -(diffT).*erf(diffT/l) );
  h1b=l/sqrt(pi)/(D^2) * ( exp(-(t1Mat/l).^2) + exp(-(t2Mat/l).^2) - exp(-(diffT/l).^2) - 1 );
  
  h2a=0.5*exp(-3*log(D) + (D*l/2)^2 + D*(diffT) + lnDiffErfs(D*l/2+t1Mat/l,D*l/2+diffT/l)) ...
      +0.5*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat - D*t1Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) ...
      -exp(-3*log(D) + (D*l/2)^2 - D*t2Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l));
  h2b=0.5*exp(-3*log(D) + (D*l/2)^2 - D*(diffT) + lnDiffErfs(D*l/2+t2Mat/l,D*l/2-diffT/l)) ...
      +0.5*exp(-3*log(D) + (D*l/2)^2 - D*t1Mat - D*t2Mat + lnDiffErfs(D*l/2,D*l/2-t1Mat/l)) ...
      -exp(-3*log(D) + (D*l/2)^2 - D*t1Mat + lnDiffErfs(D*l/2,D*l/2-t1Mat/l));
end;
