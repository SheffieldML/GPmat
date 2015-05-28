function K = nddisimXndsimKernCompute(disimKern, simKern, t1, t2)

% NDDISIMXNDSIMKERNCOMPUTE Compute a cross kernel between DISIM and SIM kernels with no decay in the SIM part.
% FORMAT
% DESC computes cross kernel terms between DISIM and SIM kernels for
% the multiple output kernel. 
% ARG disimKern : the kernel structure associated with the DISIM
% kernel.
% ARG simKern : the kernel structure associated with the SIM
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between DISIM and SIM kernels for
% the multiple output kernel. 
% ARG disimKern : the kernel structure associated with the DISIM
% kernel.
% ARG simKern : the kernel structure associated with the SIM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, disimKernParamInit, simKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007-2009
%
% COPYRIGHT : Jaakko Peltonen, 2011

% KERN


%fprintf(1,'nddisimXndsimKernCompute step1\n');
%t1
%t2


if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

if (isempty(t1)) || (isempty(t2)),
  K = zeros(size(t1,1),size(t2,1));
  return;
end;


if isfield(disimKern,'delay'),
  t1=t1-disimKern.delay;
  
  % crude way to handle times below zero, just truncate to zero
  % since the cross-kernel value at t1=0 is zero, which is the same
  % as cross-kernel values for negative t.
  I=find(t1<0);
  t1(I)=0;
end;



if disimKern.inverseWidth ~= simKern.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end
%if disimKern.di_decay ~= simKern.decay
%  error('Kernels cannot be cross combined if they have different driving input decays.');
%end
if disimKern.di_variance ~= simKern.variance
  error('Kernels cannot be cross combined if they have different driving input variances.');
end

dim1 = size(t1, 1);
dim2 = size(t2, 1);
%size(t1)
%t1
%size(t2)
%t2
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);

l = sqrt(2/disimKern.inverseWidth);
D_i = disimKern.decay;
%delta = disimKern.di_decay;

sigma = sqrt(2/simKern.inverseWidth);
if isfield(simKern, 'isNegativeS') && (simKern.isNegativeS == true)
    variancemultiplier = (simKern.sensitivity*simKern.sensitivity);
else
    variancemultiplier = simKern.variance;
end
variancemultiplier=variancemultiplier*sqrt(disimKern.variance);









k1a=-sqrt(pi)*variancemultiplier*sigma/(2*D_i*D_i)*...
    (exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat+lnDiffErfs(D_i*sigma/2+t2Mat/sigma,D_i*sigma/2+t2Mat/sigma-t1Mat/sigma))...
     +exp((D_i*sigma/2)^2-D_i*t1Mat+lnDiffErfs(D_i*sigma/2-t1Mat/sigma,D_i*sigma/2)));
k2a=-variancemultiplier*sigma*sigma/(2*D_i)*...
    (exp(-(diffT/sigma).^2)-exp(-(t1Mat/sigma).^2)+1-exp(-(t2Mat/sigma).^2));
k3a=-sqrt(pi)*variancemultiplier*sigma/(2*D_i)*...
    ((diffT-1/D_i).*erf(diffT/sigma)-(t1Mat-1/D_i).*erf(t1Mat/sigma)-(t2Mat+exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma));

%
%k1b
%k2b
%pause
K=real(k1a+k2a+k3a);
