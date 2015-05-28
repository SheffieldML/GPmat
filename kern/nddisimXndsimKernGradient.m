function [g1,g2] = nddisimXndsimKernGradient(disimKern, simKern, t1, t2, covGrad)

% NDDISIMXNDSIMKERNGRADIENT Compute a cross gradient between NDDISIM and NDSIM
% kernels with no decay in the SIM part.
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

if nargin < 5
  covGrad = t2;  
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  t1
  t2
  error('Input can only have one column');
end

dim1 = size(t1, 1);
dim2 = size(t2, 1);

if (isempty(t1)==1) || (isempty(t2)==1),
  % Complete gradient of NDDISIM parameters
  %g1=[dk_dInverseWidth dk_ddi_variance dk_dD dk_dvariance dk_ddelay];
  g1=[0 0 0 0 0];

  % gradients with respect to NDSIM kernel parameters
  g2 = [0 0];

  return;
end;





if isfield(disimKern,'delay'),
  delay=disimKern.delay;
  origt1=t1;
  origt1PosFlag=t1>delay;
  %size(origt1PosFlag)
  %dim2
  origt1PosFlag=origt1PosFlag(:, ones(1, dim2));
  t1=t1-delay;
  
  % crude way to handle times below zero, just truncate to zero
  % since the cross-kernel value at t1=0 is zero, which is the same
  % as cross-kernel values for negative t.
  I=find(t1<0);
  t1(I)=0;
end;



if disimKern.inverseWidth ~= simKern.inverseWidth
  error(sprintf('Kernels cannot be cross combined if they have different inverse widths (NDDISIM says %f, NDSIM says %f).',disimKern.inverseWidth,simKern.inverseWidth));
end
if disimKern.di_variance ~= simKern.variance
  error(sprintf('Kernels cannot be cross combined if they have different driving input variances (NDDISIM says %f, NDSIM says %f).',disimKern.di_variance,simKern.variance));
end

t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);

sigma = sqrt(2/disimKern.inverseWidth);
D_i = disimKern.decay;
variancemultiplier=disimKern.di_variance*sqrt(disimKern.variance);


[lndiff1,lndiffsigns1]=lnDiffErfs(D_i*sigma/2+t2Mat/sigma,D_i*sigma/2+t2Mat/sigma-t1Mat/sigma);
[lndiff2,lndiffsigns2]=lnDiffErfs(D_i*sigma/2-t1Mat/sigma,D_i*sigma/2);
% I=find(~isfinite(lndiff1)); lndiff1(I)=0; lndiffsigns1(I)=0;
% I=find(~isfinite(lndiff2)); lndiff2(I)=0; lndiffsigns2(I)=0;
%size(lndiff1)
%size(lndiff2)
%size(lndiffsigns1)
%size(lndiffsigns2)
%size(t1Mat)
%size(t2Mat)
%atemp1=lndiffsigns1.*exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat+lndiff1)
%atemp2=lndiffsigns2.*exp((D_i*sigma/2)^2-D_i*t1Mat+lndiff2)
k1a=-sqrt(pi)*sigma/(2*D_i*D_i)*...
    ( lndiffsigns1.*exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat+lndiff1)...
     +lndiffsigns2.*exp((D_i*sigma/2)^2-D_i*t1Mat+lndiff2) );
k2a=-sigma*sigma/(2*D_i)*...
    (exp(-(diffT/sigma).^2)-exp(-(t1Mat/sigma).^2)+1-exp(-(t2Mat/sigma).^2));
k3a=-sqrt(pi)*sigma/(2*D_i)*...
    ((diffT-1/D_i).*erf(diffT/sigma)-(t1Mat-1/D_i).*erf(t1Mat/sigma)-(t2Mat+exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma));

%
%k1b
%k2b
%pause
K=variancemultiplier*real(k1a+k2a+k3a);


% gradients with respect to NDDISIM kernel parameters:

%----------------------------------------
% Gradient with respect to inverse width
%----------------------------------------

[lndiff1,lndiffsigns1]=lnDiffErfs(D_i*sigma/2+t2Mat/sigma,D_i*sigma/2+(t2Mat-t1Mat)/sigma);
[lndiff2,lndiffsigns2]=lnDiffErfs(D_i*sigma/2-t1Mat/sigma,D_i*sigma/2);
% I=find(~isfinite(lndiff1)); lndiff1(I)=0; lndiffsigns1(I)=0;
% I=find(~isfinite(lndiff2)); lndiff2(I)=0; lndiffsigns2(I)=0;

dK_dsigma_part1a=...
    -(sqrt(pi)/2)*(1/(D_i*D_i) + (sigma^2)/2)*...
     ( lndiffsigns1.*exp((D_i*sigma/2)^2-D_i*(t1Mat-t2Mat)+lndiff1) + lndiffsigns2.*exp((D_i*sigma/2)^2-D_i*t1Mat+lndiff2) ) ...
    ...
    -1/D_i * ...
      (  (sigma/2-t2Mat/(sigma*D_i)).*exp((D_i*sigma/2)^2 - D_i*(t1Mat-t2Mat) -(D_i*sigma/2+t2Mat/sigma).^2) ...
       - (sigma/2+(t1Mat-t2Mat)/(sigma*D_i)).*exp((D_i*sigma/2)^2 - D_i*(t1Mat-t2Mat) -(D_i*sigma/2+(t2Mat-t1Mat)/sigma).^2) ) ...
    ...
    -1/D_i *  ...
     (   (sigma/2+t1Mat/(sigma*D_i)).*exp((D_i*sigma/2)^2 - D_i*t1Mat -(D_i*sigma/2-t1Mat/sigma).^2) ...
       - (sigma/2)*exp((D_i*sigma/2)^2 - D_i*t1Mat -(D_i*sigma/2)^2) ) ;

% dK_dsigma_part1a=...
%     -(sqrt(pi)/2)*(1/(D_i*D_i) + (sigma^2)/2)*...
%      lndiffsigns1.*exp((D_i*sigma/2)^2-D_i*(t1Mat-t2Mat)+lndiff1) ...
%     -1/D_i * exp((D_i*sigma/2)^2 - D_i*(t1Mat-t2Mat)) .* ...
%       (  (sigma/2-t2Mat/(sigma*D_i)).*exp(-(D_i*sigma/2+t2Mat/sigma).^2) ...
%        - (sigma/2+(t1Mat-t2Mat)/(sigma*D_i)).*exp(-(D_i*sigma/2+(t2Mat-t1Mat)/sigma).^2) ) ...
%     ...
%     -(sqrt(pi)/2)*(1/(D_i*D_i) + (sigma^2)/2)*...
%      lndiffsigns2.*exp((D_i*sigma/2)^2-D_i*t1Mat+lndiff2) ...
%     -1/D_i * exp((D_i*sigma/2)^2 - D_i*t1Mat).* ...
%      (   (sigma/2+t1Mat/(sigma*D_i)).*exp(-(D_i*sigma/2-t1Mat/sigma).^2) ...
%        - (sigma/2)*exp(-(D_i*sigma/2)^2) ) ;
  
    
% lndiff1
% lndiff2
% pause
    
    
dK_dsigma_part2a=...
    -sigma/D_i*...
    (exp(-(diffT/sigma).^2)-exp(-(t1Mat/sigma).^2)+1-exp(-(t2Mat/sigma).^2)) ...
    -1/(sigma*D_i)*...
    ((diffT.^2).*exp(-(diffT/sigma).^2)-(t1Mat.^2).*exp(-(t1Mat/sigma).^2)-(t2Mat.^2).*exp(-(t2Mat/sigma).^2));

[lndiff1,lndiffsigns1]=lnDiffErfs(t1Mat/sigma,diffT/sigma);
[lndiff2,lndiffsigns2]=lnDiffErfs(t2Mat/sigma,-diffT/sigma);
% I=find(~isfinite(lndiff1)); lndiff1(I)=0; lndiffsigns1(I)=0;
% I=find(~isfinite(lndiff2)); lndiff2(I)=0; lndiffsigns2(I)=0;
dK_dsigma_part3a=...
    -sqrt(pi)/(2*D_i)*...
    (-t1Mat.*lndiffsigns1.*exp(lndiff1) -t2Mat.*lndiffsigns2.*exp(lndiff2) ...
     +1/D_i*(erf(t1Mat/sigma)-erf(diffT/sigma)) - (exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma) ) ...
    -1/(D_i*sigma)*...
    ( (diffT-1/D_i).*(-diffT).*exp(-(diffT/sigma).^2) ...
     -(t1Mat-1/D_i).*(-t1Mat).*exp(-(t1Mat/sigma).^2) ...
     -(t2Mat+exp(-D_i*t1Mat)/D_i).*(-t2Mat).*exp(-(t2Mat/sigma).^2) );

% dK_dsigma_part1a
% dK_dsigma_part2a
% dK_dsigma_part3a
% pause


dK_dsigma = dK_dsigma_part1a + dK_dsigma_part2a + dK_dsigma_part3a;
dk_dInverseWidth = variancemultiplier*sum(sum(dK_dsigma.*covGrad))*(-disimKern.inverseWidth^(-1.5)/sqrt(2));


%----------------------------------------
% Gradient with respect to NDDISIM-level variance
%----------------------------------------
dK_dvariance = disimKern.di_variance*0.5/sqrt(disimKern.variance)*real(k1a+k2a+k3a);
dk_dvariance = sum(sum(dK_dvariance.*covGrad));


%----------------------------------------
% Gradient with respect to NDSIM-level variance
%----------------------------------------
dK_ddi_variance = sqrt(disimKern.variance)*real(k1a+k2a+k3a);
dk_ddi_variance = sum(sum(dK_ddi_variance.*covGrad));


%----------------------------------------
% Gradient with respect to delay
%----------------------------------------
[lndiff1,lndiffsigns1]=lnDiffErfs(D_i*sigma/2+t2Mat/sigma,D_i*sigma/2+t2Mat/sigma-t1Mat/sigma);
[lndiff2,lndiffsigns2]=lnDiffErfs(D_i*sigma/2-t1Mat/sigma,D_i*sigma/2);
% I=find(~isfinite(lndiff1)); lndiff1(I)=0; lndiffsigns1(I)=0;
% I=find(~isfinite(lndiff2)); lndiff2(I)=0; lndiffsigns2(I)=0;
dK_ddelay_part1a=...
    -sqrt(pi)*sigma/(2*D_i*D_i)*origt1PosFlag.*...
    ( lndiffsigns1.*exp(((D_i*sigma/2)^2)-D_i*(t1Mat-t2Mat)   +lndiff1)*D_i ...
      -2/(sigma*sqrt(pi))*exp((D_i*sigma/2)^2-D_i*(t1Mat-t2Mat) -(D_i*sigma/2+(t2Mat-t1Mat)/sigma).^2) ...
     +lndiffsigns2.*exp((D_i*sigma/2)^2-D_i*t1Mat+lndiff2)*D_i ...
     +2/(sigma*sqrt(pi))*exp((D_i*sigma/2)^2-D_i*t1Mat -(D_i*sigma/2-t1Mat/sigma).^2)  );

dK_ddelay_part2a=...
     -1/D_i*origt1PosFlag .* ( exp(-(diffT/sigma).^2).*diffT - exp(-(t1Mat/sigma).^2).*t1Mat );

[lndiff1,lndiffsigns1]=lnDiffErfs(t1Mat/sigma,diffT/sigma);
% I=find(~isfinite(lndiff1)); lndiff1(I)=0; lndiffsigns1(I)=0;
dK_ddelay_part3a=...
    -sqrt(pi)*sigma/(2*D_i)*origt1PosFlag.* ...
    ( lndiffsigns1.*exp(lndiff1) -exp(-D_i*t1Mat).*erf(t2Mat/sigma) ) ...
    -1/D_i*origt1PosFlag.* ...
    ( -(diffT-1/D_i).*exp(-(diffT/sigma).^2) +(t1Mat-1/D_i).*exp(-(t1Mat/sigma).^2) );
        
dK_ddelay = dK_ddelay_part1a + dK_ddelay_part2a + dK_ddelay_part3a;    
dk_ddelay = variancemultiplier*sum(sum(dK_ddelay.*covGrad));


%----------------------------------------
% Gradient with respect to decay
%----------------------------------------

[lndiff1,lndiffsigns1]=lnDiffErfs(D_i*sigma/2+t2Mat/sigma,D_i*sigma/2+t2Mat/sigma-t1Mat/sigma);
[lndiff2,lndiffsigns2]=lnDiffErfs(D_i*sigma/2-t1Mat/sigma,D_i*sigma/2);
% I=find(~isfinite(lndiff1)); lndiff1(I)=0; lndiffsigns1(I)=0;
% I=find(~isfinite(lndiff2)); lndiff2(I)=0; lndiffsigns2(I)=0;
dK_dD_part1a=...
   +2*sqrt(pi)*sigma/(2*(D_i^3))*...
    ( lndiffsigns1.*exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat+lndiff1)...
     +lndiffsigns2.*exp((D_i*sigma/2)^2-D_i*t1Mat+lndiff2) ) ...
   -sqrt(pi)*sigma/(2*D_i*D_i)*...
    ( lndiffsigns1.*exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat+lndiff1) .*(D_i*(sigma^2)/2 -t1Mat+t2Mat) ...
     +( exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat -(D_i*sigma/2+t2Mat/sigma).^2) ...
       -exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat -(D_i*sigma/2+t2Mat/sigma-t1Mat/sigma).^2) )*sigma/sqrt(pi) ...
     ...
     +lndiffsigns2.*exp((D_i*sigma/2)^2-D_i*t1Mat+lndiff2) .*(D_i*(sigma^2)/2 -t1Mat) ...
     + (  exp((D_i*sigma/2)^2-D_i*t1Mat -(D_i*sigma/2-t1Mat/sigma).^2) ...
        - exp((D_i*sigma/2)^2-D_i*t1Mat -(D_i*sigma/2).^2) )*sigma/sqrt(pi)   );
    
dK_dD_part2a=...
    +sigma*sigma/(2*(D_i^2))*...
    (exp(-(diffT/sigma).^2)-exp(-(t1Mat/sigma).^2)+1-exp(-(t2Mat/sigma).^2));
    
[lndiff1,lndiffsigns1]=lnDiffErfs(t1Mat/sigma,diffT/sigma);
% I=find(~isfinite(lndiff1)); lndiff1(I)=0; lndiffsigns1(I)=0;
dK_dD_part3a=...
    +sqrt(pi)*sigma/(2*(D_i^2))*...
    ( (diffT-1/D_i).*erf(diffT/sigma)-(t1Mat-1/D_i).*erf(t1Mat/sigma)-(t2Mat+exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma) ) ...
    -sqrt(pi)*sigma/(2*D_i)*...
    (-1/(D_i^2)*lndiffsigns1.*exp(lndiff1) ...
     +((t1Mat+1/D_i).*exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma) );
    
dK_dD = dK_dD_part1a + dK_dD_part2a + dK_dD_part3a;
dk_dD = variancemultiplier*sum(sum(dK_dD.*covGrad));


%----------------------------------------
% Complete gradient of NDDISIM parameters
%----------------------------------------
g1=[dk_dInverseWidth dk_ddi_variance dk_dD dk_dvariance dk_ddelay];



%----------------------------------------
% gradients with respect to NDSIM kernel parameters: note that
% gradients of shared parameters are placed in the
% NDDISIM-parameter gradient only. Since NDDISIM contains all
% parameters of NDSIM, the NDSIM gradient gets only zeros.
%----------------------------------------

g2 = [0 0];

