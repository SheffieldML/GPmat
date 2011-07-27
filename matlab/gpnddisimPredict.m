function [priormeans,posteriormeans,covmatrix,rbfposteriormeans,rbfcovmatrix] = gpnddisimPredict(model,predtimes,predict_rna);

% GPASIMPREDICT Compute predictions (means and a covariance matrix)
% of POL2 and RNA values for the GPASIM model.
%
% FORMAT
%---------------------------------
% DESC computes predictions for the asynchronous Gaussian
% process single input motif model.
%
% ARG model : the model for which the gradient is computed.
%
% ARG pol2times : the time points where predictions for POL2 are needed
%
% ARG rnatime : the time points where predictions for RNA are needed
%
% RETURN means : the predicted mean values, first the POL2
% predictions and then the RNA predictions.
%
% RETURN covmatrix : the covariance matrix between the
% predictions; for example, the diagonal values are the variances
% of each prediction.
%---------------------------------
%
% SEEALSO : gpasimCreate
%
% COPYRIGHT : Jaakko Peltonen, 2011

% GPASIMPREDICT


numGenes=model.numGenes;
if predict_rna==0,
  numGenes=0;
end;


% compute prior means 
pol2priormeans=ones(size(predtimes,1),1)*model.simMean;
%model.pol2mean;

if numGenes>0,
  % Mean for the mRNA is nonconstant over time and depends on the
  % B,D,S parameters and on the POL2 mean
  Bj=model.B(1);
  Dj=model.D(1);
  Sj=model.S(1);
  if model.use_disimstartmean==1,
    disimStartMean=model.disimStartMean(1);
  end;  
end;
  
if numGenes>0,
  % compute the RNA mean curve
%  if model.use_disimstartmean==1,
%    rnapriormeans=(Bj+model.simMean*Sj)/Dj+(disimStartMean-(Bj+model.simMean*Sj)/Dj)*exp(Dj*(-predtimes));
%  else    
%    rnapriormeans=(Bj+model.simMean*Sj)/Dj*ones(size(predtimes));  
%  end;

  rnapriormeans=[];
  tempind1=1;
  for k=1:numGenes,
    nt=length(predtimes);
    rnapriormeans=[rnapriormeans;nan*ones(nt,1)];
    if (model.use_disimstartmean==1),
      tempt=predtimes-model.delay(k);
      I=find(tempt<0);
      tempt(I)=0;
      rnapriormeans(tempind1:tempind1+nt-1)=...
          model.disimStartMean(k)*exp(model.D(k)*(-predtimes)) ...
          +(model.B(k)/model.D(k))*(1-exp(-model.D(k)*predtimes)) ...
          +(model.simMean*model.S(k)/model.D(k))*(1-exp(-model.D(k)*tempt));
    else
      tempt=predtimes-model.delay(k);
      I=find(tempt<0);
      tempt(I)=0;
      rnapriormeans(tempind1:tempind1+nt-1)=...
          ((model.B(k)+model.simMean*model.S(k))/model.D(k))*exp(model.D(k)*(-predtimes))...
          +((model.B(k)+model.simMean*model.S(k))/model.D(k))*(1-exp(-model.D(k)*tempt));
    end;
    tempind1=tempind1+nt;
  end;
  size(rnapriormeans)
  size(pol2priormeans)
end;

  
if 1,
K_new=kernCompute(model.kern, predtimes);
K_new_old=kernCompute(model.kern, predtimes, model.t);
K_old=model.K;
K_old_new=K_new_old';
end;

if 0,
K_old_ndsim=ndsimKernCompute(model.kern.comp{1}.comp{1},model.t);
K_old_nddisim=nddisimKernCompute(model.kern.comp{1}.comp{2},model.t);
K_old_nddisimXndsim=nddisimXndsimKernCompute(model.kern.comp{1}.comp{2},model.kern.comp{1}.comp{1},model.t);
K_old=[K_old_ndsim K_old_nddisimXndsim';K_old_nddisimXndsim K_old_nddisim];
K_old=real(K_old);

K_new_ndsim=ndsimKernCompute(model.kern.comp{1}.comp{1},predtimes);
K_new_nddisim=nddisimKernCompute(model.kern.comp{1}.comp{2},predtimes);
K_new_nddisimXndsim=nddisimXndsimKernCompute(model.kern.comp{1}.comp{2},model.kern.comp{1}.comp{1},predtimes);
K_new=[K_new_ndsim K_new_nddisimXndsim';K_new_nddisimXndsim K_new_nddisim];
K_new=real(K_new);

K_new_old_ndsim=ndsimKernCompute(model.kern.comp{1}.comp{1},predtimes,model.t);
K_new_old_nddisim=nddisimKernCompute(model.kern.comp{1}.comp{2},predtimes,model.t);
K_new_old_ndsimXnddisim=nddisimXndsimKernCompute(model.kern.comp{1}.comp{2},model.kern.comp{1}.comp{1},model.t,predtimes)';
K_new_old_nddisimXndsim=nddisimXndsimKernCompute(model.kern.comp{1}.comp{2},model.kern.comp{1}.comp{1},predtimes,model.t);
K_new_old=[K_new_old_ndsim K_new_old_ndsimXnddisim;K_new_old_nddisimXndsim K_new_old_nddisim];
K_new_old=real(K_new_old);

K_old_new=K_new_old';

noisekern_new=kernCompute(model.kern.comp{2}, predtimes);
noisekern_old=kernCompute(model.kern.comp{2}, model.t);
K_old=K_old+noisekern_old;
K_new=K_new+noisekern_new;
end;

tempm=model.m;
if numGenes==0,
  ot=length(model.t);
  nt=length(predtimes);
  K_old=K_old(1:ot,1:ot);
  K_new_old=K_new_old(1:nt,1:ot);
  K_old_new=K_old_new(1:ot,1:nt);
  K_new=K_new(1:nt,1:nt);
  tempm=tempm(1:ot);
end;


if numGenes>0,
  priormeans=[pol2priormeans;rnapriormeans];
else
  priormeans=pol2priormeans;
end;

%predict_rna
%K_old
%K_new_old
%pause
%figure; imagesc(K_old);
%pause

posteriormeans=priormeans+K_new_old*(K_old\tempm);
%covmatrix=K_new-K_new_old*inv(K_old)*K_old_new;

covmatrix=K_new-K_new_old*(K_old\K_old_new);
if(min(diag(covmatrix))<0),

  ot=length(model.t);
  nt=length(predtimes);

  if (numGenes>0),
    okentries=[(2:ot) (ot+2:2*ot)];
  else
    okentries=[2:ot];
  end;
  
  K_old=K_old(okentries,okentries);               
  K_new_old=K_new_old(:,okentries);
  K_old_new=K_old_new(okentries,:);
  covmatrix=K_new-K_new_old*(K_old\K_old_new);
end;



if 0,
  % Compute predictions for the driving RBF-kernel GP; only valid for
  % the SIM-DISIM model.


  % extract RBF transforms settings from DISIM part of the SIM-DISIM kernel
  disimkern=model.kern.comp{1}.comp{2};
  if isfield(disimkern, 'options') ...
        && isfield(disimkern.options, 'isNegativeS') ...
        && kern.options.isNegativeS,
    error('DISIM kern uses negative variance, cannot decide transform settings for RBF-SIM-DISIM prediction');
  else
    rbftransformsettings{1}=disimkern.transforms(1).transformsettings; % setting for inverse width
    rbftransformsettings{2}=[0 10]; % setting for RBF variance, not
                                    % used in GPNDDISIM model
  end
  

  % create a RBF kernel to match the settings in the NDSIM-NDDISIM model 
  rbfkern=kernCreate(predtimes,'rbf');
  rbfkern=kernExpandParamTransformSettings(rbfkern,rbftransformsettings);
  %[rbfpars,rbfnams]=kernExtractParam(rbfkern);
  [modelpars,modelnams]=gpnddisimExtractParam(model);
 % modelpars
 % modelnams
 % pause
  fprintf(1,'Calling gpdisimExpandParam\n');
  model=gpnddisimExpandParam(model,modelpars);
  fprintf(1,'Calling gpdisimExpandParam done\n');

  
  rbfpar_inversewidth=modelpars(1);  % assumes that inverse width is parameter 1 in the NDSIM-NDDISIM model
  rbfpar_variance=sigmoidabTransform(1, 'xtoa', rbftransformsettings{2});
  rbfpars(1)=rbfpar_inversewidth;
  rbfpars(2)=rbfpar_variance;
  rbfkern=kernExpandParam(rbfkern,rbfpars);
%  rbfkern
%  rbftransformsettings

  % create fake SIM and DISIM kernels corresponding to the NDSIM
  % and NDDISIM kernels... this is just a quick hack
  

  K_rbf=kernCompute(rbfkern,predtimes);
  K_rbf_sim=simXrbfKernCompute(model.kern.comp{1}.comp{1},rbfkern,model.t,predtimes)';
  K_rbf_disim=disimXrbfKernCompute(model.kern.comp{1}.comp{2},rbfkern,model.t,predtimes)';
  K_rbf_old=[K_rbf_sim K_rbf_disim];
  K_old_rbf=K_rbf_old';
  %size(K_rbf_old)
  %size(K_old)
  %size(model.m)
  
  K_old=model.K;
  rbfposteriormeans=0+K_rbf_old*(K_old\model.m);
  rbfcovmatrix=K_rbf-K_rbf_old*inv(K_old)*K_old_rbf;

  
%  rbfposteriormeans=real(rbfposteriormeans);
%  rbfcovmatrix=real(rbfcovmatrix);
else
  rbfposteriormeans=[];
  rbfcovmatrix=[];
end;