function [priormeans,posteriormeans,covmatrix] = gpnddisimPredictRNAConditional(model,predtimes);

% GPASIMPREDICT Compute predictions (means and a covariance matrix)
% of RNA values for the GPASIM model, conditional on the existing
% POL values in the model but not on any existing RNA values.
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
% SEEALSO : gpnddisimCreate
%
% COPYRIGHT : Jaakko Peltonen, 2011

% GPASIMPREDICT


numGenes=model.numGenes;
if numGenes==0,
  error('gpnddisimPredictRNAConditional requires a model that has parameters for RNA modeling\n');
end;


% compute prior means 
%pol2priormeans=[1:size(predtimes,1)]'*0;
if iscell(predtimes)==0,
  pol2priormeans=ones(size(predtimes,1),1)*model.simMean;
else
  pol2priormeans=ones(size(predtimes{1},1),1)*model.simMean;
end;

% Mean for the mRNA is nonconstant over time and depends on the
% B,D,S parameters and on the POL2 mean
if numGenes>0,
  Bj=model.B(1);
  Dj=model.D(1);
  Sj=model.S(1);
  if model.use_disimstartmean==1,
    disimStartMean=model.disimStartMean(1);
  end;  
end;

% compute the RNA mean curve
if numGenes>0,
  rnapriormeans=[];
  tempind1=1;
  for k=1:numGenes,
    if iscell(predtimes)==0,
      nt=length(predtimes);
    else
      nt=length(predtimes{k+1});
    end;
    rnapriormeans=[rnapriormeans;nan*ones(nt,1)];
    if (model.use_disimstartmean==1),
      if iscell(predtimes)==0,
        tempt=predtimes;
      else
        tempt=predtimes{k+1};
      end;
      delayedt=tempt-model.delay(k);
      I=find(delayedt<0);
      delayedt(I)=0;
      if length(delayedt)>0,
        rnapriormeans(tempind1:tempind1+nt-1)=...
            model.disimStartMean(k)*exp(model.D(k)*(-tempt)) ...
            +(model.B(k)/model.D(k))*(1-exp(-model.D(k)*tempt)) ...
            +(model.simMean*model.S(k)/model.D(k))*(1-exp(-model.D(k)*delayedt));
      end;
    else
      if iscell(predtimes)==0,
        tempt=predtimes;
      else
        tempt=predtimes{k+1};
      end;
      delayedt=tempt-model.delay(k);
      I=find(delayedt<0);
      delayedt(I)=0;
      if length(delayedt)>0,
        rnapriormeans(tempind1:tempind1+nt-1)=...
            ((model.B(k)+model.simMean*model.S(k))/model.D(k))*exp(model.D(k)*(-tempt))...
            +((model.B(k)+model.simMean*model.S(k))/model.D(k))*(1-exp(-model.D(k)*delayedt));
      end;
    end;
    tempind1=tempind1+nt;
  end;
end;

%if model.use_disimstartmean==1,
%  rnapriormeans=(Bj+model.simMean*Sj)/Dj+(disimStartMean-(Bj+model.simMean*Sj)/Dj)*exp(Dj*(-predtimes));
%else    
%  rnapriormeans=(Bj+model.simMean*Sj)/Dj*ones(size(predtimes));
%end;


if 1,
  K_new=kernCompute(model.kern, predtimes);
  K_new_old=kernCompute(model.kern, predtimes, model.t);

  if iscell(predtimes)==0,
    K_new=K_new(length(predtimes)+1:2*length(predtimes),length(predtimes)+1:2*length(predtimes));
    K_new_old=K_new_old(length(predtimes)+1:2*length(predtimes),1:length(model.t));
  else
    K_new=K_new(length(predtimes{1})+1:length(predtimes{1})+length(predtimes{2}),length(predtimes{1})+1:length(predtimes{1})+length(predtimes{2}));
    K_new_old=K_new_old(length(predtimes{1})+1:length(predtimes{1})+length(predtimes{2}),1:length(model.t{1}));  
  end;

  if iscell(model.t)==0,
    K_old=model.K(1:length(model.t),1:length(model.t));
  else
    K_old=model.K(1:length(model.t{1}),1:length(model.t{1}));
  end;
  K_old_new=K_new_old';
end;


if iscell(model.t)==0,
  tempm=model.m(1:length(model.t));
else
  tempm=model.m(1:length(model.t{1}));
end;


priormeans=rnapriormeans;
size(tempm)
size(K_new_old)
size(K_old)
size(priormeans)
posteriormeans=priormeans+K_new_old*(K_old\tempm);
%covmatrix=K_new-K_new_old*inv(K_old)*K_old_new;

covmatrix=K_new-K_new_old*(K_old\K_old_new);

%posteriormeans=real(posteriormeans);
%covmatrix=real(covmatrix);
