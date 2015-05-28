function [priormeans,posteriormeans,covmatrix] = gpnddisimPredictRNAOnly(model,predtimes);

% GPASIMPREDICT Compute predictions (means and a covariance matrix)
% of RNA values for the GPASIM model, conditional on the existing
% RNA values in the model but not on any existing POL2 values.
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
if numGenes==0,
  error('gpnddisimPredictRNAOnly requires a model that has parameters for RNA modeling\n');
end;


% compute prior means 
%pol2priormeans=[1:size(predtimes,1)]'*0;
pol2priormeans=ones(size(predtimes,1),1)*model.simMean;

% Mean for the mRNA is nonconstant over time and depends on the
% B,D,S parameters and on the POL2 mean
Bj=model.B(1);
Dj=model.D(1);
Sj=model.S(1);
if model.use_disimstartmean==1,
  disimStartMean=model.disimStartMean(1);
end;  


% compute the RNA mean curve

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


%if model.use_disimstartmean==1,
%  rnapriormeans=(Bj+model.simMean*Sj)/Dj+(disimStartMean-(Bj+model.simMean*Sj)/Dj)*exp(Dj*(-predtimes));
%else    
%  rnapriormeans=(Bj+model.simMean*Sj)/Dj*ones(size(predtimes));
%end;


if 1,
K_new=kernCompute(model.kern, predtimes);
K_new=K_new(length(predtimes)+1:2*length(predtimes),length(predtimes)+1:2*length(predtimes));

K_new_old=kernCompute(model.kern, predtimes, model.t);
K_new_old=K_new_old(length(predtimes)+1:2*length(predtimes),length(model.t)+1:2*length(model.t));

K_old=model.K(length(model.t)+1:2*length(model.t),length(model.t)+1:2*length(model.t));
K_old_new=K_new_old';
end;


tempm=model.m(length(model.t)+1:2*length(model.t));


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
