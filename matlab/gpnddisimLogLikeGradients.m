function g = gpnddisimLogLikeGradients(model,varargin)

% GPNDDISIMLOGLIKEGRADIENTS Compute the gradients of the log likelihood of a GPNDDISIM model.
% FORMAT
% DESC computes the gradients of the log likelihood of the given
% Gaussian process for use in a single input motif protein network.
% ARG model : the model for which the log likelihood is computed.
% RETURN g : the gradients of the parameters of the model.
% 
% SEEALSO : gpsimCreate, gpsimLogLikelihood, gpsimGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007
%
% COPYRIGHT : Jaakko Peltonen, 2011
  
% GPSIM

if length(varargin)==2,
  update_kernel=varargin{1};
  update_mean=varargin{2};
else
  update_kernel=1;
  update_mean=1;
end;

covGrad = -model.invK + model.invK*model.m*model.m'*model.invK;
covGrad = 0.5*covGrad;


if (update_kernel==1),
  if isfield(model, 'proteinPrior') && ~isempty(model.proteinPrior)
    g = kernGradient(model.kern, model.timesCell, covGrad);
  else
    g = kernGradient(model.kern, model.t, covGrad);
  end
  
%  fprintf('gpdisimLogLikeGradients: g before any priors\n');
%   g'
%  pause
  
  % In case we need priors in.
  % Add contribution of any priors 
  if isfield(model, 'bprior'),
    g = g + kernPriorGradient(model.kern);
  end
else
  g = zeros(1,model.kern.nParams);
end;

  




if (update_mean==1),

  gmuFull = model.m'*model.invK;

  % if isfield(model, 'proteinPrior') && ~isempty(model.proteinPrior)
  %   if model.includeNoise
  %     ind = model.kern.comp{1}.diagBlockDim{1} + (1:model.kern.comp{1}.diagBlockDim{2});
  %     gmu = zeros(size(1, model.numGenes));
  
  %     for i = 1:model.numGenes
  %       gmu(i) = sum(gmuFull(ind));
  %       ind = ind + model.kern.comp{1}.diagBlockDim{i+1};
  %     end
  %   else
  %     ind = model.kern.diagBlockDim{1} + (1:model.kern.diagBlockDim{2});
  %     gmu = zeros(size(1, model.numGenes));
  
  %     for i = 1:model.numGenes
  %       gmu(i) = sum(gmuFull(ind));
  %       ind = ind + model.kern.diagBlockDim{i+1};
  %     end
  %   end
  
  % else
  
  if (model.numGenes>0) && (model.use_disimstartmean==1),
    %--------------------------------
    % Version with NDDISIM start mean
    %--------------------------------
    
    % PART1 Compute gradient for basal transcription rates
    gb=zeros(1, model.numGenes);    
    if iscell(model.t)==0, tempt=model.t; else tempt=model.t{1}; end;
    indStart=length(tempt)+1;
    for k=1:model.numGenes,
      % note that delay in the DISIM model does not affect the
      % contribution of the basal transcription, thus delays are
      % not applied to model.t here.
      if iscell(model.t)==0, delayedt=model.t; else delayedt=model.t{k+1}; end;
      % indStart=length(model.t)*k + 1;
      % indEnd=indStart+length(model.t)-1;
      indEnd=indStart+length(delayedt)-1;
      %size(delayedt)
      %indStart
      %indEnd
      if size(delayedt,1)>0,
        gb(k)=gmuFull(indStart:indEnd)*((1-exp(-model.D(k)*delayedt))/model.D(k));
      else
        gb(k)=0;
      end;
      indStart=indStart+size(delayedt,1);
    end;

    % Add contribution of prior on B if it exists.
    if isfield(model, 'bprior');
      if model.numGenes>0,
	gb = gb + priorGradient(model.bprior, model.B);
      end;
    end

    % Multiply by factors from parameter transformations
    if model.numGenes>0,
      for k=1:length(gb),
	gb(k) = gb(k)*doTransform(model.B(k), 'gradfact', model.bTransform(k));
      end;
    end;
    
    % PART2 Compute gradient for DISIM start mean
    gdisimstartmean=zeros(1, model.numGenes);
    if iscell(model.t)==0, tempt=model.t; else tempt=model.t{1}; end;
    timeshift = min(tempt);
    indStart=length(tempt)+1;
    for k=1:model.numGenes,
      % note that delay in the DISIM model does not affect the
      % decay of the starting RNA concentration, thus delays are
      % not applied to model.t here.
      if iscell(model.t)==0, delayedt=model.t; else delayedt=model.t{k+1}; end;
      %indStart=length(model.t)*k + 1;
      %indEnd=indStart+length(model.t)-1;
      indEnd=indStart+size(delayedt,1)-1;
      if size(delayedt,1)>0,
        gdisimstartmean(k)=gmuFull(indStart:indEnd)* ...
            exp(-model.D(k)*(delayedt - timeshift));
      else
        gdisimstartmean(k)=0;
      end;
      indStart=indStart+size(delayedt,1);
    end;

    % Add contribution of prior on disimStartMean if it exists.
    if isfield(model, 'disimStartMeanPrior');
      if model.numGenes>0,
	gdisimstartmean = gdisimstartmean + ...
	    priorGradient(model.disimStartMeanPrior, model.disimStartMean);
      end;
    end

    % Multiply by factors from parameter transformations
    if model.numGenes>0,
      for k=1:length(gdisimstartmean),
	gdisimstartmean(k) = gdisimstartmean(k)*doTransform(model.disimStartMean(k), 'gradfact', model.disimStartMeanTransform(k));
      end;
    end;
    
    % PART3 Compute gradient for DISIM-level decays
    gd=zeros(1, model.numGenes);
    if iscell(model.t)==0, tempt=model.t; else tempt=model.t{1}; end;
    timeshift = min(tempt);
    indStart=length(tempt)+1;
    for k=1:model.numGenes,
      % note that delay in the DISIM model does not affect the
      % decay of the starting RNA concentration, or the
      % contribution of the basal transcription rate, thus delays are
      % not applied to model.t for those parts. However, delays do
      % affect the contribution of the POL2 mean ("simMean"), so
      % delays must be applied when computing the gradient with
      % respect to decay for that part of the mean function.
      if iscell(model.t)==0, tempt=model.t; else tempt=model.t{k+1}; end;
      delayedt=tempt-model.delay(k);      
      I=find(delayedt<0);
      delayedt(I)=0;
      
      indEnd=indStart+size(delayedt,1)-1;
      if size(delayedt,1)>0,
        gd(k)=gmuFull(indStart:indEnd)*...
  	    ( -model.B(k)/(model.D(k)*model.D(k))*(1-exp(-model.D(k)*tempt)) ...
              +model.B(k)/model.D(k)*exp(-model.D(k)*tempt).*tempt ...
	      +model.disimStartMean(k)* ...
              exp(-model.D(k)*(tempt - timeshift)).*(-tempt + timeshift) ...
              -model.simMean*model.S(k)/(model.D(k)*model.D(k))*(1-exp(-model.D(k)*delayedt)) ...
              +model.simMean*model.S(k)/model.D(k)*exp(-model.D(k)*delayedt).*delayedt ...
              );
      else
        gd(k)=0;
      end;
      indStart=indStart+size(delayedt,1);
    end;
    % Apply factors from transformations, and add gradient of the
    % decays to the main decay-gradient from the kernel,
    if model.numGenes>0,
      decayIndices = model.disimdecayindices;
      for k=1:length(decayIndices),
	g(decayIndices(k)) = g(decayIndices(k)) ...
	    + gd(k)*doTransform(model.D(k), 'gradfact',model.disimdecaytransformation(k));
      end;
    end;

    % PART4 Compute gradient for SIM-level mean
    gsimmean=zeros(1, 1);
    indStart=1;
    if iscell(model.t)==0, tempt=model.t; else tempt=model.t{1}; end;
    indEnd=indStart+length(tempt)-1;
    if size(tempt,1)>0,
      gsimmean=gmuFull(indStart:indEnd)*(ones(indEnd-indStart+1,1));
    else
      gsimmean=0;
    end;
    indStart=indStart+length(tempt);
    for k=1:model.numGenes,
      % note that delay in the DISIM model affects the contribution 
      % of the POL2 mean ("simMean"), so delays must be applied to time 
      % points when computing the gradient with respect to simMean.
      if iscell(model.t)==0, tempt=model.t; else tempt=model.t{k+1}; end;
      delayedt=tempt-model.delay(k);
      I=find(delayedt<0);
      delayedt(I)=0;
      
      %indStart=length(model.t)*k + 1;
      %indEnd=indStart+length(model.t)-1;
      indEnd=indStart+size(delayedt,1)-1;
      if size(delayedt,1)>0,
        gsimmean=gsimmean+gmuFull(indStart:indEnd)*((1-exp(-model.D(k)*delayedt))*model.S(k)/model.D(k));
      else
        gsimmean=gsimmean+0;
      end;
      indStart=indStart+size(delayedt,1);
    end;

    % Add contribution of prior on simMean if it exists.
    if isfield(model, 'simMeanPrior');
      if model.numGenes>0,
	gsimmean = gsimmean + ...
	    priorGradient(model.simMeanPrior, model.simMean);
      end;
    end

    % Multiply by factors from parameter transformations
    gsimmean = gsimmean*doTransform(model.simMean, 'gradfact', model.simMeanTransform);

    % PART5 Compute gradient for DISIM-level variance
    gdisimvar=zeros(1, model.numGenes);
    if iscell(model.t)==0, tempt=model.t; else tempt=model.t{1}; end;
    indStart=length(tempt)+1;
    for k=1:model.numGenes,
      % note that delay in the DISIM model affects the contribution 
      % of the POL2 mean ("simMean"), so delays must be applied to time 
      % points when computing the gradient with respect to DISIM variance,
      % for the part of the mean function related to simMean.
      if iscell(model.t)==0, tempt=model.t; else tempt=model.t{k+1}; end;
      delayedt=tempt-model.delay(k);
      I=find(delayedt<0);
      delayedt(I)=0;

      % indStart=length(model.t)*k + 1;
      % indEnd=indStart+length(model.t)-1;
      indEnd=indStart+size(delayedt,1)-1;
      %size(gmuFull)
      %indStart
      %indEnd
      %size(delayedt)
      if size(delayedt,1)>0,
        gdisimvar(k)=gmuFull(indStart:indEnd)*...
              ((1-exp(-model.D(k)*delayedt))*model.simMean/model.D(k)*0.5/model.S(k));
      else
        gdisimvar(k)=0;
      end;
      indStart=indStart+size(delayedt,1);
    end;    
    % Multiply by factors from parameter transformations, and add gradient of the
    % DISIM-variances to the main DISIM-variance gradient from the kernel
    if model.numGenes>0,
      disimvarIndices = model.disimvarianceindices;
      for k=1:length(disimvarIndices),
        g(disimvarIndices(k)) = g(disimvarIndices(k)) ...
            + gdisimvar(k)*doTransform(model.S(k)*model.S(k), 'gradfact',model.disimvariancetransformation(k));
      end;
    end;    

    % PART6 Compute gradient for DISIM-level delay
    gdisimdelay=zeros(1, model.numGenes);
    if iscell(model.t)==0, tempt=model.t; else tempt=model.t{1}; end;
    indStart=length(tempt)+1;
    for k=1:model.numGenes,
      % note that delay in the DISIM model affects the contribution 
      % of the POL2 mean ("simMean"), so delays must be applied to time 
      % points when computing the gradient with respect to DISIM variance,
      % for the part of the mean function related to simMean.
      if iscell(model.t)==0, tempt=model.t; else tempt=model.t{k+1}; end;
      delayedt=tempt-model.delay(k);      
      I=find(delayedt<0);
      delayedt(I)=0;

      % indStart=length(model.t)*k + 1;
      % indEnd=indStart+length(model.t)-1;
      indEnd=indStart+size(delayedt,1)-1;
      if size(delayedt,1)>0,
        gdisimdelay(k)=gmuFull(indStart:indEnd)*...
              (-exp(-model.D(k)*delayedt).*(delayedt>0)*model.simMean*model.S(k));
      else
        gdisimdelay(k)=0;
      end;
      indStart=indStart+size(delayedt,1);
    
%          +(model.simMean*model.S(k)/model.D(k))*(1-exp(-model.D(k)*tempt));
    
    end;    
    % gdisimdelay

    % Multiply by factors from parameter transformations, and add gradient of the
    % DISIM-variances to the main DISIM-variance gradient from the kernel
    if model.numGenes>0,
      disimdelayIndices = model.disimdelayindices;
      for k=1:length(disimdelayIndices),
        g(disimdelayIndices(k)) = g(disimdelayIndices(k)) ...
            + gdisimdelay(k)*doTransform(model.delay(k), 'gradfact',model.disimdelaytransformation(k));
      end;
    end;    
    
  else
    %--------------------------------
    % Version without DISIM start mean
    % TODO: this version of the code does not yet take into account delays in the model!
    % TODO: it also does not allow cell-format time indices, assumes RNA and POL2 have the 
    % observation times
    %--------------------------------
    gdisimstartmean=[];
    
    numData = size(model.t, 1);
    ind = 1:numData;
    ind = ind + numData;
    %  gmu = zeros(size(1, model.numGenes));
    gmu = zeros(1, model.numGenes);
    for i = 1:model.numGenes
      gmu(i) = sum(gmuFull(ind));
      ind = ind + numData;
    end
    %end
    
    if model.numGenes>0,
      gb = gmu./model.D;
    end;
    
    % In case we need priors in.
    % Add prior on B if it exists.
    if isfield(model, 'bprior');
      if model.numGenes>0,
	gb = gb + priorGradient(model.bprior, model.B);
      end;
    end
  
    if model.numGenes>0,
      for k=1:length(gb),
	gb(k) = gb(k)*doTransform(model.B(k), 'gradfact', model.bTransform(k));
      end;
    end;
  
    % Account for decay in mean.
    % This is a nasty hack to add the influence of the D in the mean to
    % the gradient already computed for the kernel. This is all very
    % clunky and sensitive to changes that take place elsewhere in the
    % code ...
    if model.numGenes>0,
      gd = -gmu.*(model.B+model.simMean)./(model.D.*model.D);
    end;

  
    % Apply transformations for decay-gradient in mean, and add to
    % main decay-gradient. Warning: only tested for 1 decay
    % parameter, the indices here might be slightly
    % incorrect in the general case!
    if model.numGenes>0,
      %decayIndices = [5];
      %for i = 3:model.kern.numBlocks
      %  decayIndices(end+1) = decayIndices(end) + 2;
      %end 
      decayIndices = model.disimdecayindices;
      for k=1:length(decayIndices),
	g(decayIndices(k)) = g(decayIndices(k)) ...
	    + gd(k)*doTransform(model.D(k), 'gradfact',model.disimdecaytransformation(k));  
      end;
    end;
  
    % Compute gradient for SIM-level mean
    gsimmean=zeros(1, 1);
    indStart=1;
    indEnd=indStart+length(model.t)-1;
    gsimmean=gmuFull(indStart:indEnd)*(ones(indEnd-indStart+1,1));
    for k=1:model.numGenes,
      indStart=length(model.t)*k + 1;
      indEnd=indStart+length(model.t)-1;
      gsimmean=gsimmean+gmuFull(indStart:indEnd)*(ones(indEnd-indStart+1,1)*model.S(k)/model.D(k));
    end;

    % Add contribution of prior on simMean if it exists.
    if isfield(model, 'simMeanPrior');
      if model.numGenes>0,
	gsimmean = gsimmean + ...
	    priorGradient(model.simMeanPrior, model.simMean);
      end;
    end

    % Multiply by factors from parameter transformations
    gsimmean = gsimmean*doTransform(model.simMean, 'gradfact', model.simMeanTransform);

    % Compute gradient for DISIM-level variance
    gdisimvar=zeros(1, model.numGenes);
    for k=1:model.numGenes,
      indStart=length(model.t)*k + 1;
      indEnd=indStart+length(model.t)-1;
      gdisimvar(k)=gmuFull(indStart:indEnd)*...
            (ones(indEnd-indStart+1,1)*model.simMean/model.D(k)*0.5/model.S(k));
    end;    
    % Multiply by factors from parameter transformations, and add gradient of the
    % DISIM-variances to the main DISIM-variance gradient from the kernel
    if model.numGenes>0,
      disimvarIndices = model.disimvarianceindices;
      for k=1:length(disimvarIndices),
        g(disimvarIndices(k)) = g(disimvarIndices(k)) ...
            + gdisimvar(k)*doTransform(model.S(k)*model.S(k), 'gradfact',model.disimvariancetransformation(k));
      end;
    end;    
    
  end;  
  
else
  gb = zeros(1,model.numGenes);
  if (model.use_disimstartmean==1),
    gdisimstartmean = zeros(1,model.numGenes);
  else
    gdisimstartmean = [];
  end;  
  gsimmean = 0;
end;

%fprintf(1,'Gradient after modifications from mean terms:\n');
%g
%pause

if model.numGenes>0,
  g = [g gb gdisimstartmean gsimmean];
else
  g = [g gsimmean];  
end;

if isfield(model, 'fix')
  for i = 1:length(model.fix)
    g(model.fix(i).index) = 0;
  end
end
