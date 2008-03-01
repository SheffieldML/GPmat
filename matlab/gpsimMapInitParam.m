function options = gpsimMapInitParam(y, yvar, options);
  
% GPSIMMAPINITPARAM Creates a set of options for a GPSIMMAP model as the initial parameters.
% FORMAT
% DESC returns a default initialisation for a GPSIMMAP model.
% ARG y :  the values of each gene at the different time points.
% ARG yvar :  the variances of each gene at the different time points.
% ARG options : the values of each gene at the different time points.
% RETURN options : the options structure.
% 
% SEEALSO : gpsimMapCreate
%
% COPYRIGHT : Pei Gao, Neil D. Lawrence and Magnus Rattary, 2008

% GPSIM  
  
 
numGenes = length(options.S);
options.S = ones(1, numGenes);
initOptions = options;     
mu = mean(y, 1);
optionDset = [0.01 0.1 0.4 1];
rand('seed',1);

for j = 1:numGenes
  initOptions.S = options.S(j);
  initOptions.gParam = options.gParam(j);     
  %    figure(j);
  %    plot(options.times, y(:,j), 'rx');
  %    hold on;
  for di = 1:length(optionDset)
    initOptions.D = optionDset(di);
    initOptions.B = initOptions.D*mu(j);
    
    model.gene{j}.init{di} = gpsimMapCreate(1, 1, options.times, y(:,j), ...
                                            yvar(:,j), initOptions);
    if strcmp(options.kern, 'mlp')
      model.gene{j}.init{di}.kern.weightVariance = 30;
      model.gene{j}.init{di}.kern.biasVariance = 20;
      % This forces kernel recompute.
      params = gpsimMapExtractParam(model.gene{j}.init{di});
      model.gene{j}.init{di} = ...
          gpsimMapExpandParam(model.gene{j}.init{di}, params);
    elseif strcmp(options.kern, 'rbf')
      model.gene{j}.init{di}.kern.inverseWidth = 0.001;
      params = gpsimMapExtractParam(model.gene{j}.init{di});
      model.gene{j}.init{di} = ...
          gpsimMapExpandParam(model.gene{j}.init{di}, params);
    end    
    
    updateFoptions = defaultOptions;
    model.gene{j}.init{di} = gpsimMapUpdateF(model.gene{j}.init{di}, ...
                                             updateFoptions);
    %       fmean =  mean(model.gene{j}.init{di}.f);
    %       fsd = sqrt(var(model.gene{j}.init{di}.f));
    %       model.gene{j}.init{di}.f = model.gene{j}.init{di}.f - fmean;
    %       model.gene{j}.init{di}.f = model.gene{j}.init{di}.f/fsd;
    %       model.gene{j}.init{di}.S = model.gene{j}.init{di}.S/exp(fmean);
    %       model.gene{j}.init{di}.gParam = model.gene{j}.init{di}.gParam/exp(fmean);
    ypredMean = mean(model.gene{j}.init{di}.ypred);
    ypredScale = sqrt(var(model.gene{j}.init{di}.ypred(model.gene{j}.init{di}.times_index)));
    model.gene{j}.init{di}.B = (model.gene{j}.init{di}.B - initOptions.D* ...
                                (ypredMean - ypredScale*mean(model.gene{j}.init{di}.y)))/ ...
        ypredScale;
    %       if isfield(model.gene{j}.init{di}, 'alpha')
    %         model.gene{j}.init{di}.alpha = model.gene{j}.init{di}.alpha/ypredScale;
    %       end
    if model.gene{j}.init{di}.B < 0 && isfield(model.gene{j}.init{di}, ...
                                               'alpha')
      model.gene{j}.init{di}.alpha = model.gene{j}.init{di}.alpha/ypredScale;
      model.gene{j}.init{di}.alpha = model.gene{j}.init{di}.alpha + ...
          model.gene{j}.init{di}.B;
      model.gene{j}.init{di}.B = 1e-6;
    end
    model.gene{j}.init{di}.S = model.gene{j}.init{di}.S/ypredScale; 
    params = gpsimMapExtractParam(model.gene{j}.init{di});  
    model.gene{j}.init{di} = gpsimMapExpandParam(model.gene{j}.init{di}, ...
                                                 params);
    f = gpsimMapFunctionalExtractParam(model.gene{j}.init{di});
    model.gene{j}.init{di} = gpsimMapFunctionalExpandParam(model.gene{j}.init{di}, f);
    model.gene{j}.ll(di) = ...
        gpsimMapLogLikelihood(model.gene{j}.init{di});
    %      plot(model.gene{j}.init{di}.mapt, model.gene{j}.init{di}.ypred);
  end
  [value, index] = max(model.gene{j}.ll);
  %    plot(model.gene{j}.init{index}.mapt, model.gene{j}.init{index}.ypred, 'm');    
  
  options.D(j) = optionDset(index);
  options.S(j) = model.gene{j}.init{index}.S;
  %     if model.gene{j}.init{index}.B < 0 
  %       if isfield(model.gene{j}.init{index},'alpha')
  %         options.alpha(j) = model.gene{j}.init{index}.alpha + ...
  %             model.gene{j}.init{di}.B;
  %       end
  %       options.B(j) = 1e-6;
  %     else
  if isfield(model.gene{j}.init{index},'alpha')
    options.alpha(j) = model.gene{j}.init{index}.alpha;  
  end      
  %     end
  options.B(j) = model.gene{j}.init{index}.B;     
  options.gParam(j) = model.gene{j}.init{index}.gParam;
end
options.S = rand(1,numGenes);
options.B = options.D.*mu;
  