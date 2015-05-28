function Obs = gpComputeObservationLogLikelihood(model, X, Y, verbose)

% GPCOMPUTEOBSERVATIONLOGLIKELIHOOD  
% FORMAT
% DESC 
% ARG model : the model
% ARG X :
% ARG Y : 
% ARG verbose :
% RETURN Obs :
%
% COPYRIGHT : Carl Henrik Ek, 2010

% GP
  
  if(verbose)
    handle_waitbar = waitbar(0,'Computing Observation Loglikelihood');
  end
  if(~iscell(X))
    % stationary label set
    Obs = zeros(size(X,1),size(Y,1));
    for(label = 1:1:size(X,1))
      for(observation = 1:1:size(Y,1))
        Obs(label,observation) = gpPointLogLikelihood(model,X(label,:),Y(observation,:));
      end
      if(verbose)
        waitbar(label/size(X,1));
      end
    end
  else
    % non-stationary label set
    Obs = zeros(size(X{1},1),size(Y,1));
    for(label = 1:1:size(X{1},1))
      for(observation = 1:1:size(Y,1))
        Obs(label,observation) = gpPointLogLikelihood(model,X{observation}(label,:),Y(observation,:));
      end
      if(verbose)
        waitbar(label/size(X{1},1));
      end
    end
  end
  if(verbose)
    close(handle_waitbar);
  end

  return;
end
