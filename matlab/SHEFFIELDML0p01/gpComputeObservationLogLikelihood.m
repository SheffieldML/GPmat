function Obs = gpComputeObservationLogLikelihood(model, X, Y, verbose)

% GPCOMPUTEOBSERVATIONLOGLIKELIHOOD
%
%	Description:
%
%	OBS = GPCOMPUTEOBSERVATIONLOGLIKELIHOOD(MODEL, X, Y, VERBOSE) 
%	 Returns:
%	  OBS - 
%	 Arguments:
%	  MODEL - the model
%	  X - 
%	  Y - 
%	  VERBOSE - 


%	Copyright (c) 2010 Carl Henrik Ek

  
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