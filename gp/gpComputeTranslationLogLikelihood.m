function Trans = gpComputeTranslationLogLikelihood(model, Xlabel, verbose);

% GPCOMPUTETRANSLATIONLOGLIKELIHOOD  
% FORMAT
% DESC  
% ARG model : 
% ARG Xlabel : 
% ARG verbose :
% RETURN Trans :
%
% COPYRIGHT : Carl Henrik Ek, 2010

% GP


if(iscell(Xlabel)&&size(Xlabel{1},1)==1)
  Xlabel = cell2mat(Xlabel);
end

if(verbose)
  handle_waitbar = waitbar(0,'Computing Translation Probabilities');
end
if(~iscell(Xlabel))
  Trans = zeros(size(Xlabel,1),size(Xlabel,1));
  for(i = 1:1:size(Xlabel,1))
    for(j = 1:1:size(Xlabel,1))
      Trans(i,j) = modelPointLogLikelihood(model,Xlabel(i,:),Xlabel(j,:));
    end
    if(verbose)
      waitbar(i/model.N);
    end
  end
else
  Trans = zeros(size(Xlabel{1},1),size(Xlabel{1},1),size(Xlabel,1)-1);
  for(i = 1:1:size(Xlabel{1},1))
    for(j = 1:1:size(Xlabel{1},1))
      for(k = 1:1:size(Xlabel,1)-1)
        Trans(i,j,k) = modelPointLogLikelihood(model,Xlabel{k}(i,:),Xlabel{k+1}(j,:));
      end
    end
    if(verbose)
      waitbar(i/size(Xlabel{1},1));
    end
  end
end
if(verbose)
  close(handle_waitbar);
end
  
return;
