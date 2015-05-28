function [passParams, passNames] = multimodelExtractParam(model)

% MULTIMODELEXTRACTPARAM Extract parameters from the MULTIMODEL model structure.
%
%	Description:
%
%	PARAM = MULTIMODELEXTRACTPARAM(MODEL) extracts parameters from the
%	multi-task learning wrapper model structure into a vector of
%	parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the model.
%	 Arguments:
%	  MODEL - the model structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = MULTIMODELEXTRACTPARAM(MODEL) extracts parameters
%	and parameter names from the multi-task learning wrapper model
%	structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the model.
%	  NAMES - cell array of strings containing names for each parameter.
%	 Arguments:
%	  MODEL - the model structure containing the parameters to be
%	   extracted.
%	
%	
%	
%
%	See also
%	% SEEALSO MULTIMODELCREATE, MULTIMODELEXPANDPARAM, MODELEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2007, 2008 Neil D. Lawrence


%	With modifications by Mauricio Alvarez 2008, 2009

  passParams = zeros(1, model.numParams);
  passNames = cell(1, model.numParams);
  if nargout > 1
    [receiveParams, receiveNames] = modelExtractParam(model.comp{1});
  else
    receiveParams = modelExtractParam(model.comp{1});
  end
  
  if numel(model.outputDim) == 1      
      %endVal = model.numParams - model.numSep;
      endVal = model.numParams - model.numSep*model.numModels;
      passParams(1:endVal) = receiveParams(model.sharedIndices);
      if nargout > 1
          % MAURICIO : I think this is wrong. But I didn't change because I'm not
          % sure
          % passNames{1:endVal} = receiveNames{model.sharedIndices};
          passNames(1:endVal) = receiveNames(model.sharedIndices);
      end
      if ~isempty(model.separateIndices)
          startVal = endVal + 1;
          endVal = endVal + model.numSep;
          passParams(startVal:endVal) = receiveParams(model.separateIndices);
          if nargout > 1
              passNames(startVal:endVal) = receiveNames(model.separateIndices);
              for j=startVal:endVal
                  passNames{j} = [model.type ' 1 '  passNames{j}];
              end
          end
          for i = 2:length(model.comp)
              startVal = endVal+1;
              endVal = endVal + model.numSep;
              if nargout > 1
                  [receiveParams, receiveNames] = modelExtractParam(model.comp{i});
              else
                  receiveParams = modelExtractParam(model.comp{i});
              end
              passParams(startVal:endVal) = receiveParams(model.separateIndices);
              if nargout > 1
                  passNames(startVal:endVal) = receiveNames(model.separateIndices);
                  for j=startVal:endVal
                      passNames{j} = [model.type ' ' num2str(i) ' '  passNames{j}];
                  end
              end
          end
      end
  else
      endVal = model.comp{1}.nParams;
      passParams(1:endVal) = receiveParams;
      if nargout > 1
          passNames(1:endVal) = receiveNames;
          for j=1:endVal
              passNames{j} = [model.type ' 1 '  passNames{j}];
          end
      end
      for i = 2:length(model.comp)
          startVal = endVal+1;
          endVal = endVal + model.comp{i}.nParams;
          if nargout > 1
              [receiveParams, receiveNames] = modelExtractParam(model.comp{i});
          else
              receiveParams = modelExtractParam(model.comp{i});
          end
          passParams(startVal:endVal) = receiveParams;
          if nargout > 1
              passNames(startVal:endVal) = receiveNames;
              for j=startVal:endVal
                  passNames{j} = [model.type ' ' num2str(i) ' '  passNames{j}];
              end
          end
      end
  end
