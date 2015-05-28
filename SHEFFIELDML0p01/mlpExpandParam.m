function model = mlpExpandParam(model, params)

% MLPEXPANDPARAM Update mlp model with new vector of parameters.
%
%	Description:
%
%	MODEL = MLPEXPANDPARAM(MODEL, PARAMS) takes a vector of MLP weights
%	and places them in their respective positions in the MLP model. For
%	single hidden layer neural networks the function is a wrapper for
%	the mlpunpak command.
%	 Returns:
%	  MODEL - the model with the weights distributed in the correct
%	   places.
%	 Arguments:
%	  MODEL - the model in which the weights are to be placed.
%	  PARAMS - a vector of the weights to be placed in the model.
%	
%
%	See also
%	MLPUNPAK, MLPCREATE, MLPEXTRACTPARAM


%	Copyright (c) 2006, 2007 Neil D. Lawrence


if length(model.hiddenDim) == 1
  model = mlpunpak(model, params);
else
  startVal = 1;
  endVal = model.inputDim*model.hiddenDim(1);
  model.w{1} = reshape(params(startVal:endVal, model.inputDim, ...
                            model.hiddenDim(1)));
  startVal = endVal + 1;
  endVal = endVal + model.hiddenDim(1);
  model.b{1} = params(startVal:endVal);
  for i = 2:length(model.hiddenDim)
    startVal = endVal + 1;
    endVal = endVal + model.hiddenDim(i-1)*model.hiddenDim(i);
    model.w{i} = reshape(params(startVal:endVal), model.hiddenDim(i-1), ...
                                model.hiddenDim(i));
    startVal = endVal + 1;
    endVal = endVal + model.hiddenDim(i);
    model.b{i} = params(startVal:endVal);
  end
  i = length(model.hiddenDim);
  startVal = endVal + 1;
  endVal = endVal + model.hiddenDim(i)*model.outputDim;
  model.w{i+1} = resphape(params(startVal:endVal), model.hiddenDim(i), ...
                          model.outputDim);
  startVal = endVal + 1;
  endVal = endVal + model.outputDim;
  model.b{i+1} = params(startVal:endVal);
end
