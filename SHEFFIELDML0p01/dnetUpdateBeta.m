function model = dnetUpdateBeta(model)

% DNETUPDATEBETA Do an M-step (update parameters) on an Density Network model.
%
%	Description:
%
%	MODEL = DNETUPDATEBETA(MODEL) updates the parameters of a Density
%	Network model.
%	 Returns:
%	  MODEL - the model with updated parameters.
%	 Arguments:
%	  MODEL - the model which is to be updated.
%	
%
%	See also
%	DNETCREATE, DNETUPDATEOUTPUTWEIGHTS


%	Copyright (c) 2008 Neil D. Lawrence


  Ypred = dnetOut(model);
  tot = 0;
  if model.N > model.M
    for i = 1:model.M
      diffY = model.y - repmat(Ypred(i, :), model.N, 1);
      diffY =diffY.*diffY.*repmat(model.w(:,i), 1, model.d);
      tot = tot + sum(sum(diffY));
    end
  else
    for i = 1:model.N
      diffY = repmat(model.y(i, :), model.M, 1) - Ypred;
      diffY = diffY.*diffY.*repmat(model.w(i,:)', 1, model.d);
      tot = tot + sum(sum(diffY));
    end
  end
  model.beta = (model.N*model.d)/tot;
