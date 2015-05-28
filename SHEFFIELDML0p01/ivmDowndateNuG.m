function model = ivmDowndateNuG(model, index)

% IVMDOWNDATENUG Downdate nu and g parameters associated with noise model.
%
%	Description:
%
%	MODEL = IVMDOWNDATENUG(MODEL, INDEX) removes a given data point from
%	the nu and g representations.
%	 Returns:
%	  MODEL - the model with the point removed.
%	 Arguments:
%	  MODEL - the model from which the point is to be removed.
%	  INDEX - the index of the point to be removed.
%	
%
%	See also
%	IVMREMOVEPOINT, IVMEPUPDATEPOINT, IVMDOWNDATEM


%	Copyright (c) 2005 Neil D. Lawrence


if isfield(model.kern, 'whiteVariance')
  subTerm = model.varSigma(index, :)-model.kern.whiteVariance;
else
  subTerm = model.varSigma(index, :);
end
subTerm = model.varSigma(index, :);
model.nu(index, :) = 1./(1./model.beta(index, :) - subTerm);
model.g(index, :) = model.nu(index, :).*(model.mu(index, :) - model.m(index, :));

