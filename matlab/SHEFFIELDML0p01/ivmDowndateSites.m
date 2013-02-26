function model = ivmDowndateSites(model, index)

% IVMDOWNDATESITES Downdate site parameters.
%
%	Description:
%
%	MODEL = IVMDOWNDATESITES(MODEL, INDEX) remove a given data point
%	from the site parameters.
%	 Returns:
%	  MODEL - the model with the given point removed.
%	 Arguments:
%	  MODEL - the IVM structure from which the point is to be removed.
%	  INDEX - the index of the point to be removed.
%	
%
%	See also
%	IVMEPUPDATEM, IVMREMOVEPOINT


%	Copyright (c) 2005 Neil D. Lawrence


model.m(index, :) = 0;
model.beta(index, :) = 0;