function [kern, noise, ivmInfo] = ivmDeconstruct(model)

% IVMDECONSTRUCT break IVM in pieces for saving.
%
%	Description:
%
%	[KERN, NOISE, IVMINFO] = IVMDECONSTRUCT(MODEL) takes an IVM model
%	structure and breaks it into component parts for saving.
%	 Returns:
%	  KERN - the kernel component of the IVM model.
%	  NOISE - the noise component of the IVM model.
%	  IVMINFO - a structure containing the other information from the
%	   IVM: what the active set is, what the inactive set is and what the
%	   site parameters are.
%	 Arguments:
%	  MODEL - the model that needs to be saved.
%	
%
%	See also
%	IVMRECONSTRUCT


%	Copyright (c) 2005 Neil D. Lawrence


kern = rmfield(model.kern, 'Kstore');
kern = rmfield(kern, 'diagK');
noise = model.noise;
ivmInfo.I = model.I;
ivmInfo.J = model.J;
ivmInfo.m = model.m;
ivmInfo.beta = model.beta;
ivmInfo.type = 'ivm';
