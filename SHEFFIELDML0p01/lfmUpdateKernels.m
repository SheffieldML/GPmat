function model = lfmUpdateKernels(model)

% LFMUPDATEKERNELS Updates the kernel representations in the LFM structure.
%
%	Description:
%
%	MODEL = LFMUPDATEKERNELS(MODEL) updates any representations of the
%	kernel in the model structure, such as invK, logDetK or K.
%	 Returns:
%	  MODEL - the model structure with the kernels updated.
%	 Arguments:
%	  MODEL - the model structure for which kernels are being updated.
%	
%
%	See also
%	% SEEALSO LFMEXPANDPARAM, LFMCREATE


%	Copyright (c) 2007 % COPYRIGHT Neil D. Lawrence



model.K = real(kernCompute(model.kern, model.t)); 
if isfield(model, 'yvar')
  model.K = model.K + diag(model.yvar);
end
[model.invK, U, jitter] = pdinv(model.K);
if jitter>1e-4
  fprintf('Warning: lfmUpdateKernels added jitter of %2.4f\n', jitter)
end
model.logDetK = logdet(model.K, U);