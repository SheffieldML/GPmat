function model = spectralUpdateX(model)
 
% SPECTRALUPDATEX Update the latent representation for spectral model.
% FORMAT
% DESC updates the latent represenation given the stiffness matrix for a
% spectral model (Laplacian matrix is stored in model.L).
% ARG model : the model to be updated (with the stiffness matrix
% computed).
% RETURN model : the model with latent positions X updated.
%
% COPYRIGHT : Neil D. Lawrence, 2009
%
% SEEALSO : leOptimise, lleOptimise, mvuOptimise

% MLTOOLS
  
  options.disp = 0; 
  options.isreal = 1; 
  options.issym = 1; 
  noeigs = false;
  if isoctave
    noeigs = true;
    %verStr= ver('octave').Version
    %if str2num(verStr(1:3))<3.2
    %  noeigs = true;
    %end
  end
  if noeigs
    warning('No eigs function in Octave');
    % Nasty hack for eigenvalue problem in Octave.
    [m, v] = eig(model.L);
    [v, order] = sort(diag(v));
    v = diag(v(1:model.q+1));
    m = m(:, order);
    m = m(:, 1:model.q+1);
  else
    [m, v] = eigs_r11(model.L, model.q+1, 'sm', options);
  end
  v = diag(v);
  if ~isfield(model, 'discardLowest') || model.discardLowest
    [void, ind] = min(v);
  else
    [void, ind] = max(v);
  end

  % Multiplying by square root ensures latent covariance of identity.
  model.X = m(:, [1:(ind-1) (ind+1):end])*sqrt(model.N);
  
end
