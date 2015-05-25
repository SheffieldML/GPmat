function model = fgplvmExpandParam(model, params)

% FGPLVMEXPANDPARAM Expand a parameter vector into a GP-LVM model.
% FORMAT
% DESC takes an FGPLVM structure and a vector of parameters, and
% fills the structure with the given parameters. Also performs any
% necessary precomputation for likelihood and gradient
% computations, so can be computationally intensive to call.
% ARG model : the FGPLVM structure to put the parameters in.
% ARG params : parameter vector containing the parameters to put in
% the FGPLVM structure.
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2009
%
% MODIFICATION : Carl Henrik Ek, 2009
% 
% SEEALSO : fgplvmCreate, fgplvmExtractParam, modelExpandParam

% FGPLVM


startVal = 1;
if isfield(model, 'back') & ~isempty(model.back)
  % update modelParameters
  endVal = model.back.numParams;
  model.back = modelExpandParam(model.back, params(startVal:endVal));

  % update latent locations
  tmp = modelOut(model.back,model.y);
  tmp_dim = 1;
  for(i = 1:1:model.q)
    if(length(find(model.back.indexOut==i))~=0)
      model.X(:,model.back.indexOut(tmp_dim)) = tmp(:,tmp_dim);
      tmp_dim = tmp_dim + 1;
    else
      startVal = endVal + 1;
      endVal = endVal + model.N;
      model.X(:,i) = reshape(params(startVal:endVal),model.N,1);
    end
  end
  clear tmp tmp_dim;
else
  endVal = model.N*model.q;
  model.X = reshape(params(startVal:endVal), model.N, model.q);
end
startVal = endVal+1;
endVal = endVal + model.kern.nParams;

switch model.approx
 case 'ftc'
  endVal = endVal;
 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
  if model.fixInducing
    % account for beta attached to the end.
    endVal = endVal + 1; 
    % X_u values are taken from X values.
    model.X_u = model.X(model.inducingIndices, :);
  else
    % Parameters include inducing variables and beta.
    endVal = endVal + model.q*model.k + 1;
  end

 otherwise
  error('Unknown approximation type.')
end
if model.learnScales
  endVal = endVal + model.d;
end
model = gpExpandParam(model, params(startVal:endVal));


% Give parameters to dynamics if they are there.
if isfield(model, 'dynamics') & ~isempty(model.dynamics)
  startVal = endVal + 1;
  endVal = length(params);

  % Fill the dynamics model with current latent values.
  model.dynamics = modelSetLatentValues(model.dynamics, model.X);

  % Update the dynamics model with parameters (thereby forcing recompute).
  model.dynamics = modelExpandParam(model.dynamics, params(startVal:endVal));
end

% Constraints
if(isfield(model,'constraints')&&~isempty(model.constraints))
  for(i = 1:1:model.constraints.numConstraints)
    model.constraints.comp{i} = constraintExpandParam(model.constraints.comp{i},model.X);
  end
end
