function model = fgplvmAddConstraint(model,options)

% FGPLVMADDCONSTRAINT Add latent constraints to FGPLVM model
% FORMAT
% DESC Adds constraint structure to FGPLVM model
% ARG model : fgplvm model
% ARG options : options strucure as returnded from
% constraintOptions
% RETURN model : the GP-LVM model.
%
% COPYRIGHT : Carl Henrik Ek, 2009
%
% SEEALSO : constraintOptions

% FGPLVM


%
if(isfield(model,'constraints')&&~isempty(model.constraints))
  % constraints init
  model = constraintCreate(model,options);
else
  % no constraints
  model.constraints = [];
  model.constraints.q = model.q;
  model.constraints.N = model.N;
  model.constraints.id = [];
  model.constraints.numConstraints = 0;
  model.constraints.comp = {};
  model = constraintCreate(model,options);
end

model.constraints.numConstraints = model.constraints.numConstraints + 1;
model.constraints.id = [model.constraints.id; false*ones(1,model.q)];
model.constraints.id(end,options.dim) = true;