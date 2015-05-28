function robOneDynamicsDisplay(model, spaceNum)

% ROBONEDYNAMICSDISPLAY Display the robot dynamics model. 

% FGPLVM

if nargin > 1
  spacing = repmat(32, 1, spaceNum);
else
  spaceNum = 0;
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('First Tailored dynamics model for robot:\n')
fprintf(spacing);
fprintf('  Length distribution value of a: %2.4f\n', model.a)
fprintf(spacing);
fprintf('  Length distribution value of b: %2.4f\n', model.b)
fprintf(spacing);
fprintf('  Length distribution mix coefficient: %2.4f\n', model.mixR);
fprintf(spacing);
fprintf('  Angle distribution sigma2: %2.4f\n', model.sigma2);
fprintf(spacing);
fprintf('  Angle distribution mix coefficient %2.4f\n', model.mixTheta);
