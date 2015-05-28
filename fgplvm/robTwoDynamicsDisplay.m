function robTwoDynamicsDisplay(model, spaceNum)

% ROBTWODYNAMICSDISPLAY Display the robot dynamics model. 

% FGPLVM

if nargin > 1
  spacing = repmat(32, 1, spaceNum);
else
  spaceNum = 0;
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Second Tailored dynamics model for robot:\n')
fprintf(spacing);
fprintf('  Length distribution value of a: %2.4f\n', model.a)
fprintf(spacing);
fprintf('  Length distribution value of b: %2.4f\n', model.b)
fprintf(spacing);
fprintf('  Length distribution mix coefficient: %2.4f\n', model.mixR);
fprintf(spacing);
fprintf('  Angle distribution component 1 sigma2: %2.4f\n', model.sigma2One);
fprintf(spacing);
fprintf('  Angle distribution component 1 mix coefficient: %2.4f\n', model.mixThetaOne);
fprintf(spacing);
fprintf('  Angle distribution component 2 sigma2: %2.4f\n', model.sigma2Two);
fprintf(spacing);
fprintf('  Angle distribution component 2 mix coefficient: %2.4f\n', model.mixThetaTwo);
