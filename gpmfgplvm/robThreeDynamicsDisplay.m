function robThreeDynamicsDisplay(model, spaceNum)

% ROBTHREEDYNAMICSDISPLAY Display the robot dynamics model. 

% FGPLVM

if nargin > 1
  spacing = repmat(32, 1, spaceNum);
else
  spaceNum = 0;
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Third Tailored dynamics model for robot:\n')
fprintf(spacing);
fprintf('  Lambda value: %2.4f\n', model.lambda)
fprintf(spacing);
fprintf('  Sigma value: %2.4f\n', model.sigma2)
