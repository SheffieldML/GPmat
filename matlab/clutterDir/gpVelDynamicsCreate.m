function model = gpVelDynamicsCreate(q, d, latentVals, kern, approx, numActive, ...
                                  diff, learn, learnScales)

% GPDYNAMICSCREATE Create the dynamics model. 


if nargin < 9
  learnScales = 0;
  if nargin < 8
    learn = 0;
    if nargin < 7
      diff = 1;
      if nargin < 6
        numActive = 100;
        if nargin < 5
          approx = 'ftc';
          if nargin < 4
            kern = {'rbf', 'white'};
          end
        end
      end
    end
  end
end

X = latentVals(1:end-1, :);
y = latentVals(2:end, :);
if diff
  y = y - X;
end
model = gpCreate(q, d, X, y, kern, approx, numActive);

model.diff = diff;
model.learn = learn;
model.learnScales = learnScales;
model.type = 'gpDynamics';

