function rotMat = rotationMatrixGradient(xangle, yangle, zangle, order, gradientRespect)

% ROTATIONMATRIXGRADIENT Compute the gradient of rotation with respect to one angle.
% FORMAT
% DESC computes the gradient of the rotation matrix with respect to one
% of the rotation angles.  
% RETURN g : the gradient of the rotation matrix with respect to the
% given angle.
% ARG xangle : angle of rotation around x-axis.
% ARG yangle : angle of rotation around y-axis.
% ARG zangle : angle of rotation around z-axis.
% ARG order : order in which rotations is applied.
% ARG gradientRespect : angle which gradient is with respect to.
%
% COPYRIGHT : Romain Floyrac, 2009
%
% MODIFICATIONS : Neil D. Lawrence, 2009

% MOCAP 
  
  if nargin < 4
    order = 'zxy';
  end
  if isempty(order)
    order = 'zxy';
  end
  % Here we assume we rotate z, then x then y.
  c1 = cos(xangle); % The x angle
  c2 = cos(yangle); % The y angle
  c3 = cos(zangle); % the z angle
  s1 = sin(xangle);
  s2 = sin(yangle);
  s3 = sin(zangle);
  
  % see http://en.wikipedia.org/wiki/Rotation_matrix for
  % additional info.
  
  rotMat = eye(3);
  for i = 1:length(order)
    switch order(i)
     case 'x'
      if (gradientRespect == 'x')
        rotMat = [0 0 0
                  0  -s1 c1
                  0 -c1 -s1]*rotMat;
        
      else
        rotMat = [1 0 0
                  0  c1 s1
                  0 -s1 c1]*rotMat;
      end
     case 'y'
      if (gradientRespect == 'y')
        rotMat = [-s2 0 -c2
                  0 0 0
                  c2 0 -s2]*rotMat;
        
      else
        rotMat = [c2 0 -s2
                  0 1 0
                  s2 0 c2]*rotMat;
      end
     case 'z' 
      if (gradientRespect == 'z')
        rotMat = [-s3 c3 0
                  -c3 -s3 0
                  0 0 0]*rotMat;       
        
      else
        rotMat = [c3 s3 0
                  -s3 c3 0
                  0 0 1]*rotMat;
      end
      
    end
    
  end
end