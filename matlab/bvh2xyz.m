function xyz = bvh2xyz(bvhStruct, channels)

% BVHCOMPUTEXYZ Compute XYZ values given structure and channels.

% MOCAP

for i = 1:length(bvhStruct)  
  if ~isempty(bvhStruct(i).posInd)
    xpos = channels(bvhStruct(i).posInd(1));
    ypos = channels(bvhStruct(i).posInd(2));
    zpos = channels(bvhStruct(i).posInd(3));
  else
    xpos = 0;
    ypos = 0;
    zpos = 0;
  end
  xyzStruct(i) = struct('rotation', [], 'xyz', []); 
  if nargin < 2 | isempty(bvhStruct(i).rotInd)
    xangle = 0;
    yangle = 0;
    zangle = 0;
  else
    xangle = deg2rad(channels(bvhStruct(i).rotInd(1)));
    yangle = deg2rad(channels(bvhStruct(i).rotInd(2)));
    zangle = deg2rad(channels(bvhStruct(i).rotInd(3)));
  end
  thisRotation = rotationMatrix(xangle, yangle, zangle, bvhStruct(i).order);
  thisPosition = [xpos ypos zpos];
  if ~bvhStruct(i).parent
    xyzStruct(i).rotation = thisRotation;
    xyzStruct(i).xyz = bvhStruct(i).offset + thisPosition;
  else
    xyzStruct(i).xyz = ...
        bvhStruct(i).offset*xyzStruct(bvhStruct(i).parent).rotation ...
        + xyzStruct(bvhStruct(i).parent).xyz + thisPosition;
    xyzStruct(i).rotation = thisRotation*xyzStruct(bvhStruct(i).parent).rotation;
    
  end
end
xyz = reshape([xyzStruct(:).xyz], 3, length(bvhStruct))';



function theta = deg2rad(omega)

theta = omega/180*pi;