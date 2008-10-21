function xyzankurAnim(X,fid,fps)

% XYZANKURANIM
%
% COPYRIGHT : Carl Henrik Ek and Neil Lawrence, 2008

% MOCAP


if(nargin<3)
  fps = 24;
  if(nargin<2)
    fid = 1;
    if(nargin<1)
      error('To Few Arguments');
    end
  end
end

for(i = 1:1:size(X,1))
  if(i==1)
    handle = xyzankurVisualise(X(i,:),1);
  else
    xyzankurModify(handle,X(i,:));
  end
  pause(1/fps);
end