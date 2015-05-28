function xyzankurAnimCompare(X, X2, fps)

% XYZANKURANIMCOMPARE Animate a prediction and ground truth for stick man from Agarwal & Triggs dataset.
%
%	Description:
%
%	XYZANKURANIMCOMPARE(Y, Y2, FPS) animates a matrix of x,y,z point
%	clound positions representing the motion of the figure used to
%	generate the silhouttes for Agarwal & Triggs silhouette data.
%	 Arguments:
%	  Y - the data to animate.
%	  Y2 - other data to animate.
%	  FPS - the number of frames per second to animate (defaults to 24).
%	
%
%	See also
%	XYZANKURVISUALISE, XYZANKURMODIFY


%	Copyright (c) 2008 Carl Henrik Ek and Neil D. Lawrence


  
  if(nargin<3)
    fps = 24;
    if(nargin<2)
      fid = 1;
      if(nargin<1)
        error('Too few arguments');
      end
    end
  end
  clf
  for(i = 1:1:size(X,1))
    if(i==1)
      subplot(1, 2, 1)
      handle = xyzankurVisualise(X(i,:), 1);
      subplot(1, 2, 2)
      handle2 = xyzankurVisualise(X2(i,:), 1);
    else
      xyzankurModify(handle, X(i,:));
      xyzankurModify(handle2, X2(i,:));
    end
    pause(1/fps);
  end
end