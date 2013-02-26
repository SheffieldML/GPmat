function xyzankurAnimCompareMultiple(X, X2, fps, names)

% XYZANKURANIMCOMPAREMULTIPLE Animate many predictions and ground truth for
%
%	Description:
%	stick man from Agarwal & Triggs dataset.
%
%	XYZANKURANIMCOMPAREMULTIPLE(X, X2, FPS) animates a matrix of x,y,z
%	point cloud positions representing the motion of the figure used to
%	generate the silhouttes for Agarwal & Triggs silhouette data.
%	 Arguments:
%	  X - the ground truth to animate.
%	  X2 - Cell of other predictions to animate.
%	  FPS - the number of frames per second to animate (defaults to 24).
%	
%	
%
%	See also
%	XYZANKURVISUALISE, XYZANKURMODIFY


%	Copyright (c) 2008 Carl Henrik Ek and Neil D. Lawrence


%	With modifications by Alfredo A. Kalaitzis 2011


nDatasets = length(X2);
more_handles = zeros(nDatasets,19);
if(nargin < 3)
    fps = 24;
    if(nargin < 1)
        error('Too few arguments');
    end
end

clf
for i = 1:1:size(X,1)
    if (i == 1)
        subplot(1, nDatasets+1, 1), title(names{1})
        handle = xyzankurVisualise(X(i,:), 1);
        for j = 1:nDatasets
            subplot(1, nDatasets+1, j+1), title(names{j+1})
            more_handles(j,:) = xyzankurVisualise(X2{j}(i,:), 1);
        end
    else
        xyzankurModify(handle, X(i,:));
        for j = 1:nDatasets
            xyzankurModify(more_handles(j,:), X2{j}(i,:));
        end
    end
    pause(1/fps);
end
