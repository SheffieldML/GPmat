function xyz = skel2xyz(skel, channels)

% SKEL2XYZ Compute XYZ values given skeleton structure and channels.

% MOCAP

fname = str2func([skel.type '2xyz']);
xyz = fname(skel, channels);