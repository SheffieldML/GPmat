MOCAP software
Version 0.11		Monday 05 Sep 2005 at 15:19
Copyright (c) 2005 Neil D. Lawrence

This toolbox allows MATLAB to read in and write bvh files. There are also routines for visualising the files in MATLAB.


Version 0.11 Release Notes
--------------------------

There were some missing files in the previous version that have been added in this version.

Version 0.1 Release Notes
-------------------------

First release of the toolbox to coincide with release of C++ GPLVM code.


MATLAB Files
------------

Matlab files associated with the toolbox are:

bvh2xyz.m: Compute XYZ values given structure and channels.
bvhConnectionMatrix.m: Compute the connection matrix for the structure.
bvhModify.m: Helper code for visualisation of bvh data.
bvhPlayFile.m: Play motion capture data from a bvh file.
bvhReadFile.m: Reads a bvh file into a tree structure.
bvhVisualise.m: For updating a bvh representation of 3-D data.
bvhWriteFile.m: Write a bvh file from a given structure and channels.
mocapResultsCppBvh.m: Load results from cpp file and visualise as a bvh format.
rotationMatrix.m: Compute the rotation matrix for an angle in each direction.
stickModify.m: Helper code for visualisation of a stick man.
stickVisualise.m: For updateing a stick representation of 3-D data.
bvhPlayData.m: Play bvh motion capture data.
