ReadMe file for the MOCAP toolbox version 0.1 Friday, Jul 29, 2005 at 13:29:08
Written by Neil D. Lawrence.

License Info
------------

This software is free for academic use. Please contact Neil Lawrence if you are interested in using the software for commercial purposes.

This software must not be distributed or modified without prior permission of the author.


This toolbox allows MATLAB to read in and write bvh files. There are also routines for visualising the files in MATLAB.

Version 0.1 Release Notes
-------------------------

First release of the toolbox to coincide with release of C++ GPLVM code.

File Listing
------------

bvh2xyz.m: Compute XYZ values given structure and channels.
bvhConnectionMatrix.m: Compute the connection matrix for the structure.
bvhModify.m: Helper code for visualisation of bvh data.
bvhPlayData.m: Play bvh motion capture data.
bvhPlayFile.m: Play motion capture data from a bvh file.
bvhReadFile.m: Reads a bvh file into a tree structure.
bvhVisualise.m: For updating a bvh representation of 3-D data.
bvhWriteFile.m: Write a bvh file from a given structure and channels.
rotationMatrix.m: Compute the rotation matrix for an angle in each direction.
stickModify.m: Helper code for visualisation of a stick man.
stickVisualise.m: For updateing a stick representation of 3-D data.
