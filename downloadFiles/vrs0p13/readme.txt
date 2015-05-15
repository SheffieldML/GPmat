MOCAP software
Version 0.13		Wednesday 22 Feb 2006 at 08:22
Copyright (c) 2006 Neil D. Lawrence

This toolbox allows MATLAB to read in and write bvh files. There are also routines for visualising the files in MATLAB.

Version 0.13 Release Notes
--------------------------

Added the ability to read motion capture files in the Acclaim format (asf and amc) to facilitate using the CMU Motion Capture database.

Version 0.12 Release Notes
--------------------------

mocapResultsCppBvh now uses the FGPLVM toolbox rather than the GPLVM toolbox.

Version 0.11 Release Notes
--------------------------

There were some missing files in the previous version that have been added in this version.

Version 0.1 Release Notes
-------------------------

First release of the toolbox to coincide with release of C++ GPLVM code.


MATLAB Files
------------

Matlab files associated with the toolbox are:

acclaim2xyz.m: Compute XYZ values given skeleton structure and channels.
acclaimLoadChannels.m: Load the channels from an AMC file.
acclaimPlayFile.m: Play motion capture data from a asf and amc file.
acclaimReadSkel.m: Reads an ASF file into a skeleton structure.
bvh2xyz.m: Compute XYZ values given structure and channels.
bvhConnectionMatrix.m: Compute the connection matrix for the structure.
bvhModify.m: Helper code for visualisation of bvh data.
bvhPlayData.m: Play bvh motion capture data.
bvhPlayFile.m: Play motion capture data from a bvh file.
bvhReadFile.m: Reads a bvh file into a tree structure.
bvhVisualise.m: For updating a bvh representation of 3-D data.
bvhWriteFile.m: Write a bvh file from a given structure and channels.
mocapResultsCppBvh.m: Load results from cpp file and visualise as a bvh format.
rotationMatrix.m: Compute the rotation matrix for an angle in each direction.
skel2xyz.m: Compute XYZ values given skeleton structure and channels.
skelConnectionMatrix.m: Compute the connection matrix for the structure.
skelModify.m: Helper code for visualisation of skel data.
skelPlayData.m: Play skel motion capture data.
skelReverseLookup.m: Return the number associated with the joint name.
skelVisualise.m: For updating a skel representation of 3-D data.
smoothAngleChannels.m: Try and remove artificial discontinuities associated with angles.
stickModify.m: Helper code for visualisation of a stick man.
stickVisualise.m: For updateing a stick representation of 3-D data.
treeFindChildren.m: Given a tree that lists only parents, add children.
treeFindParents.m: Given a tree that lists only children, add parents.
treeSwapNode.m: Swap two nodes in the tree structure array.
