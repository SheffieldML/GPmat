MOCAP software
Version 0.135		Tuesday 21 Oct 2008 at 10:52

This toolbox allows MATLAB to read in and write bvh files and read acclaim files. There are also routines for visualising the files in MATLAB.

Version 0.135 Release Notes
---------------------------

Added visualisation files written by Carl Henrik Ek for the Human Eva data and for Ankur Agarwal's Poser generated silhouette data.

Version 0.134 Release Notes
---------------------------

Bug fix release, a bug in bvh2xyz meant that if a position was included in the bvh skeleton structure for non-root nodes, the xyz positions were computed incorrectly. Thanks to Richard Widgery and Christopher Hulbert for identifying this problem. 

horse.bvh removed due to copyright reasons. To obtain a license for this file, and plenty of other motion capture data of horses, please contact Richard Widgery of Kinetic Impulse.

Version 0.133 Release Notes
---------------------------

Bug fix release, to deal with bugs in mocapResultsCppBvh, thanks to Cedric Vanaken for pointing out the problem.

Version 0.132 Release Notes
---------------------------

Moved tree handling code into NDLUTIL toolbox, version 0.156.

Version 0.131 Release Notes
---------------------------

The axis display for visualising the data gave the limbs on the wrong side (this is a problem of definition of the axis system). It was fixed by placing 

axis ij 

in the skeVisualise.m file. Thanks to Heike Vallery for pointing out this error.

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

xyzhumanevaJoint2pos.m:
bvhModify.m: Helper code for visualisation of bvh data.
xyzankurVisualise.m:
xyzhumanevaDraw.m:
stickVisualise.m: For drawing a stick representation of 3-D data.
rotationMatrix.m: Compute the rotation matrix for an angle in each direction.
bvhReadFile.m: Reads a bvh file into a tree structure.
xyzhumanevaVisualise3d.m:
skelPlayData.m: Play skel motion capture data.
xyzhumanevaHeadingAngle.m: 
xyzhumanevaAnim.m:
mocapToolboxes.m: Toolboxes required by the MOCAP toolbox.
stickModify.m: Helper code for visualisation of a stick man.
xyzankurVisualise2.m:
acclaimReadSkel.m: Reads an ASF file into a skeleton structure.
xyzhumanevaGenerateMovie.m:
acclaimPlayFile.m: Play motion capture data from a asf and amc file.
bvhVisualise.m: For updating a bvh representation of 3-D data.
skel2xyz.m: Compute XYZ values given skeleton structure and channels.
bvhPlayFile.m: Play motion capture data from a bvh file.
skelReverseLookup.m: Return the number associated with the joint name.
xyzhumanevaModify2.m:
xyzhumanevaVisualise2.m:
skelVisualise.m: For drawing a skel representation of 3-D data.
skelModify.m: Update visualisation of skeleton data.
bvhWriteFile.m: Write a bvh file from a given structure and channels.
skelConnectionMatrix.m: Compute the connection matrix for the structure.
xyzhumanevaVisualise.m:
mocapResultsCppBvh.m: Load results from cpp file and visualise as a bvh format.
xyzhumanevaRemovePart.m:
bvhConnectionMatrix.m: Compute the connection matrix for the structure.
xyzankurAnim.m:
xyzankurModify.m:
xyzhumanevaModify.m:
xyzankur2joint.m:
xyzankurDraw.m:
acclaimLoadChannels.m: Load the channels from an AMC file.
bvhPlayData.m: Play bvh motion capture data.
xyzhumaneva2joint.m:
smoothAngleChannels.m: Try and remove artificial discontinuities associated with angles.
bvh2xyz.m: Compute XYZ values given structure and channels.
acclaim2xyz.m: Compute XYZ values given skeleton structure and channels.
xyzhumanevaVisualiseModes.m:
xyzhumanevaAlign.m:
