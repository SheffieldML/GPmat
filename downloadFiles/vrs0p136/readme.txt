MOCAP software
Version 0.136		Tuesday 06 Oct 2009 at 17:26

This toolbox allows MATLAB to read in and write bvh files and read acclaim files. There are also routines for visualising the files in MATLAB.

Version 0.136 Release Notes
---------------------------

Missing a file for reading the poser data.

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

bvhPlayData.m: Play bvh motion capture data.
xyzpoppeDraw.m: Helper function for drawing data from Poppe.
xyzhumanevaVisualise3d.m:
xyzankurModify.m:  Helper function for modifying the point cloud from Agarwal and Triggs data.
smoothAngleChannels.m: Try and remove artificial discontinuities associated with angles.
bvhConnectionMatrix.m: Compute the connection matrix for the structure.
acclaimPlayFile.m: Play motion capture data from a asf and amc file.
xyzhumanevaJoint2pos.m:
xyzhumanevaAnim.m:
acclaimGradient.m: computes the gradient of x,y,z locations wrt angles.
xyzankurAnim.m: Animate point cloud of stick man from Agarwal & Triggs dataset.
bvhReadFile.m: Reads a bvh file into a tree structure.
xyzhumanevaAlign.m:
xyzpoppeVisualise.m: Draw the Poppe figure return the graphics handle.
bvhModify.m: Helper code for visualisation of bvh data.
xyzankurDraw.m: Helper function for drawing the point cloud from Agarwal and Triggs data.
xyzankurVisualise2.m:
xyzhumanevaModify.m:
xyzhumanevaDraw.m:
bvhPlayFile.m: Play motion capture data from a bvh file.
skelReverseLookup.m: Return the number associated with the joint name.
acclaimLoadChannels.m: Load the channels from an AMC file.
mocapToolboxes.m: Toolboxes required by the MOCAP toolbox.
acclaim2xyz.m: Compute XYZ values given skeleton structure and channels.
xyzankur2joint.m: Converts data to xyz positions for each joint.
xyzankurAnimCompare.m: Animate a prediction and ground truth for stick man from Agarwal & Triggs dataset.
skelVisualise.m: For drawing a skel representation of 3-D data.
skelConnectionMatrix.m: Compute the connection matrix for the structure.
xyzhumanevaVisualiseModes.m:
stickModify.m: Helper code for visualisation of a stick man.
rotationMatrix.m: Compute the rotation matrix for an angle in each direction.
xyzhumanevaVisualise2.m:
rotationMatrixGradient.m: Compute the gradient of rotation with respect to one angle.
skel2xyz.m: Compute XYZ values given skeleton structure and channels.
mocapResultsCppBvh.m: Load results from cpp file and visualise as a bvh format.
xyzpoppe2joint.m:
xyzpoppeAnim.m: Animate point cloud of stick man from Poppe dataset.
skelPlayData.m: Play skel motion capture data.
bvhVisualise.m: For updating a bvh representation of 3-D data.
xyzpoppeModify.m:
xyzhumanevaModify2.m:
bvh2xyz.m: Compute XYZ values given structure and channels.
xyzhumanevaHeadingAngle.m: 
xyzankurVisualise.m: Draw the Agarwal & Triggs figure return the graphics handle.
xyzhumaneva2joint.m:
acclaimReadSkel.m: Reads an ASF file into a skeleton structure.
xyzhumanevaGenerateMovie.m:
xyzhumanevaRemovePart.m:
bvhWriteFile.m: Write a bvh file from a given structure and channels.
skelModify.m: Update visualisation of skeleton data.
stickVisualise.m: For drawing a stick representation of 3-D data.
xyzhumanevaVisualise.m:
