This toolbox allows MATLAB to read in and write bvh files. There are also routines for visualising the files in MATLAB.

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
