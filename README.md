|----------------------------------|----------------------------------------------------------|------------------------------------------------------------------|
| [Neil's homepage](/~neill/)  |   | [research group](http://intranet.cs.man.ac.uk/mlo/)  |   | [MLO Researchers](http://intranet.cs.man.ac.uk/mlo/people/)  |   |

|-----|---------------------------------------------------------------------------------------------------------|
|     | ![University Logo](http://www.cs.manchester.ac.uk/ai/pictures/icons/est1824.gif "Manchester Est. 1824") |

MATLAB Motion Capture Toolbox
=============================

The MATLAB motion capture toolbox allows loading and playing of BVH and acclaim files in MATLAB.

Warning this toolbox seems to be affected by a possible bug in MATLAB 7.4, see [here](../matlabBug.html) for details.

The MOCAP software can be downloaded [here](/neill-bin/software/downloadForm.cgi?toolbox=mocap).

Release Information
-------------------

**Current release is 0.136**.

As well as downloading the MOCAP software you need to obtain the toolboxes specified below. **These can be downloaded using the *same* password you get from registering for the MOCAP software.**

|---------------------------------------------------|-------------|
| **Toolbox**                                       | **Version** |
| [NDLUTIL](/~neill/ndlutil/downloadFiles/vrs0p161) | 0.161       |

Missing a file for reading the poser data.

#### Release 0.135

Added visualisation files written by Carl Henrik Ek for the Human Eva data and for Ankur Agarwal's Poser generated silhouette data.

#### Release 0.134

Bug fix release, a bug in bvh2xyz meant that if a position was included in the bvh skeleton structure for non-root nodes, the xyz positions were computed incorrectly. Thanks to Richard Widgery and Christopher Hulbert for identifying this problem.

horse.bvh removed due to copyright reasons. To obtain a license for this file, and plenty of other motion capture data of horses, please contact Richard Widgery of Kinetic Impulse.

#### Release 0.133

This release is a bug fix release, to deal with bugs in mocapResultsCppBvh, thanks to Cedric Vanaken for pointing out the problem.

#### Release 0.132

Moved tree handling code into NDLUTIL toolbox, version 0.156.

#### Release 0.131

This release fixes a bug where the left was displayed on the right and vice versa.

#### Release 0.13

New in this release is the ability to load data in the acclaim ASF/AMC format. See the example below for details.

Examples
--------

Once downloaded you can try loading a BVH data set from the examples directory.
`>> [skel, channels, frameLength] = bvhReadFile('examples/Swagger.bvh');  >>`

This motion capture data was taken from Ohio State University's [ACCAD](http://accad.osu.edu/research/mocap/mocap_data.htm) centre.

You can now play the data using the command

`>> skelPlayData(skel, channels, frameLength);`

You can also download and read data in the Acclaim format (asf/amc). In the example below we load a nd play data from the [CMU Graphics Lab Motion Capture Database](http://mocap.cs.cmu.edu). We use the tenth example from the 86th subject. We have assumed that you have placed the files in the `examples` subdirectory.

`>> skel = acclaimReadSkel('examples/86.asf');  >> [channels, skel] = acclaimLoadChannels('examples/86_10.amc', skel);  >> skelPlayData(skel, channels, 1/120);  >> `

Where the frames per second is given on the CMU site as 120.

Page updated on Tue Oct 6 17:26:28 2009
| [Disclaimer](http://www.manchester.ac.uk/aboutus/documents/disclaimer/ "Disclaimer") | [Privacy](http://www.manchester.ac.uk/aboutus/documents/privacy/ "Privacy") | [Copyright notice](http://www.manchester.ac.uk/aboutus/documents/copyright/ "Copyright Notice") | [Accessibility](http://www.manchester.ac.uk/aboutus/documents/accessibility/ "Accessibility") | [Freedom of information](http://www.manchester.ac.uk/aboutus/documents/foi/ "Freedom of information") | [Feedback](http://www.manchester.ac.uk/aboutus/contact/feedback/ "Feedback") |

Please contact <webmaster.cs@manchester.ac.uk> with comments and suggestions


