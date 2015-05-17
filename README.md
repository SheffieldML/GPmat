
Matlab Prior Toolbox
====================

This page describes the prior toolbox.

Release Information
-------------------

**Current release is 0.22**.

As well as downloading the PRIOR software you need to obtain the toolboxes specified below. 

| **Toolbox**                               | **Version** |
|-------------------------------------------|-------------|
| [OPTIMI](/optimi/downloadFiles/vrs0p132)  | 0.132       |
| [NDLUTIL](/ndlutil/downloadFiles/vrs0p16) | 0.16        |

Release Notes
-------------

### Current Release

Changes to the code for reading in priors written by C++ code.

#### Release 0.13

This toolbox allows computation of several different prior distributions and their gradients. You can add distributions to the toolbox by creating versions of the relevant files. Once added, they can be tested using the `priorTest`. For example you can test the Gaussian prior

```matlab
>> priorTest('gaussian')
```

There are several prior models implemented, the mains ones that are being maintained for the latest release are:

`gaussian  gamma` These prior models are Gaussian (zero mean with precision 1), gamma (parameterised by a and b).

Page updated on Fri Apr 17 17:43:23 2009
