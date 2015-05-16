
Matlab Noise Models Toolbox
===========================


Release Information
-------------------

**Current release is 0.141**.

Added noiseReadFromFile and noiseReadFromFID for compatability with CPP releases.

#### Version 0.14

As of version 0.14, the NCNM noise model has been incorporated into this toolbox.

This toolbox allows computation of several different noise models and their gradients. You can add noise models to the toolbox by creating versions of the relevant files. Once added, they can be tested using the `noiseTest`. For example you can test the probit noise model by writing

```
>> noiseTest('probit')
```

There are several noise models implemented, the mains ones that are being maintained for the latest release are:

`gaussian  ordered  probit  ncnm`

These noise models are Gaussian (standard Gaussian noise), probistic regression for binary classificaiton and ordered categorical for ordered categories (ordinal regression).

Page updated on Tue Oct 6 17:28:55 2009


