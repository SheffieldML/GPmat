#
# g95opts.sh	Shell script for configuring MEX-file creation script,
#               mex.  These options were tested with the specified compiler.
#
# usage:        Do not call this file directly; it is sourced by the
#               mex shell script.  Modify only if you don't like the
#               defaults after running mex.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# SELECTION_TAGs occur in template option files and are used by MATLAB
# tools, such as mex and mbuild, to determine the purpose of the contents
# of an option file. These tags are only interpreted when preceded by '#'
# and followed by ':'.
#
#SELECTION_TAG_MEX_OPT: Template Options file for building Fortran 90 MEX-files via the system ANSI compiler
#
# Copyright 1984-2004 The MathWorks, Inc.
# $Revision: 1.4.4.6 $  $Date: 2004/04/25 21:30:51 $
#----------------------------------------------------------------------------
#
    TMW_ROOT="$MATLA"B
    MFLAGS=''
    if [ "$ENTRYPOINT" = "mexLibrary" ]; then
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat -lmwservices -lut"
    else  
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat"
    fi
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
    G95LIBDIR='/opt/g95-install/lib/gcc-lib/i686-pc-linux-gnu/4.1.2/'
    FC='g95'
    FFLAGS='-ffree-form'
    FLIBS="$MLIBS -lm -L$G95LIBDIR -lf95 -lcl -lc -lisamstub"
    FOPTIMFLAGS='-O'
    FDEBUGFLAGS='-g'
    
    LD="$COMPILER"
    LDEXTENSION='.mexglx'
    LDFLAGS="-pthread -shared -m32 -Wl,--version-script,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
    LDOPTIMFLAGS='-O'
    LDDEBUGFLAGS='-g'
#
    POSTLINK_CMDS=':'
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           FC="$FC"
#           FFLAGS="$FFLAGS"
#           FOPTIMFLAGS="$FOPTIMFLAGS"
#           FDEBUGFLAGS="$FDEBUGFLAGS"
#           FLIBS="$FLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#----------------------------------------------------------------------------
#############################################################################


