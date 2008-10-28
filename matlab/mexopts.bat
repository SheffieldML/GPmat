@echo off
rem CVF61OPTS.BAT
rem
rem    Compile and link options used for building MEX-files
rem    using the Compaq Visual Fortran compiler version 6.1
rem
rem    $Revision $  $Date: 2006/01/13 20:56:10 $
rem
rem ********************************************************************
rem General parameters
rem ********************************************************************
set MATLAB=c:\Archiv~1\Matlab\R2006a
set DF_ROOT=C:\Program Files\Microsoft Visual Studio
set VCDir=%DF_ROOT%\VC98
set MSDevDir=%DF_ROOT%\Common\msdev98
set DFDir=%DF_ROOT%\DF98
set PATH=%MSDevDir%\bin;%DFDir%\BIN;%VCDir%\BIN;%PATH%
set INCLUDE=%DFDir%\INCLUDE;%INCLUDE%
set LIB=%DFDir%\LIB;%VCDir%\LIB;%LIB%

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
set COMPILER=df
set COMPFLAGS=/fpp:"/m /S%MATLAB%\extern\include" -c -G5 -nologo -DMATLAB_MEX_FILE
set OPTIMFLAGS=/MD -Ox -DNDEBUG
set DEBUGFLAGS=/MD -Zi
set NAME_OBJECT=/Fo

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
set LIBLOC=%MATLAB%\extern\lib\win32\microsoft
set LINKER=link
set LINKFLAGS=/DLL /EXPORT:_MEXFUNCTION@16 /LIBPATH:"%LIBLOC%" libmx.lib libmex.lib libmat.lib /implib:%LIB_NAME%.lib /NOLOGO kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /MACHINE:IX86
set LINKOPTIMFLAGS=
set LINKDEBUGFLAGS=/debug
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT="/out:%OUTDIR%%MEX_NAME%%MEX_EXT%"
set RSP_FILE_INDICATOR=@

rem ********************************************************************
rem Resource compiler parameters
rem ********************************************************************
set RC_COMPILER=rc /fo "%OUTDIR%mexversion.res"
set RC_LINKER=
