@echo off
rem f77.bat
rem
rem    Compile and link options used for building Fortran MEX-files
rem    using the mingw compiler.
rem
rem    Options file based on gnumex generated file. 
rem    Author: Henning Schmidt, henning@sbtoolbox2.org
rem    1st of January, 2008
rem
rem ********************************************************************
rem Copyright 2000 Free Software Foundation, Inc.
rem This program is free software; you can redistribute it and/or modify
rem it under the terms of the GNU General Public License as published by
rem the Free Software Foundation; either version 2 of the License, or
rem (at your option) any later version.
rem ********************************************************************

rem ********************************************************************
rem Path parameters
rem ********************************************************************
set MATLAB=C:\PROGRA~1\MATLAB\R2006a
set GM_PERLPATH=C:\PROGRA~1\MATLAB\R2006a\sys\perl\win32\bin\perl.exe
set GM_UTIL_PATH=E:\WORK\TOOLBO~1\WORKIN~1\SBTOOL~1\AUXILI~1\windows\gnumexSB
set PATH=E:\WORK\TOOLBO~1\WORKIN~1\SBTOOL~1\AUXILI~1\windows\mingwSB\bin;%PATH%

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
set GM_ADD_LIBS=-lg2c
set GM_MEXTYPE=mex
set GM_MEXLANG=f
set COMPILER=gcc
set COMPFLAGS=-c -DMATLAB_MEX_FILE -mrtd -fcase-upper -fno-underscoring -fleading-underscore 
set OPTIMFLAGS=-O3 -malign-double -fno-exceptions -march=pentium4 -ffast-math -funroll-all-loops
set DEBUGFLAGS=-g
set NAME_OBJECT=-o

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
set LINKER=%GM_PERLPATH% %GM_UTIL_PATH%\linkmex.pl
set LINKFLAGS=
set LINKOPTIMFLAGS=-s
set LINKDEBUGFLAGS=-g  -Wl,--image-base,0x28000000\n
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT=-o %OUTDIR%%MEX_NAME%.mexw32

rem ********************************************************************
rem Resource compiler parameters
rem ********************************************************************
set RC_COMPILER=%GM_PERLPATH% %GM_UTIL_PATH%\rccompile.pl  -o %OUTDIR%mexversion.res
set RC_LINKER=

