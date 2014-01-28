% buildMEXfilesWINDOWS: Builds the lipsol mex files on Windows systems
% In order to build the MEX files on Windows you need to install the MinGW
% compiler.
%
% Please note that usually on Windows the MEX files do not have to be rebuild.

% Information:
% ============
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org
% GNU GPL >= 2

disp('Building lipsol MEX functions.');

% build the MEX files using the mexfSB function included in the SBTOOLBOX2
cd srcwindows
mexfSB('ordmmd.f ordmmdg.f');
mexfSB('symfct.f  symfctg.f');
mexfSB('inpnv.f   inpnvg.f');
mexfSB('blkfct.f  blkfctg.f');
mexfSB('blkslv.f  blkslvg.f');

% clean up
!move *.mexw32 ../mex/
cd ..

