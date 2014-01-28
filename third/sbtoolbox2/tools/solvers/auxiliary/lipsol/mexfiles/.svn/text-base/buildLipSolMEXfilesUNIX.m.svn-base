% buildMEXfilesUNIX: Builds the lipsol mex files on Unix systems
%
% Please note that usually on Unix the MEX files do have to be rebuild.

% Information:
% ============
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org
% GNU GPL >= 2

disp('Building lipsol MEX functions.');

cd srcunix
mex  ordmmd.f  ordmmdg.f
mex  symfct.f  symfctg.f
mex  inpnv.f   inpnvg.f
mex  blkfct.f  blkfctg.f
mex  blkslv.f  blkslvg.f

eval(sprintf('!mv ordmmd ../mex/ordmmd.%s',mexext));
eval(sprintf('!mv symfct ../mex/symfct.%s',mexext));
eval(sprintf('!mv inpnv  ../mex/inpnv.%s',mexext));
eval(sprintf('!mv blkfct ../mex/blkfct.%s',mexext));
eval(sprintf('!mv blkslv ../mex/blkslv.%s',mexext));
cd ..

