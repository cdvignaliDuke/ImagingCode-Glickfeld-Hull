% buildQuadProgWindows: Builds the quadprog mex files on Windows systems
% In order to build the MEX files on Windows you need to install the MinGW
% compiler.
%
% Please note that usually on Windows the MEX files do not have to be rebuild.

% Information:
% ============
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org
% GNU GPL >= 2

disp('Building quadprog MEX functions.');

% build the MEX files using the mexfSB function included in the SBTOOLBOX2
cd src
mexfSB('aind.f aindg.f');
mexfSB('solveQP.f solveQPg.f util.f');
mexfSB('solveQPcompact.f solveQPcompactg.f util.f');

% clean up
!move *.mexw32 ../mex/
cd ..
