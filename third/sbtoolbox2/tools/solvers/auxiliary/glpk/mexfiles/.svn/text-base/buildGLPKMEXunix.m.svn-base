% This script generates a static library file for the GNU GLPK library
% for Unix systems

% Information:
% ============
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org
% GNU GPL >= 2

disp('Building GLPK MEX function.');

% Build the library
cd glpk-4.17/src
mex -c -O -I../include *.c
!ar rc glpklib.a *.o
delete *.o
!mv *.a ../..
cd ../..

% Compile the MEX file
mex -O *.cpp -Iglpk-4.17/include glpklib.a
delete *.a