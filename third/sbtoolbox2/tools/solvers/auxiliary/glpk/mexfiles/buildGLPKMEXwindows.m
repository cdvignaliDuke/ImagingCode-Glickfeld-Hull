% This script generates a static library file for the GNU GLPK library
% for Windows systems

% Information:
% ============
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org
% GNU GPL >= 2

disp('Building GLPK MEX function.');

% Build the library
cd glpk-4.17/src
mexcSB('-c *.c -I../include');
% Strange construct but not all systems find ar.exe ... strange!
eval(sprintf('!%s rc glpklib.lib *.obj',shortpath(which('mingwSB\bin\ar.exe'))));
delete *.obj
!move *.lib ../..
cd ../..

% Compile the MEX file
mexcSB('*.cpp -Iglpk-4.17/include glpklib.lib');
delete *.lib