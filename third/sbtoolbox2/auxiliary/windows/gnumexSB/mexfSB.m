function [] = mexfSB(options)
% mexfSB: compile Fortran MEX files under Windows using the MinGW compiler.
% The MinGW compiler is already included in the SBTOOLBOX2 distribution.
% 
% Then you just go to the folder in which you have your Fortran MEX source
% files you want to compile and type:
%
% >> mexfSB(optionsstring)
%
% optionsstring: all the stuff you would write as options to the standard
%       mex call. Usually you will only specify the fortran files to be
%       compiled.
%
% Example 1:      mexfSB('ordmmd.f ordmmdg.f')
% Example 2:      In the folder SBTOOLBOX/tools/solvers/auxiliary/lipsolSB/mexfiles 
%                 you can have a look at the script:
%                 buildLipSolMEXfilesWINDOWS.m to see a working example.

% Information:
% ============
% Copyright (C) 2008 Henning Schmidt, henning@sbtoolbox2.org
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
% USA.

% under vista the mingw compiler does not work
if ~useMingwSB(),
    error(sprintf('Sorry, but using the LCC compiler you can not compile C++ MEX functions.'));
end

% Check if Windows. Otherwise error
if isunix,
    error('This function makes only sense on Windows computers.');
end

% Create the f77.bat file using the gnumexfSB function
gnumexfSB();

% Construct string with the paths to the libraries that need to be linked
shortmatlabroot = [shortpath(matlabroot) '\extern\lib\win32\microsoft'];
libs = sprintf('"%s\\libmat.lib" "%s\\libmex.lib" "%s\\libmx.lib"',shortmatlabroot,shortmatlabroot,shortmatlabroot);

% Constructing the mex call
mexcallcmd = sprintf('mex -f f77.bat %s %s',options,libs);

% Run mex
eval(mexcallcmd);

% Delete the options file again
delete f77.bat
