function [] = mexcSB(options)
% mexcSB: compile C/C++ MEX files under Windows using the MinGW compiler.
% The MinGW compiler is already included in the SBTOOLBOX2 distribution.
% 
% Then you just go to the folder in which you have your C/C++ MEX source
% files you want to compile and type:
%
% >> mexcSB(optionsstring)
%
% optionsstring: all the stuff you would write as options to the standard
%       mex call. Usually you will only specify the C files to be
%       compiled.
%
% Example:      mexcSB('test.c')

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

% Check if Windows. Otherwise error
if isunix,
    error('This function makes only sense on Windows computers.');
end

if useMingwSB(),
    % MinGW is used
    % Create the gccopt.bat file using the gnumexcSB function
    gnumexcSB();
    % Construct string with the paths to the libraries that need to be linked
    shortmatlabroot = [shortpath(matlabroot) '\extern\lib\win32\microsoft'];
    libs = sprintf('"%s\\libmat.lib" "%s\\libmex.lib" "%s\\libmx.lib"',shortmatlabroot,shortmatlabroot,shortmatlabroot);
    % Constructing the mex call
    mexcallcmd = sprintf('mex -f gccopt.bat %s %s',options,libs);
    % Run mex
    eval(mexcallcmd);
    % Delete the options file again
    delete gccopt.bat
else
    % LCC is used
    % Check if C++ files are to be compiled ... then error
    if ~isempty(strfind(lower(options),'.cpp')),
        error(sprintf('Sorry, but using the LCC compiler you can not compile C++ MEX functions.'));
    end
    % Run the rest
    mexcallcmd = sprintf('mex %s',options);
    % Run mex
    eval(mexcallcmd);    
end
