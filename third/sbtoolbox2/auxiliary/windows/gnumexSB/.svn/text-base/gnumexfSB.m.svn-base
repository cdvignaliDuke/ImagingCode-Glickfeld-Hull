function varargout = gnumexfSB(directories)
% gnumexfSB: an adapted version of the gnumex package
% (https://sourceforge.net/projects/gnumex). It allows everyone to build
% Fortran MEX files under Windows using the freely available MinGW 
% compiler (http://www.mingw.org/). 
%
% gnumexfSB constructs a MEX options file for windows that is used by the
% MEX function to determine the compiler, linker, options, etc.
% 
% This file does not need to be run by the user. Instead have a
% look at the mexfSB function included in this folder.

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
    error('This function only is required under Windows.');
end

% Check if MinGW is present
if exist('mingw32') ~= 7,
    error('The mingw compiler (found at: http://www.mingw.org/) is not installed or not present in the MATLAB path.');
end

% get the paths
mingwpath = shortpath(getPath('gcc.exe','mingwSB'));
optfile = 'f77.bat'; % create options file in current folder
gnumexpath = shortpath(getPath('gnumexfSB.m','gnumexSB'));
matlabpath = shortpath(matlabroot);

% write out the options file
fid = fopen(optfile,'w');
fprintf(fid,'@echo off\n');
fprintf(fid,'rem f77.bat\n');
fprintf(fid,'rem\n');
fprintf(fid,'rem    Compile and link options used for building Fortran MEX-files\n');
fprintf(fid,'rem    using the mingw compiler.\n');
fprintf(fid,'rem\n');
fprintf(fid,'rem    Options file based on gnumex generated file. \n');
fprintf(fid,'rem    Author: Henning Schmidt, henning@sbtoolbox2.org\n');
fprintf(fid,'rem    1st of January, 2008\n');
fprintf(fid,'rem\n');
fprintf(fid,'rem ********************************************************************\n');
fprintf(fid,'rem Copyright 2000 Free Software Foundation, Inc.\n');
fprintf(fid,'rem This program is free software; you can redistribute it and/or modify\n');
fprintf(fid,'rem it under the terms of the GNU General Public License as published by\n');
fprintf(fid,'rem the Free Software Foundation; either version 2 of the License, or\n');
fprintf(fid,'rem (at your option) any later version.\n');
fprintf(fid,'rem ********************************************************************\n');
fprintf(fid,'\n');
fprintf(fid,'rem ********************************************************************\n');
fprintf(fid,'rem Path parameters\n');
fprintf(fid,'rem ********************************************************************\n');
fprintf(fid,'set MATLAB=%s\n',matlabpath);
fprintf(fid,'set GM_PERLPATH=%s\\sys\\perl\\win32\\bin\\perl.exe\n',matlabpath);
fprintf(fid,'set GM_UTIL_PATH=%s\n',gnumexpath);
fprintf(fid,'set PATH=%s\\bin;%%PATH%%\n',mingwpath);
fprintf(fid,'\n');
fprintf(fid,'rem ********************************************************************\n');
fprintf(fid,'rem Compiler parameters\n');
fprintf(fid,'rem ********************************************************************\n');
fprintf(fid,'set GM_ADD_LIBS=-lg2c\n');
fprintf(fid,'set GM_MEXTYPE=mex\n');
fprintf(fid,'set GM_MEXLANG=f\n');
fprintf(fid,'set COMPILER=gcc\n');
fprintf(fid,'set COMPFLAGS=-c -DMATLAB_MEX_FILE -mrtd -fcase-upper -fno-underscoring -fleading-underscore \n');
fprintf(fid,'set OPTIMFLAGS=-O3 -malign-double -fno-exceptions -march=pentium4 -ffast-math -funroll-all-loops\n');
fprintf(fid,'set DEBUGFLAGS=-g\n');
fprintf(fid,'set NAME_OBJECT=-o\n');
fprintf(fid,'\n');
fprintf(fid,'rem ********************************************************************\n');
fprintf(fid,'rem Linker parameters\n');
fprintf(fid,'rem ********************************************************************\n');
fprintf(fid,'set LINKER=%%GM_PERLPATH%% %%GM_UTIL_PATH%%\\linkmex.pl\n');
fprintf(fid,'set LINKFLAGS=\n');
fprintf(fid,'set LINKOPTIMFLAGS=-s\n');
fprintf(fid,'set LINKDEBUGFLAGS=-g  -Wl,--image-base,0x28000000\n');
fprintf(fid,'set LINK_FILE=\n');
fprintf(fid,'set LINK_LIB=\n');
fprintf(fid,'set NAME_OUTPUT=-o %%OUTDIR%%%%MEX_NAME%%.mexw32\n');
fprintf(fid,'\n');
fprintf(fid,'rem ********************************************************************\n');
fprintf(fid,'rem Resource compiler parameters\n');
fprintf(fid,'rem ********************************************************************\n');
fprintf(fid,'set RC_COMPILER=%%GM_PERLPATH%% %%GM_UTIL_PATH%%\\rccompile.pl  -o %%OUTDIR%%mexversion.res\n');
fprintf(fid,'set RC_LINKER=\n');
fprintf(fid,'\n');
fclose(fid);

% final message
%disp('Fortran options file created');
return

function [elements] = explodePCSB(text,varargin)
if nargin == 1,
    separatorCharacter = ',';
elseif nargin == 2,
    separatorCharacter = varargin{1};
else
    error('Incorrect number of input arguments.');
end
elements = {};
openParenthesis = 0;
lastIndex = 1;
elementIndex = 1;
for k2 = 1:length(text),
    if text(k2) == '(',
        openParenthesis = openParenthesis + 1;
    end
    if text(k2) == ')',
        openParenthesis = openParenthesis - 1;
    end
    if text(k2) == separatorCharacter && openParenthesis == 0,
        elements{elementIndex} = strtrim(text(lastIndex:k2-1));
        elementIndex = elementIndex + 1;
        lastIndex = k2+1;
    end
end
elements{elementIndex} = strtrim(text(lastIndex:end));
return

% subroutines
function s = sepcat(strs, sep)
% returns cell array of strings as one char string, separated by sep
if nargin < 2
    sep = ';';
end
if isempty(strs)
    s = '';
    return
end
strs = strs(:)';
strs = [strs; repmat({sep}, 1, length(strs))];
s = [strs{1:end-1}];
return

function result = getPath(searchname,lastname)
fullpath = which(searchname);
parts = explodePCSB(fullpath,'\');
result = '';
for k=1:length(parts),
    if k>1,
        result = [result '\' parts{k}];
    else
        result = parts{k};
    end
    if strcmp(parts{k},lastname),
        % if lastname detected then return
        break;
    end
end
return
