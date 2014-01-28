function [] = installSB(varargin)
% installSB
% Installation function for the SBTOOLBOX2. 
% Edit the data below to match your system and run it.
%
%       installSB
%
% This function can be called with the optional syntax:
%
%       installSB('quick')
%
% This adds the SBTOOLBOX2 and all subdirectories to the MATLAB path.
% If you are in a computer environment in which MATLAB is not able to 
% store the path settings you can run this quick installation function 
% each time when starting MATLAB. 
%
% Note that the quick installation does not build the necessary libraries.
% This is only done by calling installSB without the optional argument. Use
% this the first time when installing a new version of the toolbox.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT THE FOLLOWING VARIABLES TO MATCH YOUR SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBML TOOLBOX (NOT FOR WINDOWS USERS!!!)
% =======================================
% Add the path to the folder in which the TranslateSBML and the OutputSBML 
% MEX functions are located. Only needed for Unix/Linux/Mac users, since
% the Windows MEX functions are included in the SBTOOLBOX2.
PATH_SBMLTOOLBOX = ''; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BELOW NO MANUAL CHANGES ARE REQUIRED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check previous installations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SBTpresent = ver('SBTOOLBOX');
if length(SBTpresent) > 0,
    error('You seem to still have installed the old generation SBTOOLBOX. Please delete it before installing the SBTOOLBOX2 to avoid undersired effects.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
quickinstall = 0;
if nargin == 1,
    if strcmp(varargin{1},'quick'),
        quickinstall = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean up desktop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version of the toolbox for display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VERSION = 'Version 2.1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that installSB is started in the right folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentDir = pwd;
installSBDir = fileparts(which('installSB.m'));
if ~strcmp(currentDir,installSBDir),
    error('Run the ''installSB'' script from the folder where it is located.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that correct local path (network paths are not allowed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(currentDir(1:2),'\\'),
    error(sprintf('The installation can not be run from a network path (\\\\...).\nPlease run the installation from a local path.'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the MATLAB version 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mver = ver('matlab'); mver = str2num(mver.Version);
if mver < 7.1,
    error(sprintf('Your version of MATLAB is to old!\nPlease consider an update to at least R2006a (Version 7.2).'));
end
if ~quickinstall && mver < 7.2,
    disp(sprintf('Your version of MATLAB is older than R2006a.\nThe majority of the functions will still work but there is no guarantee that all functions do.'));
    disp('Press any key to continue ...'); 
    pause;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply path settings and save path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PATH_SBTOOLBOX2 = pwd;
warning off;
addpath(genpath(PATH_SBTOOLBOX2))
addpath(genpath(PATH_SBMLTOOLBOX))
addpath(tempdirSB)
result = savepath;
warning on;
if result == 1 && ~quickinstall,
    disp(' ');
    disp('Your MATLAB installation does not allow the saving of the updated MATLAB path variable.');
    disp('This means that you have to add the SBTOOLBOX2 folder and all its subfolders to the MATLAB');
    disp(sprintf('path each time you start MATLAB. For simplifying this task you can call installSB(''quick''),'));
    disp('which does that without rebuilding all the libraries.');
    disp(' ');
    disp('Press any key to continue.');
    disp(' ');
    pause;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile and install the needed packages 
% Compilation is done for Unix AND Windows systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~quickinstall,
    if ~ispc,
        % UNIX AND MAC
        try
            % quadprogSB
            cd(fileparts(which('buildQuadProgUNIX.m')));
            disp('PLEASE IGNORE THE WARNING MESSAGES! ... THEY ARE OK');
            buildQuadProgUNIX
            cd(PATH_SBTOOLBOX2)
            disp('Don''t worry, the last few warnings are totally ok :-)');
            disp(' ');        
        catch end
        try
            % lipsolSB
            cd(fileparts(which('buildLipSolMEXfilesUNIX.m')));
            disp('PLEASE IGNORE THE WARNING MESSAGES! ... THEY ARE OK');
            buildLipSolMEXfilesUNIX
            cd(PATH_SBTOOLBOX2)
            disp('Don''t worry, the last few warnings are totally ok :-)');
            disp(' ');        
        catch end
        try
            % gplkSB
            olddir = pwd();
            cd(fileparts(which('buildGLPKMEXunix.m')));
            buildGLPKMEXunix
            cd(PATH_SBTOOLBOX2)
        catch end
    else
        % PC - WINDOWS
        olddir = pwd();
        % shortpath
        cd(fileparts(which('shortpath.c')));
        mex shortpath.c
        cd(PATH_SBTOOLBOX2)
        if useMingwSB(),
            if mver > 7.1,
                % only build for versions > 7.1
                try
                    % quadprog
                    cd(fileparts(which('buildQuadProgWindows.m')));
                    disp('PLEASE IGNORE THE WARNING MESSAGES! ... THEY ARE OK');
                    buildQuadProgWindows
                    cd(PATH_SBTOOLBOX2)
                    clc;
                    disp('Don''t worry, the last few warnings are totally ok :-)');
                    disp(' ');
                catch end
                try
                    % lipsolSB
                    cd(fileparts(which('buildLipSolMEXfilesWINDOWS.m')));
                    disp('PLEASE IGNORE THE WARNING MESSAGES! ... THEY ARE OK');
                    buildLipSolMEXfilesWINDOWS
                    cd(PATH_SBTOOLBOX2)
                    clc;
                    disp('Don''t worry, the last few warnings are totally ok :-)');
                    disp(' ');
                catch end
            end
            try
                % gplkSB
                olddir = pwd();
                cd(fileparts(which('buildGLPKMEXwindows.m')));
                buildGLPKMEXwindows
                cd(PATH_SBTOOLBOX2)
            catch end
            try
                % TranslateSBML / OutputSBML
                cd(fileparts(which('buildSBMLFCTSwindows.m')));
                buildSBMLFCTSwindows
                cd(PATH_SBTOOLBOX2)            
            catch end
        else
            % if mingw is not used then do not recompile the lipsol and
            % glpk functions
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output license information, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
disp(sprintf('Systems Biology Toolbox 2 for MATLAB - %s\n',VERSION));
disp(sprintf('Developed by Henning Schmidt, henning@sbtoolbox2.org.'));
disp(sprintf('The SBTOOLBOX2 contains several third party packages, therefor'));
disp(sprintf('copyright statements are present in the individual functions.'));
disp(' ');
disp(sprintf('This is free software, distributed under the GNU General Public License.'));
disp(sprintf('You are welcome to redistribute it under certain conditions; type'));
disp(sprintf('''open copying.html'' for details. The toolbox comes with ABSOLUTELY NO'));
disp(sprintf('WARRANTY; for details type ''open warranty.html''. The complete GNU GPL'));
disp(sprintf('can be found at http://www.gnu.org/licenses/gpl.html.'));
disp(sprintf('The newest version can always be downloaded from http://www.sbtoolbox2.org\n'));
if quickinstall,
    disp(sprintf('installSB quickinstall: Note that the quick installation function does\nnot build the necessary libraries. This is only done by installSB()\nwhich should be used when installing SBTOOLBOX2 the first time.\n(WINDOWS: libraries are only rebuild if MinGW is the chosen compiler).'));
end
if ispc,
    if useMingwSB,
        disp(sprintf('\nMEX compiler: MinGW\n'));
    else
        disp(sprintf('\nMEX compiler: LCC\n'));
    end
end