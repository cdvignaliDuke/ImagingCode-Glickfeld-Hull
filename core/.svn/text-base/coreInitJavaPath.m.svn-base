function coreInitJavaPath(svnroot, ijroot)
%COREINITJAVAPATH Initializes java path for Reid Lab Core
% COREINITJAVAPATH(SVNROOT, IJROOT) where is SVNROOT is the Reid Lab SVN
% repository root directory and IJROOT and the ImageJ root directory. Due
% to an unresolved MATLAB bug (R14), JAVAINITPATH function clears the global variables.
% $Id: coreInitJavaPath.m 204 2008-04-30 19:07:26Z histed $

svncore = fullfile(svnroot,'core');

if nargin < 1, error('must specify location of reidlab/core directory'); end
if nargin < 2, ijroot = 'D:\Program Files (x86)\ImageJ'; end

%% java path (WARNING JAVAADDPATH CLEARS GLOBAL VARIABLES)
javaaddpath(ijroot);
javaaddpath(fullfile(ijroot,'ij.jar'));
javaaddpath(fullfile(ijroot,'plugins','Stacks'));
javaaddpath(fullfile(ijroot,'plugins','Hypervolumes'));
javaaddpath(fullfile(ijroot,'plugins','turboreg'));
javaaddpath(fullfile(svncore, 'java'));

% some required code might sit there
% javaaddpath('E:\Shared_Java');

disp('COREINITJAVAPATH: Added java classpath dirs (global vars cleared!!)');

return;
