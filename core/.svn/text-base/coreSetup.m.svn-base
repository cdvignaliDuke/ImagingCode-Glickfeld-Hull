function coreSetup(svnroot,ijroot,user)
%CORESETUP Set java and matlab paths, initializes global variables. 
%   SETUPCORE(SVNROOT,IJROOT,USERNAME) specifies the Reid Lab SVN root
%   directory, the ImageJ root directory, and the user name for user-specific
%   settings. Call SETUPCORE from your startup.m file to setup your
%   enviromnent.
%
%   Example:  (put in startup.m)
%       svnroot = 'i:/users/histed/reidlab-svn';
%       ijdir = 'D:\Program Files (x86)\ImageJ'; 
%       % if java plugins in different dir, use this and add manually
%       warning('off', 'MATLAB:javaclasspath:invalidFile');
%
%       cd(fullfile(svnroot, 'core'));  % must cd to find coreSetup.m
%       coreSetup(svnroot, ijdir);
%       cd('../users/yourname');        % set working dir for yourself
%
%$Id: coreSetup.m 204 2008-04-30 19:07:26Z histed $


if nargin < 1, error('must specify location of Reid Lab SVN root directory'); end
if nargin < 2, ijroot = 'D:\Program Files (x86)\ImageJ'; end
if nargin < 3, user = '' ; end;

[pFirst pLast] = fileparts(svnroot);
if strcmp(lower(pLast), 'core')
    error('svnroot must point to directory above core (svn top dir)');
end

% common to all users
coreInitJavaPath(svnroot, ijroot); % MUST BE CALLED FIRST. CLEARS GLOBAL VARIABLES
coreInitMatlabPath(svnroot, ijroot);

% remove warnings for casts permanently; used by kinds of calcium/fastrig code
warning('off', 'MATLAB:intConvertNonIntVal'); 

% user-specific settings
switch (user)
    case 'vincent'
        setupvincent(svnroot,ijroot);
end

return;

function setupvincent(svnroot,ijroot);

global DIRS;
DIRS.root = 'i:\users\vincent';
DIRS.svn = fullfile(DIRS.root,'svn');
DIRS.data = fullfile(DIRS.root,'data');
DIRS.images = fullfile(DIRS.root,'images');
DIRS.analysis = fullfile(DIRS.root,'analysis');
DIRS.multicore = fullfile(DIRS.root,'multicore'); 

disp('SETUPCORE: Initialized global variables.');

return;
