function dirs = frGetDirs
% FRGETDIRS Returns directory structure based on user name and host name.
%
% by Vincent Bonin based on directories.m by Mark Histed
% $id$

global USER_PROFILE;

if isempty(USER_PROFILE)
    user = username;
else
    user = USER_PROFILE;
end;

switch user
    case {'vincent'}
        dirs.root = 'i:\users\vincent';
        dirs.svn = fullfile(dirs.root,'svn');
        dirs.data = 'j:\data';
        dirs.images = 'k:\images';
        dirs.registered = 'k:\registered';
        dirs.pca = 'k:\pca';
        dirs.analysis = 'i:\users\vincent\analysis';
        dirs.multicore = fullfile(dirs.root,'multicore');
        dirs.stimulation = fullfile(dirs.root,'stimulation');
        dirs.presentation = fullfile(dirs.stimulation,'presentation');
        dirs.plexon = fullfile(dirs.root,'plexon');
        dirs.figsOut = 'i:\figures';
        dirs.backup = '\\zoloto\bigstorlab\fastrig';
       
    case {'histed'}
        dirs.root = 'i:\users\vincent';
        dirs.svn = fullfile(dirs.root,'svn');
        dirs.data = 'j:\data';
        %dirs.data = '//zmey/storlab2/fastrig2';
        dirs.images = 'k:\images';
        dirs.pca = 'k:\pca';
        dirs.registered = 'k:\registered';
        dirs.analysis = 'i:\analysis';
        dirs.multicore = fullfile(dirs.root,'multicore'); 
        dirs.stimulation = fullfile(dirs.root,'stimulation');
        dirs.presentation = fullfile(dirs.stimulation,'presentation');
        dirs.scrapbook = fullfile(dirs.root,'scrapbooks');
        dirs.plexon = fullfile(dirs.root, 'plexon');
        dirs.figsOut = 'i:/users/histed/output-images';
        
    case {'kenichi'}
        dirs.root = 'i:\users\vincent';
        dirs.svn = fullfile(dirs.root,'svn');
        %dirs.data = 'j:\data';
        dirs.data = 'J:\data';
%        dirs.data = '//zmey/storlab2/fastrig/data';
        dirs.images = 'k:\images';
        dirs.registered = 'k:\registered';
        dirs.analysis = 'i:\analysis';
        dirs.multicore = fullfile(dirs.root,'multicore'); 
        dirs.stimulation = fullfile(dirs.root,'stimulation');
        dirs.presentation = fullfile(dirs.stimulation,'presentation');
        dirs.scrapbook = fullfile(dirs.root,'scrapbooks');
        dirs.plexon = fullfile(dirs.root, 'plexon');
        dirs.figsOut = 'i:/users/histed/output-images';

    case {'lindsey'}
        dirs.root = 'G:\users\lindsey'; % replace with your own home directory
        dirs.svn = fullfile(dirs.root,'svn');
        dirs.data = 'G:\users\lindsey\dataLG';
        dirs.images = 'G:\users\lindsey\dataLG'; % location of the raw data
        %dirs.images = '\\Zmey\storlab2\data\Lindsey\fastrig\';
        dirs.registered = 'G:\users\lindsey\regLG'; % location of the registered data
        dirs.pca = 'G:\users\lindsey\pcaLG';  
        dirs.analysis = 'G:\users\lindsey\analLG'; % location of derived data
        %dirs.multicore = fullfile(dirs.root,'multicore');
        dirs.stimulation = fullfile(dirs.root,'stimulation');
        dirs.presentation = fullfile(dirs.stimulation,'presentation');
        %dirs.plexon = fullfile(dirs.root,'plexon');
        %dirs.figsOut = 'i:\figures';
        dirs.backup = '\\zoloto\bigstorlab\fastrig';

    otherwise
        error('unknown user name');
end

return;
