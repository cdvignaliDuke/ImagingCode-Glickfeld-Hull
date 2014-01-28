function dirs = GetDirs (root,analysis)
% GETDIRS Returns directory structure based on user name and host name.
%
% based on directories.m by Mark Histed
% modified by Kenichi Ohki 7/14/2009
% $id$

switch username

    case {'kenichi'}
        dirs.tif=fullfile(root,'tif');
        dirs.rtif=fullfile(root,'rtif');
        dirs.label=fullfile(analysis,'labelimg');
        dirs.tcourse=fullfile(analysis,'tcourse');
        dirs.stats=fullfile(analysis,'stats');
        dirs.oristats=fullfile(analysis,'oristats');
        MakeDirs(dirs);

    otherwise
        dirs.tif=fullfile(root,'tif');
        dirs.rtif=fullfile(root,'rtif');
        dirs.label=fullfile(analysis,'labelimg');
        dirs.tcourse=fullfile(analysis,'tcourse');
        dirs.stats=fullfile(analysis,'stats');
        MakeDirs(dirs);
end

return;

end

function MakeDirs (dirs)
    dirlist=struct2cell(dirs);
    for i=1:length(dirlist)
        if ~exist(dirlist{i})
            mkdir(dirlist{i});
        end
    end
end
            

