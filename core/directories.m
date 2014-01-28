function dirs = directories
%DIRECTORIES 
%
%   This is an example!!!  Edit this yourself to customize directory locations.
%
%$Id: directories.m 146 2008-03-26 16:16:37Z histed $

error(['This copy of directories.m is an example - edit it yourself and put ' ...
       'it on the path, ideally in your matlab dir (same dir as startup.m)']);

hName = hostname;

dirs.dataRoot = 'i:/users/vincent/data';

switch hName
  case 'mambo'
    dirs.toolsSnapshots = '/Users/histed/shared/snapshots';
  case 'zquad'
    dirs.toolsSnapshots = 'i:/users/histed/snapshots';
end



