function fnames=filenameManager(baseFilename, prefix)

% give filenames
% prefix: default='' (do not change). If you want to add 'r' '+r'. 
%   If you want to delete 'r', '-r'
%
%   5.13.2008. Kenichi Ohki

if nargin < 2
    prefix ='';
end

if ~isempty(prefix)
    if prefix(1)=='+'
        prefix(1)=[];
        baseFilename=[prefix,baseFilename];
    else if prefix(1)=='-'
            prefix(1)=[];
            if ~strcmpi(prefix, baseFilename(1:length(prefix)))
                display('prefix does not match');
                fnames='';
            end
            baseFilename(1:length(prefix))=[];
        else 
            display('first character should be + or -');
            fnames='';
        end
    end
end

fnames.data=[baseFilename,'.tif'];
fnames.mat=[baseFilename,'.mat'];
fnames.avg=[baseFilename,'_avg.tif'];
fnames.target=[baseFilename,'_target.tif'];
fnames.labelimg=[baseFilename,'_labelimg.mat'];
fnames.shadecells=[baseFilename,'_shadecells.tif'];
fnames.tcourse=[baseFilename,'_tcourse.mat'];
fnames.stats=[baseFilename,'_stats.mat'];
fnames.Oristats=[baseFilename,'_Oristats.mat'];
