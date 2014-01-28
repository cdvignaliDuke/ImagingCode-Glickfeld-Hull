function fnames=GetFnames(fname, root_dir, analysis_dir, prefix, suffix)

% give filenames
% prefix: default='' (do not change). If you want to add 'r' '+r'. 
%   If you want to delete 'r', '-r'
% suffix = '', or '_tcourse' etc., which will be removed from the file
% names
%   5.13.2008. Kenichi Ohki

if nargin < 4
    prefix ='';
end

if nargin < 5
    suffix ='';
end

dirs=GetDirs(root_dir, analysis_dir);
fname=GetBaseFname(fname, prefix, suffix);

if isempty(fname)
    fnames='';
    return;
end

switch username

    case {'kenichi'}

fnames.data=fullfile(dirs.tif,[fname,'.tif']);  %original stack
fnames.mat=fullfile(dirs.tif,[fname,'.mat']);  %mat file which contains scan parameters
fnames.rdata=fullfile(dirs.rtif,[fname,'.tif']);    %realigned stack
fnames.rparams=fullfile(dirs.rtif,[fname,'_rparams.mat']); %alignemnt parameters
fnames.avg=fullfile(dirs.rtif,[fname,'_avg.tif']);  % avg img after alignment 
fnames.target=fullfile(dirs.rtif,[fname,'_target.tif']);    %target img used for alignment
fnames.labelimg=fullfile(dirs.label,[fname,'_labelimg.mat']);   %cell mask (mat file)
fnames.shadecells=fullfile(dirs.label,[fname,'_shadecells.tif']);   %cell mask displayed on avg img
fnames.tcourse=fullfile(dirs.tcourse,[fname,'_tcourse.mat']);   %time courses of cells
fnames.stats=fullfile(dirs.stats,[fname,'_stats.mat']); %basic statistics
fnames.Oristats=fullfile(dirs.oristats,[fname,'_Oristats.mat']); %orientation statistics

    otherwise

fnames.data=fullfile(dirs.tif,[fname,'.tif']);  %original stack
fnames.mat=fullfile(dirs.tif,[fname,'.mat']);  %mat file which contains scan parameters
fnames.rdata=fullfile(dirs.rtif,[fname,'.tif']);    %realigned stack
fnames.rparams=fullfile(dirs.rtif,[fname,'_rparams.mat']); %alignemnt parameters
fnames.avg=fullfile(dirs.rtif,[fname,'_avg.tif']);  % avg img after alignment 
fnames.target=fullfile(dirs.rtif,[fname,'_target.tif']);    %target img used for alignment
fnames.labelimg=fullfile(dirs.label,[fname,'_labelimg.mat']);   %cell mask (mat file)
fnames.shadecells=fullfile(dirs.label,[fname,'_shadecells.tif']);   %cell mask displayed on avg img
fnames.tcourse=fullfile(dirs.tcourse,[fname,'_tcourse.mat']);   %time courses of cells
fnames.stats=fullfile(dirs.stats,[fname,'_stats.mat']); %basic statistics
fnames.Oristats=fullfile(dirs.oristats,[fname,'_Oristats.mat']); %orientation statistics

end

end
