function update_images(newdir,varargin)
%UPDATE_IMAGES Reconstruct fast scanner images of newly added directories 
% UPDATE_IMAGES checks entire raw directory
% UPDATE_IMAGES(NEWDIR) specifies raw directory
% UPDATE_IMAGES(... , VARARGIN) where VARARGIN is one of more 
% of the property pairs 
%   Overwrite - Overwrites existing image directory [ true | false ]
% 	FilesN - Index vector specifying which dad files to process
%	Channels - Channels to reconstruct  [ 0|{1} , 0|{1}]
%   AverageN - Downsample by this many frames {1}
%   Limits - Pixel range {[2048,4096]}
%   FixDelay - Fix delay (typically ranges from 15 to 30 pixels)
%   Depth - Bit depth of reconstructed images [8,{16}]
%   RecalculateDelay - Recalculate delay at every frame [ true | {false} ]

% image directory be overwritten (default FALSE).
%
% by Vincent Bonin

defaultopts = {'Overwrite',false};
options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);

for iarg = 1:2:length(varargin)
    options.(varargin{iarg}) =varargin{iarg+1}
end

dirs = setupdirs;

%% recursively process children directories
rawdir = fullfile(dirs.data,newdir);
imdir = fullfile(dirs.images,newdir);

if exist(rawdir) ~= 7
    str = sprintf('raw directory ''%s'' does not exist',rawdir);
    error(str);
end

list = dir(rawdir);list(1:2)=[]; % deletes curent and parent dir
sel = find([list.isdir]);

if length(sel) 
    for index = 1 : length(sel)
        update_images(fullfile(newdir,list(index).name),varargin{:});
    end
    return;
end    

%% process only if reconstruction does not exist
if options.Overwrite || ~exist(imdir) 
    fprintf(1,'Source: %s\n',rawdir);
    fprintf(1,'Target: %s\n',imdir);

    ijdad2tiffseq(rawdir,imdir,varargin{:});
else
   fprintf(1,'Ignoring directory %s\n',rawdir);
end


return;