function frReconstruct(newdir,varargin)
%FRRECONSTRUCT Reconstruct fast scanner images of newly added directories 
% FRRECONSTRUCT checks entire raw directory
% FRRECONSTRUCT(NEWDIR) specifies raw directory
% FRRECONSTRUCT(... , VARARGIN) where VARARGIN is one or more 
% of the property pairs 
%   Overwrite - Overwrites existing image directory [ true | {false} ]
% 	FilesN - Index vector specifying which dad files to process
%	Channels - Channels to reconstruct  [ 0|{1} , 0|{1}]
%   AverageN - Downsample by this many frames {1}
%   Limits - Pixel range {[2048,4096]}
%   FixDelay - Fix delay (typically ranges from 15 to 30 pixels)
%   Depth - Bit depth of reconstructed images [8,{16}]
%   RecalculateDelay - Recalculate delay at every frame [ true | {false} ]
%   DoRecurse - Recurse through child dirs or just do this dir [ {true} | false ]
%   WaitForDads - start reconstruction even if dad dir is empty, wait for new dads
%       to arrive (number to wait for calculated from scan log [ true | {false} ]
%   Engine - Specify algorithm [ {'ijdad2tiffseq'} | 'ijdad2tiffseq2' ]
%
%   Automatically detects z-stack and sets options for ijdad2tiffseq
%
% by Vincent Bonin
% 080105 MH: add DoRecurse, WaitForDads option; allow use with zstacks and
%    for online reconstruction

if ~verLessThan('matlab','7.5')
    prevn = maxNumCompThreads(1); % marginal gain beyond 1 thread
end

defaultopts = {'Overwrite', false, 'DoRecurse', true, 'WaitForDads', false, 'Engine', 'ijdad2tiffseq' };
options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);
for iarg = 1:2:length(varargin)
    options.(varargin{iarg}) =varargin{iarg+1};
end

dirs = frGetDirs;

%% process and check directories
if exist(newdir)==7
    rawdir = newdir;
else
    rawdir = fullfile(dirs.data,newdir);
end

imdir = fullfile(dirs.images,newdir);

if exist(rawdir) ~= 7
    str = sprintf('raw directory ''%s'' does not exist',rawdir);
    error(str);
end

list = dir(rawdir);list(1:2)=[]; % deletes curent and parent dir
sel = find([list.isdir]);

if options.DoRecurse && length(sel) > 0 
    for index = 1:length(sel)
        frReconstruct(fullfile(newdir,list(index).name),varargin{:});
    end
    return;
end    

%% process only if reconstruction does not exist, or Overwrite true
if options.Overwrite || ~exist(imdir) 
    fprintf(1,'Source: %s\n',rawdir);
    fprintf(1,'Target: %s\n',imdir);
    
    %% load scan log
    sl = frParseScanLog(rawdir);

    %% get parameters for ijdad
    waitForDads = false;
    filesN = [];
    if ~isempty(sl)
        % check if zstack
        if sl.cycle_Do && sl.cycle_NSteps > 1
            filesPerLevel = ceil(sl.cycle_NFramesPrStep / sl.DAQ_dad_size_in_frames);

% VB cycle_Do is what determines whether a stack, otherwise crashes when
% reconstructing normal streams after having acquired stacks
%            assert(sl.cycle_Do == 1, 'bug: log file inconsistent?');
        else
            filesPerLevel = [];
        end
        
        if options.WaitForDads
            waitForDads = true;
            % detect filesN
            nDads = ceil(sl.cycle_NFramesPrStep / sl.DAQ_dad_size_in_frames);
            if sl.cycle_Do
                nDadPerLevel = ceil(sl.cycle_NFramesPrStep / sl.DAQ_dad_size_in_frames);
                nDads = nDadPerLevel .* sl.cycle_NSteps;
            end
            filesN = 1:nDads; %0:(nDads-1);
        end            
    else
        fprintf(1, 'No scan log found: old data files?\n');
        assert(options.WaitForDads == false, ...
               'Need scan log for WaitForDads: call ijdad directly');
    end    
    
    % user can overwrite FilesN with input argument
    if isfield(options,'FilesN') & length(options.FilesN)
        filesN = options.FilesN;
    end
        
    f = str2func(options.Engine);
    fprintf(1,'Calling %s\n',options.Engine);
    feval(f,rawdir, imdir, varargin{:}, ...
                  'ContinuousData', filesPerLevel, ...
                  'WaitForDads', waitForDads, ...
                  'FilesN', filesN);   
    
else
   fprintf(1,'Ignoring directory %s\n',rawdir);
end

if ~verLessThan('matlab','7.5')
    maxNumCompThreads(prevn); % restore default
end

return;
