function [periods,stats]=ijdad2tiffseq(sourcepath,targetpath, varargin)
%IJDAD2TIFFSEQ Converts fast scanner raw dad files into sequences of TIFF files
% [PERIODS]=IJDAD2TIFFSEQ(SOURCEPATH,TARGETPATH, VARARGIN)
% VARARGIN are one or more of the following options as property pairs
% 	FilesN - Index vector specifying which dad files to process
%	Channels - Channels to reconstruct  [ 0|{1} , 0|{1}]
%   AverageN - Downsample by this many frames {1}
%   Limits - Pixel range {[2048,4096]}
%   FixDelay - Fix delay (typically ranges from 15 to 30 pixels)
%   Depth - Bit depth of reconstructed images [8,{16}]
%   RecalculateDelay - Recalculate delay at every frame [ true | {false} ]
%   ContinuousData - How many files till restart of acquisition, for
%       z-stacks.  Empty means data is continuous  [integer or { [] } ]
%   ForcedNumberOfDataFrames - forces to have N frames in each dad file
%       (truncates/pads with zeros if necessary)  Used only if 'ContinuousData'=[]
%       Empty means do not force { [] }
%   RawReconstruction - reconstruct image without cos rescaling [ true | {false} ]
%   WaitForDads - Start reconstruction even if not all files available, reconstruct them 
%       as they arrive in directory.  Allows online use.  Must set FilesN. [true|{false}]
%
% by Sergey Yurgenson, Vincent Bonin
%
% e.g. [periods]=ijdad2tiffseq('c:\sourcedir','d:\targetdir', 'FilesN',1)
%

% Change log
% 04/28/07 vb let specify channels
% 05/03/07 vb optimized. faster zero crossing calculation.
%             replaced wrapper imwrite with private function wtifc
% 05/12/07 vb averages N frames before writing to tif file
% 05/24/07 vb reorganized args, now reads file according to mod date
% 05/29/07 vb now creates both 16-bit 8-bit tiff files
%          vb reconstructs both channels
%          vb chanstr takes 'green', 'red' or 'both'
%          vb reconstruct straight to 8 bit
% 10/23/07 vb ported to imagej library
% 03/19/08 sy delay is two-element matrix, first element: delay (or NaN),
%             second element - delay recalculation flag: 0 - calculate
%             ones, 1 - calculate for each data file separetely
% 03/21/08 vb corrected filelist bug, default
%          vb redesigned input arguments to allow for unlimited options
% 04/02/08 sy changed zero crossing detection algorithm to work with noisy
%             data, overload finddaddelay to accept data matix itself
% 04/07/08 sy additional modification of zero crossing detection
%             (lin. regression), channels for delay calculation now determined
%             by 'Channels'
% 04/10/08 SY new parameters: 'Width' - image width, 'OverSampl' - maximum
%             number of data points to be averaged for each image pixel
% 05/30/08 SY new parameters: 'ForcedNumberOfDataFrames' - forces to have N frames in each dad file
%             Cuts dad if necessary, pads zeros at the file end if
%             necessary. Parameter used only if 'ContinuousData'=false
% 06/10/08 SY new parameter 'RawReconstruction'- true/false to reconstruct
%             image without cos rescaling
% 06/23/08 SY changed delay recalculation, 'RecalculateDelay' is integer
%             now; -1 - no recalculation, 0 - at the start of each file,
%             n - recalculate delay every n frames
% 06/24/08 SY performance optimization
% 07/22/08 SY performance optimization, use of streamlined function
%             finddaddelay2
% 10/29/08 SY finddaddelay2 now applied over window of size RecalculateDelay
% 10/29/08 VB disp options included in diary
% 01/05/08 MH reconstruct z stacks, use ContinuousData = nFilesPerZLevel
% 01/05/08 MH add WaitForDads option to speed up online use
% 03/23/09 SY add SamplingRate option [MHz]

%% default arguments
defaultopts = {'FilesN',[],'Channels',[1 1],'AverageN',1,'FixDelay',[],'Limits',[2048,4096],'Depth',16,'RecalculateDelay',-1,'StackSize',1000 ...
    'DiaryMode', 'on','ContinuousData',[],'Width',256,'OverSampl',1,'ForcedNumberOfDataFrames',[],'RawReconstruction', false,'Overwrite',false,'WaitForDads',false,'SamplingRate',3.3};

options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);
for iarg = 1:2:length(varargin)
    options.(varargin{iarg}) =varargin{iarg+1};
end

height = 239;
width = options.Width;
format compact

[dPath dSeriesName dExt] = fileparts(sourcepath);


currentsingleperiod=845.5*options.SamplingRate/3.3; %initial estimation


%% deal with files / directories
% each if option sets nfiles, filelist, indices
if ~options.WaitForDads
    filelistS = dir(fullfile(sourcepath,'*.dad'));
    nfiles = length(filelistS);
    fprintf(1,'Found %i dad files to process\n',nfiles);
    if nfiles==0
        diary off;
        error('No dad files found at %s', sourcepath);
    end

    if isempty(options.FilesN); options.FilesN = 1:nfiles; end;

    filelist = {filelistS.name};
    [vals,indices]=sort(datenum({filelistS.date})); % sort by modification date    
else
    assert(~isempty(options.FilesN), 'Must set FilesN if WaitForDads == true');
    assert(min(options.FilesN) == 1, 'FilesN starts numbering from one');
    options.FilesN = sort(options.FilesN);
    
    % make full file list from numbers and series name
    nfiles = length(options.FilesN);
    for iF = 1:nfiles
        filelist{iF} = sprintf('%s_%06d.dad', dSeriesName, options.FilesN(iF));
    end
    indices = 1:nfiles;
end
    
% save some information for debugging
stats.options = options;
stats.nFoundFiles = nfiles;
stats.FilesN = options.FilesN;
stats.nUsedFiles = length(options.FilesN);



if exist(targetpath)==7

    if options.Overwrite 
        str = 'y';
    else
        prompt = sprintf('Directory already exists. Delete?\n');
        str = input(prompt,'s');
    end

    if strcmp(str,'y') 
        %rmdir(targetpath,'s');
        fprintf(1,'Deleted existing directory\n');
        disp(targetpath);
    else
        return;
    end
end

[pathstr,name]=fileparts(targetpath);

%%
targetgreen = fullfile(targetpath,[name '_green']);
targetred = fullfile(targetpath,[name '_red']);

mkdir(targetpath);
mkdir(targetgreen);
mkdir(targetred);

stats.targetPath = targetpath;
stats.targetGreen = targetgreen;
stats.targetRed = targetred;

%% sends standard output to log file
if strcmp(options.DiaryMode,'on')
    diary(fullfile(targetpath,[name '.log']));
end

basefngreen = fullfile(targetgreen,name);
basefnred = fullfile(targetred,name);

options

%% initialization
periods = [];

if options.RawReconstruction %calculate image width if using raw reconstract (no sin function)
    % TODO - replace algorithm by one used in main code
    thisfn = filelist{indices(1)};
    inputfn = fullfile(sourcepath,thisfn);
    fid= fopen(inputfn,'r');
    [out, count]=fread(fid,inf,'uint16=>int16');
    fclose(fid);
    out = reshape(out,4,count/4);
    crossup=find(out(3,1:end-6)<=2048 & out(3,2:end-5)>2048 & out(3,3:end-4)>2048 & out(3,4:end-3)>2048 & out(3,5:end-2)>2048 & out(3,6:end-1)>2048 & out(3,7:end)>2048);
    dind=find(diff(crossup)<currentsingleperiod-200 | diff(crossup)>currentsingleperiod+200);
    if ~isempty(dind)
        fprintf(1,'Has discontinuity at data point # ~ %6d \n',crossup(dind)+round(currentsingleperiod/2));
    end
    for r=1:length(crossup)
        firstpoint=min(crossup(r)-1,20);
        bbb = regress_no_error_check(double(out(3,max(1,crossup(r)-firstpoint):crossup(r)+20)-2048)',[ones(firstpoint+21,1) (-firstpoint:20)']);
        crossup(r)=crossup(r)-bbb(1)/bbb(2);
    end
    if isempty(dind)
        period=(crossup(end)-crossup(1))/(length(crossup)-1);
    else
        period=(crossup(dind(1))-crossup(1))/(dind(1)-1);
    end
    width=floor(period);
    height = 120;
end

fgreen=ones(256,width*options.OverSampl,'single')*2048;
fgreensum=zeros(256,width*options.OverSampl,'single');
greenStack = ij.ImageStack(width,height);

fred=ones(256,width*options.OverSampl,'single')*2048;
fredsum=zeros(256,width*options.OverSampl,'single');
redStack = ij.ImageStack(width,height);

processedFileN=0;
processedFrameN = 0;
savedGreenFrameN = 0;
savedGreenFileN = 0;
savedRedFrameN = 0;
savedRedFileN = 0;

previous = [];


if options.Depth == 16
    a = single(2048);
    b = single(options.AverageN);
elseif options.Depth == 8
    a = single(4096-options.Limits(1));
    b = single(255/( options.Limits(2)-options.Limits(1) )/options.AverageN);
end

tic;

if ~isempty(options.FixDelay)
    delay = options.FixDelay;
end

delchan=1:2;
delchan(options.Channels==0)=[];

%% is z stack?
if ~isempty(options.ContinuousData)
    isZ = true;
    filesPerLevel = options.ContinuousData;
    assert(isZ && length(indices) == length(options.FilesN), ...
           'Z stack: must use all files'); 
else
    isZ = false;
    filesPerLevel = length(options.FilesN);
end



%% process dad files
stats.discontPoints = [];
stats.framesPerDadFile = [];
for fileIndex=options.FilesN
    thisnum = indices(fileIndex);
    thisfn = filelist{thisnum};
    inputfn = fullfile(sourcepath,thisfn);
    
    fprintf(1,'Processing file %i/%i (%s)...\n',fileIndex,options.FilesN(end),thisfn);

    % wait if requested, unlimited waiting
    iS=1;
    if options.WaitForDads == true
        while true
            if ~exist(inputfn, 'file')
                fprintf(1, '%3ds: waiting for file: %s\n', iS, inputfn);
                pause(1);
                iS=iS+1;                
            else
                % found file
                break;
            end
        end
    end

    % open file
    fid= fopen(inputfn,'r');
    if fid<0
        display(['Unable to open ' inputfn]);
        break;
    end
    % for some reason reading native uint16 from disk very slow
    [out, count]=fread(fid,inf,'uint16=>int16');
    fclose(fid);

    % zstack?
    isZFirstFile = mod(thisnum, filesPerLevel) == 1;
    isZLastFile = mod(thisnum, filesPerLevel) == 0;
    if (thisnum == 1 ...
        || (isZ && isZFirstFile) )
        % first file of a zstep.
        out = reshape(out,4,count/4);
        %fprintf(1, '** thisnum %d: started at beginning\n', thisnum);
    else        
        out = [previous reshape(out,4,count/4)];
        %fprintf(1, '** thisnum %d: continuing\n', thisnum);
    end


    %check for trigger position
    if isZFirstFile || fileIndex==options.FilesN(1)
        yPosMax=max(out(4,:));
        yPosMin=min(out(4,:));
        yPosTshld=yPosMin+(yPosMax-yPosMin)/18;
        ycross=find(out(4,:)<yPosTshld);
        if ycross(1)>1
            out=out(:,ycross(1)+1:end);
        end
    end

    %    chImd=double(out(delchan(1),:));
    chIms=single(out(delchan(1),:));
    if length(delchan)>1
        %    chImd=chImd+double(out(delchan(2),:));
        chIms=chIms+single(out(delchan(2),:));
    end



    %find all crossing up
%    firstcrossup=find(out(3,1:1000)<=2048 & out(3,2:1001)>2048 & out(3,3:1002)>2048 & out(3,4:1003)>2048 & out(3,5:1004)>2048 & out(3,6:1005)>2048 & out(3,7:1006)>2048);
% change to accomodate higth sampling rate    
firstcrossup=find(out(3,1:round(currentsingleperiod*1.2))<=2048 & out(3,2:round(currentsingleperiod*1.2)+1)>2048 & out(3,3:round(currentsingleperiod*1.2)+2)>2048 & out(3,4:round(currentsingleperiod*1.2)+3)>2048 & out(3,5:round(currentsingleperiod*1.2)+4)>2048 & out(3,6:round(currentsingleperiod*1.2)+5)>2048 & out(3,7:round(currentsingleperiod*1.2)+6)>2048);

    newcrossup=zeros(1,100000);
    newcrossup(1)=firstcrossup(1);
    cr=1;
    mmm=size(out,2);
    lastfirstpoint=100;
    while 1
        rnewcrossup=round(newcrossup(cr));
        firstpoint=min(rnewcrossup-1,20);
        %bbb = regress_no_error_check(double(out(3,max(1,rnewcrossup-firstpoint):rnewcrossup+20)-2048)',[ones(firstpoint+21,1) (-firstpoint:20)']);

        % next is a copy of [b,bint,r,rint,stats] = regress(y,X,alpha) with
        % some optimization
        y=double(out(3,max(1,rnewcrossup-firstpoint):rnewcrossup+20)-2048)';
        if firstpoint~=lastfirstpoint
            X=[ones(firstpoint+21,1) (-firstpoint:20)'];
            [Q,R,perm] = qr(X,0);
            bbb = zeros(2,1);
        end
        lastfirstpoint=firstpoint;
        bbb(perm) = R \ (Q'*y);
        newcrossup(cr)=newcrossup(cr)-bbb(1)/bbb(2);
        %if newcrossup(cr)+866< mmm
        if newcrossup(cr)+currentsingleperiod+20< mmm    
            newcrossup(cr+1)=newcrossup(cr)+currentsingleperiod;
            if cr>2
                currentsingleperiod=newcrossup(cr)-newcrossup(cr-1);
            end
            cr=cr+1;
        else
            break
        end
    end
    crossup=newcrossup(1:cr);

    dind=find(diff(crossup)<currentsingleperiod-200 | diff(crossup)>currentsingleperiod+200);
    if ~isempty(dind)
        fprintf(1,'Has discontinuity at data point # ~ %6d \n',crossup(dind)+round(currentsingleperiod/2));
        stats.discontPoints(end+1) = crossup(dind)+round(currentsingleperiod/2);
    end

    % recalculate delay for each file if it is first file or
    % RecalculateDelay=0
    if (isempty(options.FixDelay) && processedFrameN == 0) || (isempty(options.FixDelay) && options.RecalculateDelay==0)
        if processedFileN==0
%            fprintf(1,'Calculating delay... ');
            [period,delay]=finddaddelay2(1,150,0,40,1,width,crossup,chIms);
            %        [period,delay]=finddaddelay2(1,150,0,40,1,width,crossup,chImd);
        else
%            fprintf(1,'Calculating delay... ');
            [period,delay]=finddaddelay2(1,150,delay-5,delay+5,1,width,crossup,chIms);
            %        [period,delay]=finddaddelay2(1,150,delay-5,delay+5,1,width,crossup,chImd);
        end
    end
    period=(crossup(end)-crossup(1))/(length(crossup)-1);
    fprintf(1,'Initial delay is %2.1f pixels; period %6.3f\n', delay, period);

    periods(fileIndex)=period;
    if crossup(1)>round(period/4)
        linestartInd=1;
    else
        linestartInd=2; % vb changed from 2 to 1 on 080116
    end

    frameIndices=1:10000;
    if ~isempty(options.ForcedNumberOfDataFrames)
        if options.ForcedNumberOfDataFrames*128 +linestartInd < length(crossup)
            %cut data file if necessary
            crossup(options.ForcedNumberOfDataFrames*128 +2:end)=[];
        elseif options.ForcedNumberOfDataFrames*128 +linestartInd > length(crossup)
            %append data file if necessary
            extracrossups=1:(options.ForcedNumberOfDataFrames*128 +linestartInd-length(crossup));
            crossup(end+1:end+length(extracrossups))=crossup(end)+extracrossups*period;
            out(:,end:round(crossup(end)) + 1)=0;
        end
    end

    w = warning('off');

    framele=round(period*120); % frames are 120 cycles long
    if options.RecalculateDelay<1
        if options.RawReconstruction
            yda=delay:framele-1+delay;
            xii=uint32(floor(mod(yda,period)));
            xii(xii>width)=width;
            yii=uint32(yda*119.99999/(framele-1)+0.5);
            ind = yii + (xii-1)*256;
        else
            step=2*pi/period;
            yda=delay:framele-1+delay;
            %xii=uint32((1-cos(yda*step))*(width*options.OverSampl-1)/2+1);
            %yii=uint32(yda*239.99999/(framele-1)+0.5);
            %ind=sub2ind_no_error_check([256,width*options.OverSampl],yii,xii); % transposing
            ind=uint32(yda*239.99999/(framele-1)+0.5) + (uint32((1-cos(yda*step))*(width*options.OverSampl-1)/2))*256;
        end
    end
    warning(w);
    thiscrossup = [];
    thisdadframes = 0;
    for frameIndex=frameIndices

        lastcrossup = thiscrossup;
        thiscrossup = linestartInd+(frameIndex-1)*128;
        if thiscrossup + 128 > length(crossup)
       %     if fileIndex==options.FilesN(end)
             if (fileIndex==options.FilesN(end) || (isZ && isZLastFile) )
                if thiscrossup == length(crossup)
                    fprintf(1,'Processed %i frames in %2.1f seconds (%2.0f frames/seconds)\n',...
                        processedFrameN, toc, processedFrameN/toc);
                    break
                else
                    %append last data file if necessary
                    extracrossups=1:(thiscrossup + 128-length(crossup));
                    crossup(end+1:end+length(extracrossups))=crossup(end)+extracrossups*period;
                    out(:,end:round(crossup(end)) + 1)=0;
                    %    chImd=double(out(delchan(1),:));
                    chIms=single(out(delchan(1),:));
                    if length(delchan)>1
                        %    chImd=chImd+double(out(delchan(2),:));
                        chIms=chIms+single(out(delchan(2),:));
                    end
                end
            else
                previous = out(:,round(crossup( lastcrossup + 127 )) + 1 : end);
                fprintf(1,'Processed %i frames in %2.1f seconds (%2.0f frames/seconds)\n',...
                    processedFrameN, toc, processedFrameN/toc);
                break
            end
        end


        if isempty(options.FixDelay) && options.RecalculateDelay>0 && mod(frameIndex-1,options.RecalculateDelay)==0
            %recalculate delay each time
            frameIndexRange=max(0,options.RecalculateDelay-1);
            if frameIndex==1
                [period22,delay22]=finddaddelay2(max(1,frameIndex-frameIndexRange),frameIndex+frameIndexRange,delay-10,delay+10,1,width,crossup,chIms);
                %                [period22,delay22]=finddaddelay2(max(1,frameIndex-frameIndexRange),frameIndex+frameIndexRange,delay-10,delay+10,1,width,crossup,chImd);
            else
                [period22,delay22]=finddaddelay2(max(1,frameIndex-frameIndexRange),frameIndex+frameIndexRange,delay-5,delay+5,1,width,crossup,chIms);
                %                [period22,delay22]=finddaddelay2(max(1,frameIndex-frameIndexRange),frameIndex+frameIndexRange,delay-5,delay+5,1,width,crossup,chImd);
            end
            delay=delay22;
            step=2*pi/period22;
            framele=round(crossup((frameIndex-1)*128+1+120)-crossup((frameIndex-1)*128+1));
            yda=delay:framele-1+delay;
            if options.RawReconstruction
                xii=uint32(floor(mod(yda,period)));
                xii(xii>width)=width;
                yii=uint32(yda*119.99999/(framele-1)+0.5);
                ind = yii + (xii-1)*256;
            else
                ind=uint32(yda*239.99999/(framele-1)+0.5) + (uint32((1-cos(yda*step))*(width*options.OverSampl-1)/2))*256;
            end

        end

        %****** end of delay recalculation
        frameIndexMod = mod(processedFrameN,options.AverageN);
        % next section was moved up
        %{
        lastcrossup = thiscrossup;
        thiscrossup = linestartInd+(frameIndex-1)*128;
        if thiscrossup + 128 > length(crossup)
            previous = out(:,round(crossup( lastcrossup + 127 )) + 1 : end);
            fprintf(1,'Processed %i frames in %2.1f seconds (%2.0f frames/seconds)\n',...
                processedFrameN, toc, processedFrameN/toc);
            break
        end
        %}
        framestart=round(crossup(thiscrossup)-period/4);

        % update green channel
        if options.Channels(1)
            if options.OverSampl>1
                fgreen(:)=NaN;
                fgreen(ind) = out(1,framestart:framestart+framele-1);

                if frameIndexMod == 0
                    fgreensum = reshape(nanmean(reshape(fgreen,256,options.OverSampl,width),2),256,width);
                else
                    fgreensum = fgreensum + reshape(nanmean(reshape(fgreen,256,options.OverSampl,width),2),256,width);
                end
                fgreensum(isnan(fgreensum))=2048;
            else
                fgreen(:)=2048;
                fgreen(ind) = out(1,framestart:framestart+framele-1);
                if frameIndexMod == 0
                    fgreensum = fgreen;
                else
                    fgreensum = fgreensum + fgreen;
                end
            end

            if frameIndexMod==options.AverageN-1
                fgreenav = max(0,(a-fgreensum) * b);
                switch options.Depth
                    case 8
                        greenProcessor = ij.process.ByteProcessor(width,height);
                    case 16
                        greenProcessor = ij.process.ShortProcessor(width,height);
                end
                greenProcessor.setFloatArray(fgreenav(1:height,1:width)')
                greenStack.addSlice('frame', greenProcessor);
                savedGreenFrameN = savedGreenFrameN + 1 ;
                if mod(savedGreenFrameN,options.StackSize) == 0
                    greenPlus = ij.ImagePlus('GreenStack',greenStack);
                    greenFileSaver = ij.io.FileSaver(greenPlus);
                    fn = sprintf('%s_%06d.tif',basefngreen,savedGreenFileN);
                    greenFileSaver.saveAsTiffStack(fn);
                    greenStack = ij.ImageStack(width,height);
                    savedGreenFileN = savedGreenFileN + 1 ;
                    disp(fn);
                end
            end
        end
        % update red channel
        % TODO - make red cannel the same as green
        if options.Channels(2)
            fred(:)=2048;
            fred(ind) = out(2,framestart:framestart+framele-1);
            if frameIndexMod == 0
                fredsum = single(fred);
            else
                fredsum = fredsum + single(fred);
            end
            if frameIndexMod==options.AverageN-1
                fredav = max(0,(a-fredsum) * b);

                switch options.Depth
                    case 8
                        redProcessor = ij.process.ByteProcessor(width,height);
                    case 16
                        redProcessor = ij.process.ShortProcessor(width,height);
                end
                redProcessor.setFloatArray(fredav(1:height,1:width)');
                redStack.addSlice('frame', redProcessor);
                savedRedFrameN = savedRedFrameN + 1 ;

                if mod(savedRedFrameN,options.StackSize) == 0
                    redPlus = ij.ImagePlus('RedStack',redStack);
                    redFileSaver = ij.io.FileSaver(redPlus);
                    fn = sprintf('%s_%06d.tif',basefnred,savedRedFileN);
                    redFileSaver.saveAsTiffStack(fn);
                    redStack = ij.ImageStack(width,height);
                    savedRedFileN = savedRedFileN + 1 ;
                end
            end
        end
        processedFrameN = processedFrameN + 1;
        thisdadframes = thisdadframes+1;
    end % frameIndex

    out=[];
    crossup=[];
    sg=[];
    dsg=[];

    stats.framesPerDadFile(end+1) = thisdadframes;
    processedFileN=processedFileN+1;
end % fileIndex

%% save remaining frames

if mod(savedGreenFrameN,options.StackSize) ~= 0
    greenPlus = ij.ImagePlus('GreenStack',greenStack);
    greenFileSaver = ij.io.FileSaver(greenPlus);
    fn = sprintf('%s_%06d.tif',basefngreen,savedGreenFileN);
    greenFileSaver.saveAsTiffStack(fn);
    greenStack = ij.ImageStack(width,height);
    savedGreenFileN = savedGreenFileN + 1 ;
    disp(fn);
end

if mod(savedRedFrameN,options.StackSize) ~= 0
    redPlus = ij.ImagePlus('RedStack',redStack);
    redFileSaver = ij.io.FileSaver(redPlus);
    fn = sprintf('%s_%06d.tif',basefnred,savedRedFileN);
    redFileSaver.saveAsTiffStack(fn);
    redStack = ij.ImageStack(width,height);
    savedRedFileN = savedRedFileN + 1 ;
end

fprintf(1,'Saved %i green  frames to dir %s\n',savedGreenFrameN,targetgreen);
fprintf(1,'Saved %i red frames to dir %s\n',savedRedFrameN,targetred);

diary off;

% save some outputs

return;
