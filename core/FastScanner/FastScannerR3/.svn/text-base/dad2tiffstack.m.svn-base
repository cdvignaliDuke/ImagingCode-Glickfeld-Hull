function [periods]=dad2tiffstack(sourcepath,outputfn,channel,averageN,filesN,delay)
%
% [periods]=dad2tiff(sourcepath,outputfn,channel,averageN,filesN,delay)
%
% by Sergey Yurgenson, Vincent Bonin
% 
%
% change log
% 04/28/07 vb let specify channel
% 05/03/07 vb optimized. faster zero crossing calculation. 
%             replaced wrapper imwrite with private function wtifc
% 05/12/07 vb averages N frames before writing to tif file
% 05/24/07 vb reorganized args, now reads file according to mod date
% 05/29/07 vb now creates both 16-bit 8-bit tiff files
%

diary([outputfn '.log']);

format compact

filelist = dir(sourcepath);filelist = filelist(3:end);
nfiles = length(filelist);

fprintf(1,'Found %i dad files to process\n',nfiles);
[vals,indices]=sort([filelist.datenum]); % sort by modification date

outputfn16 = [outputfn '_16.tif'];

if nargin < 3;channel =1; end;
if nargin < 4 ;averageN = 1; end;
if nargin < 5; filesN = [];end;
if isempty(filesN);filesN = 1:nfiles; end;
if nargin < 6;
    fprintf(1,'Calculating delay... ');
    inputfn = fullfile(sourcepath,filelist(indices(filesN(1))).name);
    [period,delay]=finddaddelay(inputfn,1,150,2,10,40,0.5);
    fprintf(1,'Delay is %2.1f pixels\n',delay);
end

fi=zeros(256,256,'uint16');
fii=zeros(256,256,'uint32');

fmin = inf;
fmax = -inf;

overwrite = 1;
processedFrameN = 0; 
savedFrameN = 0;

for fileIndex=filesN
    thisfn = filelist(indices(fileIndex)).name;
    inputfn = fullfile(sourcepath,thisfn);
	fprintf(1,'Processing file %i/%i (%s)...\n',fileIndex,filesN(end),thisfn);
    
    tic;
    fid= fopen(inputfn,'r');
    
    if fid<0
        display(['Unable to open ' inputfn]);
        break;
    end

	[out, count]=fread(fid,inf,'uint16=>int16');
    fclose(fid);
    
    out=reshape(out,4,count/4);
    out = out - 2048;
    
    %% calculate period
    xpos = out(3,:);
    sg = xpos>=0;
    dsg=diff(sg);
    crossup=find(dsg>0); % sample before upward zero crossing
    period=(crossup(end)-crossup(1))/(length(crossup)-1);
    periods(fileIndex)=period;
    
    if crossup(1)>round(period/4)
         linestartInd=1;
    else
        linestartInd=2;
    end
    
    framele=round(period*120); % frames are 120 cycles long
    step=2*pi/period;
    yda=delay:framele-1+delay;
    xii=uint16((1-cos(yda*step))*255/2+1);
    yii=uint16(yda*239.99999/(framele-1)+0.5);
	% ind=sub2ind_no_error_check([256,256],xi,yii);
	ind=sub2ind_no_error_check([256,256],yii,xii); % transposing
    
    xii=[];
    yii=[];
    yda=[];
    
    frameIndices=1:10000;

   for frameIndex=frameIndices       
        frameIndexMod = mod(processedFrameN,averageN);
        
        if linestartInd+(frameIndex-1)*128>length(crossup)
            frames_in_file=frameIndex-1;
            fprintf(1,'Processed %i frames for a total of %i frames in %2.0f seconds \n',...
                    frames_in_file, processedFrameN, toc);
            break
        end       
        
        linestart=round(crossup(linestartInd+(frameIndex-1)*128)-period/4);

        fi(ind) = 2048 - out(channel,linestart:linestart+framele-1);
        
        if frameIndexMod == 0
            fii = uint32(fi);
        else
            fii = fii + uint32(fi);
        end

        % overwrite = fileIndex==firstFileN && frameIndex==frameIndices(1);
        
        %  wtifc(data, map, outputfn, compression, description,
        %        resolution, writemode, colorspace);    
        
        if frameIndexMod==averageN-1                        
            fii = fii / averageN ;

            fmin = min(fmin,min(fii(find(fi>0))));
            fmax = max(fmax,max(fii(:)));
            
            wtifc(uint16(fii), [], outputfn16, 'none', '', 72, overwrite, 'rgb');            
            overwrite = 0;
            savedFrameN = savedFrameN + 1 ;
        end
        %imwrite(fi,outputfn,'tif','Compression','none','WriteMode','overwrite'); 
        processedFrameN = processedFrameN + 1;
    end
    
    out=[];
    ind=[];
    crossup=[];
    
end

if averageN == inf
    fii = fii / processedFrameN;

    fmin = min(fmin,min(fii(find(fi>0))));
    fmax = max(fmax,max(fii(:)));
    
    wtifc(uint16(fii), [], outputfn16, 'none', '', 72, overwrite, 'rgb');
    overwrite = 0;
	savedFrameN = savedFrameN + 1 ;
end

fprintf(1,'Saved %i 16-bit frames to file %s\n',savedFrameN,outputfn);

fprintf(1,'Now converting to 8 bit tif\n');
outputfn8 = [outputfn '_8.tif'];
tif16totif8(outputfn16,outputfn8,[fmin,fmax]);

diary off;

return;