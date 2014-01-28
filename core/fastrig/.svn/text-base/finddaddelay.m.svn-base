function [period,mean_del2]=finddaddelay(filename,frameStart,frameEnd,ichan,mindel,maxdel,Dstep,width,tcrossup)
%
% [del,mean_del2]=finddaddelay(filename,frameStart,frameEnd,ichan,mindel,maxdel,Dstep)
%
%function reads data file (.dad) and determines best shift(delay) for image
%reconstruction
%
%example:
%[d1,d2]=finddaddelay('C:\Data\FastScanner\pipette48.dad',1,32,1,10,40,0.5);
% in example file pipette48.dad, frames from 1 to 32, image in channel 1
% were used. Shift was probed from 10 to 40 with step 0.5
% SY 04/18/2007

% changelog
% 05/03/07 vb optimized
% 05/16/07 vb corrected rounding error, renamed to finddaddelay
% 04/01/08 sy filename can be name of data file or data itself
% 04/07/08 sy ichan can be [1 2] now, use of lin. regression for zero
%          crossing calculation
% 04/10/08 sy new parameter width - image width

% 06/24/08 SY performance optimization

format compact

fi=zeros(width,256);

if nargin<2
    frameStart=1;
end

if nargin<3
    frameEnd=1;
end

if nargin<4
    ichan=1;
end
if nargin<5
    mindel=1;
end
if nargin<6
    maxdel=1;
end

if nargin<7
    Dstep=1;
end

if nargin<9
    tcrossup=[];
end



if ischar(filename)
    %open and load file
    fid= fopen(filename,'r');
    if fid<0
        error(['Unable to open ' filename]);
        return
    end
    [out, count]=fread(fid,inf,'uint16');
    fclose(fid);

    %extract image and fast mirror signal
    out=reshape(out,4,count/4);
    xpos=out(3,:);
    imgdata=out(ichan,:);
    xpos=xpos-2048;
    out=[];
else
    if isempty(tcrossup)
        xpos=double(filename(3,:)-2048);
    end
end

%imgdata=sum(imgdata,1);
if isempty(tcrossup)
    % following method works better for noisy xdata
    % TODO - make the same as in ijdad2tiffseq
    tcrossup=find(xpos(1:end-6)<=0 & xpos(2:end-5)>0 & xpos(3:end-4)>0 & xpos(4:end-3)>0 & xpos(7:end)>0);

    for r=1:length(tcrossup)
        firstpoint=min(tcrossup(r)-1,10);
        bbb = regress_no_error_check(double(xpos(tcrossup(r)-firstpoint:tcrossup(r)+10))',[ones(firstpoint+11,1) (-firstpoint:10)']);
        tcrossup(r)=tcrossup(r)-bbb(1)/bbb(2);
    end
end
period=(tcrossup((frameStart-1)*128+1+120)-tcrossup((frameStart-1)*128+1))/120;

if tcrossup(1)>round(period/4)
    linestartInd=1;
else
    linestartInd=2;
end

framele=round(period*120);
nframes = floor(length(tcrossup)-linestartInd + 1)/128; % number of complete frames
N=frameStart:min(nframes,frameEnd);
linestarts=round(tcrossup(linestartInd+(N-1)*128)-period/4);

%avg frames data together
% gi - avereged frame data

gi=zeros(1,framele);
if length(ichan)==2
%    fnm1=filename(ichan(1),:);
%    fnm2=filename(ichan(2),:);
%    fnm=double(fnm1+fnm2);
    
    %good for shorter data
    fnm1a=filename(ichan(1),linestarts(frameStart+1-frameStart):linestarts(min(nframes,frameEnd)+1-frameStart)+framele-1);
    fnm2a=filename(ichan(2),linestarts(frameStart+1-frameStart):linestarts(min(nframes,frameEnd)+1-frameStart)+framele-1);
    fnm=double(fnm1a+fnm2a);
    
    %good for longer data
%    fnmb=double(sum(filename(ichan,linestarts(frameStart+1-frameStart):linestarts(min(nframes,frameEnd)+1-frameStart)+framele-1),1));
    
else
    fnm=double(filename(ichan(1),:));
end

%fnm=double(sum(filename(ichan,:),1));

for n=N+1-frameStart
    if ~ischar(filename)
%        gi=gi+fnm(linestarts(n):linestarts(n)+framele-1);
        gi=gi+fnm(linestarts(n)-linestarts(1)+1:linestarts(n)+framele-1-linestarts(1)+1);
   else
        gi=gi+imgdata(linestarts(n):linestarts(n)+framele-1)/(frameEnd-frameStart+1);
    end
end

% for each delay reconstract image and calculate "combined vertical diff"
w = warning('off');
step=2*pi/period;
delN=0;
for delay=mindel:Dstep:maxdel
    delN=delN+1;
    yda=delay:framele-1+delay;
  %      xii=uint32((1-cos(yda*step))*(width-1)/2+1);
  %      yii=uint32(yda*239.99999/(framele-1)+0.5);
  %      ind2=sub2ind_no_error_check([width,256],xii,yii);
    ind=uint32(yda*239.99999/(framele-1)+0.5-1)*width + uint32((1-cos(yda*step))*(width-1)/2+1);
    
    xii=[];
    yii=[];
    yda=[];
    fi(ind)=gi;
    fii=reshape(fi,width,256);
%    d1(delN)=mean(mean(abs(diff(fii,1,2))));
    d1(delN)=sum(sum(abs(diff(fii,1,2))));
    ind=[];
end
warning(w);

%find minimum and assume it is best shift
%[A,I]=min(d1);
%mean_del2=mindel+(I-1)*Dstep;

% metod # 2 - find min assuming that function is symmetrical
%figure
%plot(d1)
[minD,I]=min(d1);
branch1=d1(1:I);
branch2=d1(I:end);

le=length(d1);
shiftb1=mindel+(0:I-1)*Dstep;
shiftb2=mindel+(I-1:le-1)*Dstep;

%minD=A;
maxD=min(max(branch1),max(branch2));
%maxD=minD+(maxD-minD)*0.9;
maxD=minD+(maxD-minD)*0.5;
q=linspace(minD+0.05*(maxD-minD),maxD,1000);

branch1(diff(branch1)==0)=branch1(diff(branch1)==0)+0.000001;
branch2(diff(branch2)==0)=branch2(diff(branch2)==0)-0.000001;


sh1=interp1(branch1,shiftb1,q);
sh2=interp1(branch2,shiftb2,q);

del2=(sh1+sh2)/2;
mean_del2=mean(del2);
