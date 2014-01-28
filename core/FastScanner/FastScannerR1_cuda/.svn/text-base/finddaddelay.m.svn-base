function [period,mean_del2]=finddaddelay(filename,frameStart,frameEnd,ichan,mindel,maxdel,Dstep)
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

format compact

fi=zeros(256,256);
fii=zeros(256,256);

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

if nargin<6
    Dstep=1;
end

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

sg = xpos>0;
dsg=diff(sg);
sg=[];
crossup=find(dsg>0); % sample before upward zero crossing
dsg = [];
period=(crossup(end)-crossup(1))/(length(crossup)-1);
tcrossup = crossup;
 
% % find zero crossing by mirror
% xpos(xpos==0)=0.00000001;
% sg=sign(xpos);
% dsg=diff(sg);
% sg=[];
% crossup=find(dsg==2);
% dsg=[];
% tcrossup=double(crossup)+abs(double(xpos(crossup)))./double((xpos(crossup+1)-xpos(crossup)));
% crossup=[];
% xpos=[];
% 
% %determine mirror period and initial start of lines
% period=(tcrossup(end)-tcrossup(1))/(length(tcrossup)-1);
% period

if tcrossup(1)>round(period/4)
    linestartInd=1;
else
    linestartInd=2;
end

framele=round(period*120);
nframes = floor(length(tcrossup)-linestartInd + 1)/128; % number of complete nframes
N=frameStart:min(nframes,frameEnd);
linestarts=round(tcrossup(linestartInd+(N-1)*128)-period/4);

%avg frames data together
% gi - avereged frame data
gi=zeros(1,framele);
for n=N
    if linestartInd+(n-1)*128>length(tcrossup)
        frames_in_file=n-1
        break
    end
    gi=gi+imgdata(linestarts(n):linestarts(n)+framele-1)/(frameEnd-frameStart+1);
end


% for each delay reconstract image and calculate "combined vertical diff"
step=2*pi/period;
delN=0;
for delay=mindel:Dstep:maxdel
    delN=delN+1;
    %   yda=delay:framele-1+delay;
    yda=0:framele-1;
    yda=yda+delay;
    xii=uint16((1-cos(yda*step))*255/2+1);
    yii=uint16(yda*239.99999/(framele-1)+0.5);
    ind=sub2ind_no_error_check([256,256],xii,yii);
    xii=[];
    yii=[];
    yda=[];
    fi(ind)=gi;
    fii=reshape(fi,256,256)';
    d1(delN)=mean(mean(abs(diff(fii,1,1))));
    ind=[];
    fii=zeros(256,256);
end

%find minimum and assume it is best shift
[A,I]=min(d1);
del=mindel+(I-1)*Dstep;

% metod # 2 - find min assuming that function is symmetrical
branch1=d1(1:I);
branch2=d1(I:end);

le=length(d1);
shiftb1=mindel+(0:I-1)*Dstep;
shiftb2=mindel+(I-1:le-1)*Dstep;

minD=A;
maxD=min(max(branch1),max(branch2));

q=linspace(minD+0.05*(maxD-minD),maxD,1000);

sh1=interp1(branch1,shiftb1,q);
sh2=interp1(branch2,shiftb2,q);

del2=(sh1+sh2)/2;
mean_del2=mean(del2);

% plot some graphs
%  figure;
%  plot(del2)
%  
%  figure;
%  plot(mindel:Dstep:maxdel,d1);
%  grid;
%
return;