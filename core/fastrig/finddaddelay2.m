function [period,mean_del2]=finddaddelay2(frameStart,frameEnd,mindel,maxdel,Dstep,width,tcrossup,chIm)
%[period,mean_del2]=finddaddelay2(frameStart,frameEnd,mindel,maxdel,Dstep,width,tcrossup,chIm)
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
% 07/22/08 SY created from finddaddelay, performance optimization,
%          streamlined

fi=single(zeros(width,256));
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

gi=single(zeros(1,framele));

for n=N+1-frameStart
    gi=gi+chIm(linestarts(n):linestarts(n)+framele-1);
end

% for each delay reconstract image and calculate "combined vertical diff"
w = warning('off');
step=2*pi/period;
delN=0;
ydb0=0:framele-1;
ydb=(ydb0+mindel-Dstep)*239.99999/(framele-1)+0.5-1;
cosydb=cos(ydb0*step);
sinydb=sin(ydb0*step);

for delay=mindel:Dstep:maxdel
    delN=delN+1;
    ydb=ydb+Dstep*239.99999/(framele-1);
    
    sindelayA=sin(delay*step)*(width-1)/2;
    cosdelayA=cos(delay*step)*(width-1)/2;
    widthA=(width-1)/2+1;

    ind=uint32(ydb)*width;
%    cos(a+b)=cos(a)*cos(b) - sin(a)*sin(b)
    ind=ind+uint32(widthA-cosydb*cosdelayA+sinydb*sindelayA);
    
%    sindelay=sin(delay*step);
%    cosdelay=cos(delay*step);   
%    ind2=uint32(ydb)*width;
%    ind=ind+uint32((1-cosydb*cosdelay+sinydb*sindelay)*(width-1)/2+1);

    fi(ind)=gi;
    fii=reshape(fi,width,256);
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
%{
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
%}

% find min interpolating two point around min point
% interpolate if possible
if I>1 && I<length(d1)
    mean_del2=mindel+(I-1)*Dstep + 0.5*(d1(I-1)-d1(I+1))*Dstep/max(d1(I-1)-minD,d1(I+1)-minD);
else
    if I==1
        [period,mean_del2]=finddaddelay2(frameStart,frameEnd,mindel-5,mindel+5,Dstep,width,tcrossup,chIm);
    end
    if I==length(d1)
        [period,mean_del2]=finddaddelay2(frameStart,frameEnd,maxdel-5,maxdel+5,Dstep,width,tcrossup,chIm);
    end
        
%    mean_del2 = mindel +(I-1)*Dstep ;
end
