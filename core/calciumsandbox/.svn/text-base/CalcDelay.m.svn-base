function [out,d1]=CalcDelay(currentImgData)

% for -35 to +35: reshift image (use only central part, border 40 pixels), calculate "combined vertical diff"

sz=size(currentImgData);
tstImag=zeros(sz(1),sz(2)-80);

delN=0;
mindel=-35;
maxdel=35;
del=mindel:maxdel;
for delay=del
    tstImag=zeros(sz(1),sz(2)-80);
    delN=delN+1;
    for j=1:sz(1)
        sh=delay*2*(mod(j,2)-0.5);
        tstImag(j,:)=currentImgData(j,41+sh:end-40+sh);
    end
    d1(delN)=mean(mean(abs(diff(tstImag,1,1))));
end


%find min assuming that function is symmetrical
%find minimum and assume it is best shift
[A,I]=min(d1);

branch1=d1(1:I);
branch2=d1(I:end);
branch1_delay=del(1:I);
branch2_delay=del(I:end);

minD=A;
maxD=min(max(branch1),max(branch2));

q=linspace(minD+0.05*(maxD-minD),maxD,1000);
try
    sh1=interp1(branch1,branch1_delay,q);
    sh2=interp1(branch2,branch2_delay,q);
catch
    figure;
    plot(branch1,branch1_delay);
    figure;
    plot(branch2,branch2_delay);
    return
end

del2=(sh1+sh2)/2;
mean_del2=mean(del2);

out=mean_del2;

