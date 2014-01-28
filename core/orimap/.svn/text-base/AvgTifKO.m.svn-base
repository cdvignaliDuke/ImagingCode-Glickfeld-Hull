function AvgIm=AvgTifKO(filename,period,Nrep,nbinning)

temp=add2d2dCR(double(imread(filename,'tif',1)),nbinning);
[xDim, yDim]=size(temp);
AvgIm=zeros(xDim,yDim,period);
if nargin<3
    inf=imfinfo(filename);
    le=length(inf);
else
    le=period*Nrep;
end

for j=1:period
    j
    AvgIm(:,:,j)=add2d2dCR(double(imread(filename,'tif',j)),nbinning);
end

for i=(period+1):le
    i
    AvgIm(:,:,mod(i-1,period)+1)=AvgIm(:,:,mod(i-1,period)+1)+add2d2dCR(double(imread(filename,'tif',i)),nbinning);
end

% AvgIm=(AvgIm/le)*period;
