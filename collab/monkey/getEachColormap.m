


function [colormapR, colormapG, colormapB]=getEachColormap(dF, percentile)

if nargin<2
    percentile = 99.99;
end

[Nx,Ny,Nstim]=size(dF);

redON=dF(:,:,8);
redOFF=dF(:,:,4);
greenON=dF(:,:,2);
greenOFF=dF(:,:,6);
blue=dF(:,:,9)*2;
yellow=dF(:,:,10);

top1pRG=getPercentile([redON, redOFF,greenON,greenOFF],percentile);
top1pBY=getPercentile([blue,yellow],percentile);

nredON=imScale(redON,[0,top1pRG]);
nredOFF=imScale(redOFF,[0,top1pRG]);
ngreenON=imScale(greenON,[0,top1pRG]);
ngreenOFF=imScale(greenOFF,[0,top1pRG]);
nblue=imScale(blue,[0,top1pBY]);
nyellow=imScale(yellow,[0,top1pBY]);


colormapR=zeros(Nx,Ny,3);
colormapR(:,:,1)=nredON;
colormapR(:,:,2)=nredOFF;
colormapR(:,:,3)=nredOFF;

colormapG=zeros(Nx,Ny,3);
colormapG(:,:,1)=ngreenOFF;
colormapG(:,:,2)=ngreenON;
colormapG(:,:,3)=ngreenOFF;

colormapB=zeros(Nx,Ny,3);
colormapB(:,:,1)=nyellow;
colormapB(:,:,2)=nyellow;
colormapB(:,:,3)=nblue;

