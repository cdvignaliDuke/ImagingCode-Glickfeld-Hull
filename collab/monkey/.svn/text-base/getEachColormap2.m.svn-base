


function [colormapR, colormapG, colormapB,colormapRG]=getEachColormap(dF, percentile)


Sp=[115, 97, 170]'/250;
Sm=[119, 120, 49]'/250;
Lp=[236, 118, 144]'/250;
Lm=[36, 107, 53]'/250;
Mm=[184, 30, 76]'/250;
Mp=[71, 165, 94]'/250;
LpMm=[177, 59, 97]'/250;
MpLm=[41, 153, 71]'/250;

if nargin<2
    percentile = 99.99;
end

[Nx,Ny,Nstim]=size(dF);

redON=dF(:,:,8);
redOFF=dF(:,:,4);
greenON=dF(:,:,2);
greenOFF=dF(:,:,6);
redONgreenOFF=dF(:,:,7);
greenONredOFF=dF(:,:,3);
blueON=dF(:,:,9)*2;
blueOFF=dF(:,:,10);


top1pRG=getPercentile([redON, redOFF,greenON,greenOFF,redONgreenOFF,greenONredOFF],percentile);
top1pBY=getPercentile([blueON,blueOFF],percentile);

nredON=imScale(redON,[0,top1pRG]);
nredOFF=imScale(redOFF,[0,top1pRG]);
ngreenON=imScale(greenON,[0,top1pRG]);
ngreenOFF=imScale(greenOFF,[0,top1pRG]);
nredONgreenOFF=imScale(redONgreenOFF,[0,top1pRG]);
ngreenONredOFF=imScale(greenONredOFF,[0,top1pRG]);
nblueON=imScale(blueON,[0,top1pBY]);
nblueOFF=imScale(blueOFF,[0,top1pBY]);

colormapR=zeros(Nx,Ny,3);
colormapR(:,:,1)=nredON*Lp(1)+nredOFF*Lm(1);
colormapR(:,:,2)=nredON*Lp(2)+nredOFF*Lm(2);
colormapR(:,:,3)=nredON*Lp(3)+nredOFF*Lm(3);

colormapG=zeros(Nx,Ny,3);
colormapG(:,:,1)=ngreenON*Mp(1)+ngreenOFF*Mm(1);
colormapG(:,:,2)=ngreenON*Mp(2)+ngreenOFF*Mm(2);
colormapG(:,:,3)=ngreenON*Mp(3)+ngreenOFF*Mm(3);

colormapB=zeros(Nx,Ny,3);
colormapB(:,:,1)=nblueON*Sp(1)+nblueOFF*Sm(1);
colormapB(:,:,2)=nblueON*Sp(2)+nblueOFF*Sm(2);
colormapB(:,:,3)=nblueON*Sp(3)+nblueOFF*Sm(3);

colormapRG=zeros(Nx,Ny,3);
colormapRG(:,:,1)=nredONgreenOFF*LpMm(1)+ngreenONredOFF*MpLm(1);
colormapRG(:,:,2)=nredONgreenOFF*LpMm(2)+ngreenONredOFF*MpLm(2);
colormapRG(:,:,3)=nredONgreenOFF*LpMm(3)+ngreenONredOFF*MpLm(3);



