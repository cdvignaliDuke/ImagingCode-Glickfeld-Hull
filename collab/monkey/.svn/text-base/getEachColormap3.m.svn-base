


function [nredON, nredOFF, ngreenON, ngreenOFF, nblueON, nblueOFF, nredONgreenOFF, ngreenONredOFF]=getEachColormap3(dF, percentile)


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
