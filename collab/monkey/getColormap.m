


function colormap=getColormap(dF, percentile)

if nargin<2
    percentile = 99.99;
end

[Nx,Ny,Nstim]=size(dF);

red=squeeze(sum(dF(:,:,[6:8]),3));
green=squeeze(sum(dF(:,:,[2:4]),3));
blue=dF(:,:,9)*2;
yellow=dF(:,:,10);

top1pRG=getPercentile([red,green],percentile);
top1pBY=getPercentile([blue,yellow],percentile);

nred=imScale(red,[0,top1pRG]);
ngreen=imScale(green,[0,top1pRG]);
nblue=imScale(blue,[0,top1pBY]);
nyellow=imScale(yellow,[0,top1pBY]);


colormap=zeros(Nx,Ny,3);
colormap(:,:,1)=nred+nyellow;
colormap(:,:,2)=ngreen+nyellow;
colormap(:,:,3)=nblue;

imshow(colormap);