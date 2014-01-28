


function colormap=getConemap(dF, percentile)

if nargin<2
    percentile = 99.99;
end

[Nx,Ny,Nstim]=size(dF);

red=dF(:,:,1);
green=dF(:,:,2);
blue=dF(:,:,3);

top1p=getPercentile([red,green,blue],percentile);

nred=imScale(red,[0,top1p]);
ngreen=imScale(green,[0,top1p]);
nblue=imScale(blue,[0,top1p]);


colormap=zeros(Nx,Ny,3);
colormap(:,:,1)=nred;
colormap(:,:,2)=ngreen;
colormap(:,:,3)=nblue;

imshow(colormap);