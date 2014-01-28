function rasterplot(xx,yy,raster,lims);
%RASTERPLOT(XX,YY,RASTER,LIMS);

if nargin < 3
    raster = xx;
end

inverted = -raster;

if nargin < 4
    lims = [prctile(inverted(:),2) 0];
else 
    lims = [-lims(2) lims(1)];
end

im = imScale(repmat(inverted',[1 1 3]),lims);

if nargin < 3
    imagesc(im);%colormap gray;
else
    imagesc(xx,yy,im);%colormap gray;
end

box off;

return;
