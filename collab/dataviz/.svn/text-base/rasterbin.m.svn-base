function outRaster=rasterbin(inRaster,n)
%outRaster=rasterbin(inRaster,n)

siz = size(inRaster);

outRaster = ((squeeze((sum(reshape(inRaster,[n,siz(1)/n,siz(2:end)]),1)))));

return;
