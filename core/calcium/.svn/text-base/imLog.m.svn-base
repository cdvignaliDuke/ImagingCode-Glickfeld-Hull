function out = imLog(img,type);
%IMLOG Log transform image intensities
%  out = imLog(img)
%  out = imLog(img, type), specified return type

if nargin < 2; type = 'double';end

img = double(img);
out = log(img);
sel = nnz(img(:));
out(find(out==-inf)) = min(img(sel));
out = imScale(out,type);

return