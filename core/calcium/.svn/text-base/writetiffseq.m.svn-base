function b = writetiffseq(array, pathname, typestr, len)
%WRITETIFFSEQ
%  B = WRITETIFFSEQ(ARRAY, PATHNAME, TYPESTR, LEN)
%
% see WRITETIFF

if nargin < 4;    len = 	1000; end;

if ~exist(pathname)
    mkdir(pathname);
end;

[pathstr,name,ext,versn]=fileparts(pathname);

[ny,nx,nfr]=size(array);
nChunks = ceil(nfr/len);

for iChunk = 1:nChunks
    
    start = (iChunk-1)*len + 1;    stop = min(start + len - 1,nfr);
    newfn = sprintf('%s_%06d.tif',name,iChunk);
    bb(iChunk) = writetiff(array(:,:,start:stop),fullfile(pathstr,name,newfn),typestr);
    
end;

b = all(bb);

return;

