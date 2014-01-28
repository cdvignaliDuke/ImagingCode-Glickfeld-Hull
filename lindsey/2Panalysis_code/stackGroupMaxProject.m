function out = stackGroupProject(array,ratio)
%STACKGROUPPROJECT Downsample stack array by ratio
%   out = stackGroupProject(array,ratio)

[ny,nx,nframes]=size(array);

% which frames?
grpfrns = floor(nframes/ratio)*ratio; 
if grpfrns ~= nframes
    % truncate
    desfrns = 1:grpfrns;
else
    desfrns = 1:nframes;
end
ndes = length(desfrns);

siz = ndes / ratio;
%{
c = class(array);
out = zeros(ny,nx,siz,c);

for i = 1:siz
    out(:,:,i)=mean(array(:,:,(i-1)*ratio+1:i*ratio),3);
end
%}

array=reshape(max(reshape(array(:,:,desfrns),ny,nx,ratio,siz),[],3),ny,nx,siz);

out=array;

return