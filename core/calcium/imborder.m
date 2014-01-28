function mask = imborder(siz,width,val)
%IMBORDER
%MASK = IMBORDER(SIZ,WIDTH,VAL)

if nargin < 2
    width = 5;
end

if nargin < 3
    val = 0;
end

mask = ones(siz);

mask(1:width,:)=val;
mask(end-width+1:end,:)=val;
mask(:,1:width)=val;
mask(:,end-width+1:end)=val;

return;
