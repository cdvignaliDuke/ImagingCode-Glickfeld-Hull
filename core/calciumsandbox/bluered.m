function c = bluered(m,k,l)
%BLUERED Blue-red colormap
% C = BLUERED(M) where M specifies mean luminance at origin (default:1).
% C = BLUERED(M,K) where K specifies exponent (default: 1)

if nargin < 1 
    m = 1;
end

if nargin < 2
    k = 1;
end

if nargin < 3
    l = 0;
end

if length(m)< 3
    m = m*[1 1 1];
end
    
bottom = ([linspace(m(1),1,32)'.^k, linspace(m(2),0,32)'.^k, linspace(m(3),0,32)'.^k]);
top = ([linspace(l,m(1),32)'.^k, linspace(l,m(2),32)'.^k, linspace(1,m(3),32)'.^k]);
bottom(1,:)=[];
c = [top;bottom];

return;
