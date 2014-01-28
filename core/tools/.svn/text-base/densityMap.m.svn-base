function map=densityMap(X,Y,sigmaX,Xbin,Ybin)

% map=densityMap(X,Y,sigma,Xbin,Ybin)
% demonstrate 2D scatter plots as a density map
% X: data for X axis, Nx1 or 1xN
% Y: data for Y axis, Nx1 or 1xN
% sigmaX: sigma of Gaussian kernal for smoothing, in the unit of X
%           sigmaX=0, no smoothing (default)
% Xbin: binsize for X axis: default 1
% Ybin: binsize for Y axis: default 1
%
% e.g.  map=densityMap(COvalue,OSI,0,0.01,0.01);
%
% 2010. 02. 08. Kenichi Ohki

if nargin <3
    sigmaX=0;
end

if nargin <4
    Xbin=1;
end

if nargin <5
    Ybin=1;
end

N=length(X);

Xmin=min(X(:));
Xmax=max(X(:));
Ymin=min(Y(:));
Ymax=max(Y(:));

map=zeros(round((Xmax-Xmin)/Xbin)+1, round((Xmax-Xmin)/Xbin)+1);
X2=round((X-Xmin)/Xbin)+1;
Y2=round((Y-Ymin)/Ybin)+1;
for i=1:N
    map(X2(i),Y2(i))=map(X2(i),Y2(i))+1;
end

map=map/Xbin/Ybin;

if sigmaX > 0
    sigma2=sigmaX/Xbin;
    ker=fspecial('gaussian', ceil(sigma2*5),sigma2);
    map= imFilter2 (ker, map);
end

figure
imagesc(flipud(map'))
colorbar






