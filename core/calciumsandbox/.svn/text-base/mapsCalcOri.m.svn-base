function params = mapsCalcOri(maps)
%MAPSCALCORI
% PARAMS = MAPSCALCORI(MAPS) where MAPS is array of size [ny,nx,ndir]
%
% PARAMS is set of parameters of direction (orientation) selectivity are obtained by vector averaging.
% th: preferred angle obtained by vector averaging, normalized to 0-1.
% mag: vector magnitude
% ave_change: average signal change to all directions
% max_change: max signal change to any directions
% tune: vector sum / scalar sum of singal changes to all directions, according to Bonhoeffer et al. (1995)
% this function can be used to estimate parameters of both orientation & direction selectivity.
% if input is dF, mag, ave_change and max_change will be absolute change.
% if input is ratio, mag, ave_change and max_change will be ratio change.
% th & tune will be the same, regardless the input is dF or ratio.
% when the signal change to one direction is nagative, it is replaced by zero, because tune becomes strange.
% thus, tune parameter is between 0 and 1.
%
% by Vincent Bonin based on code by Kenichi Ohki

[ny,nx,ndir]=size(maps);

params.th=zeros(ny, nx);
params.mag=zeros(ny, nx);
params.ave_change=zeros(ny, nx);
params.max_change=zeros(ny, nx);
params.tune=zeros(ny, nx);

a=zeros(ndir,1);

for x=1:ny
    for y=1:nx
        Vx=0;
        Vy=0;
        sum=0;
        % vector averageing
        for i=1:ndir
            if maps(x,y,i)<0; a(i)=0;  else a(i)=maps(x,y,i); end;
            Vx=Vx+a(i)*cos(2*(i-1)*pi/ndir);
            Vy=Vy+a(i)*sin(2*(i-1)*pi/ndir);
            sum=sum+a(i);
        end
        params.th(x,y)=atan2(Vx,Vy);
        params.mag(x,y)=(Vx^2+Vy^2)^(1/2);
        params.ave_change(x,y)=sum./ndir;
        params.max_change(x,y)=max(a);
        if sum~=0; params.tune(x,y)=params.mag(x,y)/sum; else params.tune(x,y)=0; end;
    end
end

params.th = params.th./pi/2+0.5; % normalize to 0-1.

return;