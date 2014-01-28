function [h,r,p]=drawScatter(x,y,z,s,cl)
%[H,R,P]=DRAWSCATTER(X,Y,Z,S,CL)

if nargin < 4
    s = 5;
end

if nargin < 3 | isempty(z)
    h = scatter(x,y,s,'ko','filled');
else
    h = scatter(x,y,s,z,'o','filled');
    colorbar;
end

[r,p]=corrcoef(x,y);
text(.1,.99,sprintf('r=%2.2f p=%2.0g',r(1,2),p(1,2)),'units','norm','color','k');

if nargin < 5
    cl = get(gca,'clim');
    cl = max(abs(cl))*[-1 1];
end
    
set(gca,'clim',cl);
    
return;
