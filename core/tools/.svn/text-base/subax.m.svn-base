function ax = subax(rows,cols,ind,delta)
%SUBAX Create axes in tiled positions. Faster version of SUBPLOT.
% AX = SUBAX(ROWS,COLS,IND)
% AX = SUBAX(ROWS,COLS,IND,SPACE) where SPACE is fractional spacing between
% axes (default 0.25). 
%
% See also SUBPLOT

if nargin < 4
    delta = .01; 
end

fig = gcf;
children = get(fig,'children');
ax = children(find(strcmp(get(children,'type'),'axes')));

if length(ax)
    cax = ax(end);
else
    cax = gca;
    set(cax,'pos',[0 0 1 1]);
    set(cax,'Visible','off','color','none')
end

% Create and plot into axes
% pos = get(cax,'Position');

deltax = delta; 
deltay = delta; 

width = (1-(cols+1)*deltax)/cols;
height = (1-(rows+1)*deltay)/rows;
xstep = width+deltax;
ystep = height+deltay;

xlim = zeros([rows cols 2]);
ylim = zeros([rows cols 2]);
BigAxHV = get(cax,'HandleVisibility');
BigAxParent = get(cax,'Parent');

% fill figure rowwise
[j,i]=ind2sub([cols,rows],ind);

axPos = [deltax + (j-1)*xstep, deltay+(rows-i)*ystep ...
         width height];
         
% findax = findobj(fig,'Type','axes','Position',axPos);
% if isempty(findax),
ax = axes('Position',axPos,'HandleVisibility',BigAxHV,'parent',BigAxParent);
set(ax,'visible','on');
% else
%   ax = findax(1);
% end
    
set(ax,'xlimmode','auto','ylimmode','auto','xgrid','off','ygrid','off');

return
