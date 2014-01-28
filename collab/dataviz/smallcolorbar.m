function h = smallcolorbar(varargin)
%SMALLCOLORBAR
%H = SMALLCOLORBAR(VARARGIN)

h = colorbar(varargin{:});
pos = get(h,'pos');
newpos = pos;
newpos(3) = pos(3)/2;
newpos(4) = pos(4)/4;
set(h,'pos',newpos);   

return;
