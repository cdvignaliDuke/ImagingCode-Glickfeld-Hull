function [hh,delta] = offsetbar(xx,yy,delta,os,varargin)
%TCOFFSETPLOT Plot time courses with offset
% h = tcOffsetPlot(yy) plots time courses
% h = tcOffsetPlot(xx,yy) specifies x axes. Default is sample #.
% h = tcOffsetPlot(xx,yy,delta) specifies offset between traces.
%  Default is three times the average standard deviation of the traces.
% h = tcOffsetPlot(xx,yy,delta, ...) additional arguments to plot

if nargin < 2
    yy = xx;
    [nrows,ncols]=size(yy);
    xx = [0:nrows-1];
end

if nargin < 3 | isempty(delta)
    if length(yy)==1
        delta = yy;
        yy = xx;
        [nrows,ncols]=size(yy);
        xx = [0:nrows-1];
    else
        delta = double(5*nanmean(nanstd(yy)));
    end
end

if nargin < 4
    os = 0;
end

if any(size(yy)==1)
    nrows = length(yy);
    ncols = 1;
    yy = yy(:);
else
    [nrows,ncols] = size(yy);
end;

offsets = 0:delta:(ncols-1)*delta;
args = {varargin{:}};
% args{end+1}='marker';
% args{end+1}='none';

hh=[];
for icol = 1:ncols
    this = yy(:,icol);
    sel = find(this);
    if length(sel)
        hh(icol) = bar(xx(sel),this(sel)+offsets(icol)+os,args{:});
        set(hh(icol),'BaseValue',offsets(icol)+os);
        set(hh(icol),'ShowBaseline','off','BarWidth',1);
    end    
    hold on;

end;

if ~length(varargin)
    set(hh,'edgecolor','k','facecolor','k');
end

if ncols >1 
    ylim([-delta delta*(ncols)]);
end;

box off;

return;