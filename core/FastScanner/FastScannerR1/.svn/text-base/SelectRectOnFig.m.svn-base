function p=SelectRectOnFig
global fd;


p1=get(gca,'CurrentPoint');
rec=rbbox;
p2=get(gca,'CurrentPoint');

p1 = p1(1,1:2);              % extract x and y
p2 = p2(1,1:2);
if p1(1)<1 || p1(1)>256 || p1(2)<1 || p1(2)>256
    return
end
fd.stats.region.x1=round(max(min(p1(1),p2(1)),1));
fd.stats.region.y1=round(max(min(p1(2),p2(2)),1));
fd.stats.region.x2=round(min(max(p1(1),p2(1)),256));
fd.stats.region.y2=round(min(max(p1(2),p2(2)),256));

