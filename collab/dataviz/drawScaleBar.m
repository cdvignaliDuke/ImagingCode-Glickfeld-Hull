function drawScaleBar(xlen,ylen)
% drawScaleBar(xlen,ylen)

hold on;
xl = round(xlim);
yl = round(ylim);

plot(xl(1)+[0 xlen],0.99*[1 1 ]*yl(1),'k','linewidth',3);
plot([1 1]*xl(1),yl(1)+[0 ylen],'k','linewidth',3);

text(xl(1)+xlen,yl(1),num2str(xlen),'fontweight','bold','fontsize',14,'HorizontalAlignment','left',...
                         'VerticalAlignment','top');
text(xl(1),yl(1)+ylen,num2str(ylen),'fontweight','bold','fontsize',14,'HorizontalAlignment','right',...
                         'VerticalAlignment','bottom');

return;
