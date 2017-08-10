function plotTC_CRP(tc_avg, sm, reg_mean, ICuse, frGetFrameRate, out_dir, saveData)
colord=[         0         0    1.0000
    0    0.4000         0
    1.0000         0         0
    0    0.7500    0.7500
    0.7500         0    0.7500
    0.8, 0.5, 0
    0         0    0.5
    0         0.85      0];


shift=0;
figure;
h1 = subplot(1,2,1);
imshow(mat2gray(reg_mean))
for i = ICuse
    hold on
    contour(sm(:,:,i),'Color',colord(mod(i-1,size(colord,1))+1,:))
    
end
hold off
ax=get(h1,'Position');
ax(1)=ax(1)-0.1; %or wathever
ax(2)=ax(2)-0.1; %or wathever
set(h1,'Position',ax);
axis image tight off

h2=subplot(1,2,2);
for i = ICuse
    plot(tc_avg(:,i)+shift,'Color',colord(mod(i-1,size(colord,1))+1,:)); hold on
    %     if ismember(i, [5,10,15,20,25,30,35,40,45,50])
    %         hline(shift);
    %     end
    shift = shift+18000;
end

ax=get(h2,'Position');
ax(3)=ax(3)*1.5; %or wathever
ax(4)=ax(4)*1.1; %or wathever
ax(1)=ax(1)-0.12; %or wathever
ax(2)=ax(2)-0.06; %or wathever
% ax(1)=ax(1)-0.1; %or wathever
set(h2,'Position',ax);
set(gca,'YTick',18000:18000:shift);
set(gca,'YTicklabel',(ICuse));

if saveData == 1
    savefig([out_dir, 'TC.fig']);
    print([out_dir, 'TC.eps'],'-depsc')
end