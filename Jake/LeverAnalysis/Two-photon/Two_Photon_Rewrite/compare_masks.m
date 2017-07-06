function compare_masks(reg_max, reg_max2, sm, ICuse, mouse);
colord=[         0         0    1.0000
    0    0.4000         0
    1.0000         0         0
    0    0.7500    0.7500
    0.7500         0    0.7500
    0.8, 0.5, 0
    0         0    0.5
    0         0.85      0]; 

figure;
h1 = subplot(1,2,1);
imshow(mat2gray(reg_max))
for i = ICuse
    hold on;
    contour(sm(:,:,i),'Color',colord(mod(i-1,size(colord,1))+1,:));
end
hold off
ax=get(h1,'Position');
ax(1)=ax(1)-0.1; 
ax(2)=ax(2)-0.1; 
title('session 1')
%axis image tight off

subplot(1,2,2);
imshow(mat2gray(reg_max))
for i = ICuse
    hold on;
    contour(sm(:,:,i),'Color',colord(mod(i-1,size(colord,1))+1,:));
end
hold off
ax=get(h1,'Position');
ax(1)=ax(1)-0.1;
ax(2)=ax(2)-0.1; 
title('session 2')
suptitle([mouse, 'comparison of mask overlayed on max projection from session 1 and session 2'])