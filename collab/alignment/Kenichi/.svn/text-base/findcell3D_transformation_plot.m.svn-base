
% construct the affine matirx, M, from hand-picked coordinates
% input: hand-picked coordinates from the fixed and vivo image (coordf & coordv)
%        raw data of hand-picked coordinates from the fixed and vivo image (excell)
% output: 3-by-4 affine matrix, M 
load 'D:\ReidLab rotation\DAPI_stacks\test\excell';
load 'D:\ReidLab rotation\DAPI_stacks\test\BW_fixed_z134_z301\labeled3D_centroids.mat';
load 'D:\ReidLab rotation\DAPI_stacks\test\BW_Set2_z80_z367\labeled3D_centroids.mat';
excell(:,3)=excell(:,3)+1; % excell file (from imageJ) starts from z=0, while findcell3D (matlab code) starts from z=1
i=[1:2:size(excell,1)];
j=[2:2:size(excell,1)];
coordf=excell(i,:);coordf=(coordf)';
coordf=cat(1, coordf,ones(1,size(excell,1)/2));
coordv=excell(j,:);coordv=(coordv)';
M=coordv/coordf;
em=(coordv-M*coordf).^2; 
errors=sqrt(em(1,:)+em(2,:)+em(3,:)); % errors in pixel of transformed coordinates. 

% plot transformed coordinates. Highlight hand-picked points with green.
tformed=M*centroidsf;
figure;scatter3(tformed(1,:),tformed(2,:),tformed(3,:),50,'b','o');
axis([0 600 0 600 -100 600]);
hold;
tformed2=M*coordf
scatter3(tformed2(1,:),tformed2(2,:),tformed2(3,:),50,'g','filled');
h=gcf;
set(h,'name','fixed_tformed');

% plot target coordinates (coordv). Highlight hand-picked points with green. 
figure;scatter3(centroidsv(1,:),centroidsv(2,:),centroidsv(3,:),50,'r','o');
axis([0 600 0 600 -100 600]);
hold;
scatter3(coordv(1,:),coordv(2,:),coordv(3,:),50,'g','filled');
h=gcf;
set(h,'name','vivo');



