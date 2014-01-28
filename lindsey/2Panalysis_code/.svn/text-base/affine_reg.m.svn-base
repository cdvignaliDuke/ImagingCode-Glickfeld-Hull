%load up 2 images, call them AVG and target (latter is the image you want to align to)
base = 'G:\users\lindsey\analysisLG\active mice\LG27';
expt = 'LG27_avg_target';
fn= fullfile(base, [expt '.tif']);
outDir = 'G:\users\lindsey\analysisLG\active mice\LG27';
avg_target = readtiff(fn);
AVG = avg_target(:,:,1);
target = avg_target(:,:,2);

sz14 = size(AVG);

clear('input_points');
clear('base_points');

%use matlab GUI to select 3 points in image AVG and same 3 points in target
cpselect(double(AVG)/4e3,double(target)/4e3)
mytform = maketform('affine',input_points(1:3,:), base_points(1:3,:));
registered = imtransform(double(AVG),mytform,'XData',[1 sz14(2)],'YData',[1 sz14(1)]);

figure;
imagesc(target); colormap('gray'); axis image
figure; 
imagesc(registered); colormap('gray'); axis image

fn_out_FOV = fullfile(outDir, 'FOV_reg.tif');
writetiff(registered, fn_out_FOV);

intrinsic_expt = 'LG27_110126_V1ALPM';
fn_intrinsic = fullfile(base, [intrinsic_expt '.tif']);
intrinsic = readtiff(fn_intrinsic);
reg_intrinsic = imtransform(double(intrinsic),mytform,'XData',[1 sz14(2)],'YData',[1 sz14(1)]);
fn_out_intrinsic = fullfile(outDir, 'intrinsic_reg.tif');
writetiff(reg_intrinsic, fn_out_intrinsic);