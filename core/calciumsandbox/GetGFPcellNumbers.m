
GFPref_blue_fname = 'J:\data\mouse090316\Images\mouse090316_1_z150_920nm_gfp_ch01.tif';
GFPref_green_fname = 'J:\data\mouse090316\Images\mouse090316_1_z150_920nm_gfp_ch00.tif';
avgimg_fname='J:\data\mouse090316\rtif\rvstim4_avg.tif';
labelimg_fname='J:\data\mouse090316\labelimg\rvstim4_labelimg.mat';

lp_filter=1;
hp_filter=5;


GFPref_blue=imread(GFPref_blue_fname);
GFPref_green=imread(GFPref_green_fname);
avg_img=imread(avgimg_fname);
avg_img=double(avg_img);
avg_img=avg_img./(max(avg_img(:)));
load(labelimg_fname);

% align ref images to avgimg
output = dftregistration(fft2(avg_img),fft2(double(GFPref_green)));
shifted_ref_blue=circshift(GFPref_blue(:,:,1),output([3,4]));
shifted_ref_green=circshift(GFPref_green(:,:,1),output([3,4]));

% spatial filtering of GFP image
ker_lp=fspecial('gaussian',round(lp_filter*5),lp_filter);
ker_hp=fspecial('gaussian',round(hp_filter*5),hp_filter);
blue=filter2(ker_lp,double(shifted_ref_blue));
blue=blue-filter2(ker_hp,blue);
blue=blue./max(blue(:));

%overlay images
im0(:,:,1)=shifted_ref_blue;
im0(:,:,2)=shifted_ref_green;
im0(:,:,3)=shifted_ref_blue;

% red: GFP, cyan: cell masks
bwlabel=(labelimg>0);
im(:,:,1)=blue*5;
im(:,:,2)=bwlabel;
im(:,:,3)=bwlabel;

% get cell numbers by clicking
figure;
subplot(1,2,2);
imshow(im0);
subplot(1,2,1);
imshow(im);
[xi,yi,p]=impixel;
Ncells=length(xi);
cellNumbers=zeros(Ncells,1);
for i=1:Ncells
    cellNumbers(i)=labelimg(yi(i),xi(i));
end

sort(cellNumbers)

im2=imShade(avg_img, logical(labelimg));
im2R=im2(:,:,1);
im2G=im2(:,:,2);
im2B=im2(:,:,3);
for i=1:Ncells
    im2R(find(labelimg==cellNumbers(i)))=255;
    im2G(find(labelimg==cellNumbers(i)))=0;
    im2B(find(labelimg==cellNumbers(i)))=0;
end
im2(:,:,1)=im2R;
im2(:,:,2)=im2G;
im2(:,:,3)=im2B;

figure
subplot(2,2,1);
imshow(im);
subplot(2,2,2);
imshow(im0);
subplot(2,2,3);
imshow(im2);

xy=getCellCoordinate(labelimg);
for i=1:length(cellNumbers)
    text(xy(1,cellNumbers(i))+5,xy(2,cellNumbers(i)),num2str(cellNumbers(i)));
end

