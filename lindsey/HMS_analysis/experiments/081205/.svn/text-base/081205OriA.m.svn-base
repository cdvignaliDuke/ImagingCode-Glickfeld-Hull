%% set path name

pn = 'G:\users\lindsey\dataLG\081205\081205OriA';

%% load stack in memory
stack_green = readtiff(pn,[],'ch01',true);
stacks_green = single(stack_green);

stack_red = readtiff(pn,[],'ch00',true);
stacks_red = single(stack_red);

%% plot green channel
img_green = mean(stacks_green,3);
imagesc(img_green)

img_green_adapt = imScale(stackLocalContrastAdaptation(img_green, 30, 1));
img_clear_green = img_green_adapt;
% img_clear(1:32,:)=0;
% img_clear(512-31:512,:)=0;
% img_clear(:,1:32)=0;
% img_clear(:,512-31:512)=0;

imagesc(img_clear_green);

%% plot red channel
% to remove astrocytes, first subtract red channel
img_red = mean(stacks_red,3);
imagesc(img_red)

img_red_adapt = imScale(stackLocalContrastAdaptation(img_red, 30, 1));
img_clear_red = img_red_adapt;
imagesc(img_clear_red);

%% find cell masks
% template matching code from Kenichi
img_segment = img_clear_green;
thresh = .5; % correlation threshold
shat = 16; % size of mexican hat in pixels
scell = 10;

[mask_cell] = imFindCellsTM (img_segment, shat,thresh,scell);
figure;
imagesq(imShade(img_green_adapt,mask_cell>0));

imagesq(imShade(img_green,mask_cell>0))

%% find astrocyte masks
img_segment_astro = img_clear_red;
thresh = .6; % correlation threshold
shat = 22; % size of mexican hat in pixels
scell = 10;

[mask_astro] = imFindCellsTM (img_segment_astro, shat,thresh,scell);
figure;
imagesq(imShade(img_red_adapt,mask_astro>0));

imagesq(imShade(img_red,mask_astro>0))

%% cell time courses
timeCourses = stackGetTimeCourses(stacks_green,mask_cell);

%individual cell time courses
tcOffsetPlot(timeCourses);

%averaged time courses
av = tcCycleAverage(timeCourses,64);
tcOffsetPlot(av);
imagesc(tcRemoveDC(av)')

xlabel('Time (s)');
ylabel('Cell #');

