%% set path name

pn = 'G:\users\lindsey\analysisLG\081205\OriCRegistered.tif';
pn2 = 'G:\users\lindsey\dataLG\081205\081205OriC';
%% load stack in memory
stack_green = readtiff(pn);
stacks_green = single(stack_green);

stack_red = readtiff(pn2,[],'ch00',true);
stacks_red = single(stack_red);

%% register green channel
%img_green_target = mean(stacks_green(:,:,10:40),3);
%[OUTS, stacks_green_reg] = stackRegister(stacks_green, img_green_target);

%img_green_reg = mean(stacks_green_reg,3);
%imagesc(img_green_reg)
%% plot green channel
img_green = mean(stacks_green,3);
img_green_adapt = imScale(stackLocalContrastAdaptation(img_green, 30, 1));
imagesc(img_green_adapt)

%% remove astrocytes
img_red = mean(stacks_red,3);
img_red_adapt = imScale(stackLocalContrastAdaptation(img_red, 30, 1));
img_red_adapt(:,507:512) = 0;
imagesc(img_red_adapt);
img_cell = img_green_adapt - img_red_adapt;
imagesc(img_cell);

%% find cell masks
% template matching code from Kenichi
img_segment = img_cell;
thresh = .4; % correlation threshold
shat = 23; % size of mexican hat in pixels
scell = 12;

[mask_cell] = imFindCellsTM (img_segment, shat,thresh,scell);
figure;
imagesq(imShade(img_green_adapt,mask_cell>0));

%manual_mask_cell = imCellEditInteractive(img_cell, mask_cell, [], []);

%% cell time courses
timeCourses = stackGetTimeCourses(stacks_green,mask_cell);

%individual cell time courses
%tcOffsetPlot(timeCourses);

%averaged time courses
av = tcCycleAverage(timeCourses,80);
tcOffsetPlot(tcRemoveDC(av));
%imagesc(tcRemoveDC(av)')
xlabel('Time (s)');
ylabel('Cell #');

