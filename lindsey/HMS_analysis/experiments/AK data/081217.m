%% set path name

pn = 'G:\users\lindsey\dataLG\AK data\081217_ori_1to4_1hz_psti.tif';
analpn = 'g:\users\lindsey\analysisLG\AKdata';

%% load stack in memory
stack_green = readtiff(pn);
stacks_green = single(stack_green);

stack_red = readtiff(pn2,[],'ch00',true);
stacks_red = single(stack_red);
%% file properties
on = 8;
off = 8;

%% plot green channel
img_green = mean(stacks_green,3);
img_green_adapt = imScale(stackLocalContrastAdaptation(img_green, 30, 1));
imagesq(img_green_adapt);

%% find cell masks
% cell find code from MH
bwimgcell = imCellEditInteractive(img_green,[]);
mask_cell = bwlabel(bwimgcell);
image(mask_cell);
figure;
imagesq(imShade(img_green_adapt,mask_cell>0));
save(fullfile(analpn,'081217cellmasks.mat'),'mask_cell');

%% find astrocyte masks
img_red = mean(stacks_red,3);
img_red_adapt = imScale(stackLocalContrastAdaptation(img_red, 30, 1));
bwimgastro = imCellEditInteractive(img_red_adapt,[]);
mask_astro = bwlabel(bwimgastro);
figure;
imagesq(imShade(img_green_adapt,mask_astro>0));
save(fullfile(analpn,'OriAastromasks.mat'),'mask_astro')

%% reload masks
load 'g:\users\lindsey\analysisLG\AKdata\081217cellmasks.mat';
load 'g:\users\lindsey\analysisLG\090128\OriAastromasks.mat';

%% remove astrocytes

common=(mask_astro>1).*mask_cell;
astro_id = unique(common(find(common(:)>0)));

mask_neuron = mask_cell;
for iCell = astro_id'
    mask_neuron(find(mask_cell==iCell))=0;
end
mask_neuron=bwlabel(mask_neuron>0);
figure;
imagesq(imshade(img_green_adapt, mask_neuron>0));
save(fullfile(analpn,'OriAneuronmasks.mat'),'mask_neuron')

%% cell time courses
load 'g:\users\lindsey\analysisLG\090128\OriAneuronmasks.mat';

%get time courses and remove low frequencies
allcells = max(unique(mask_cell));
timeCourses = stackGetTimeCourses(stacks_green,mask_cell);
timeCourses_lowcut = tcLowCut (timeCourses, 200, 'gaussian', 1);

%get dF/F for whole timecourse
baseline = mean(timeCourses_lowcut(off:(off+on):end,:));
dF = bsxfun(@minus,timeCourses_lowcut,baseline);
ratio_all = bsxfun(@rdivide,dF,baseline)*100;

%averaged individual cell time courses
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,1:5)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,6:10)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,11:15)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,16:20)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,21:25)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,26:30)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,31:35)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,36:40)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,41:45)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,46:50)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,51:55)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,56:60)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,61:65)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,66:70)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,71:75)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,76:80)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,81:85)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,86:90)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,91:95)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,96:100)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,101:105)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,106:110)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,111:115)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,116:120)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,121:125)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,126:130)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,131:135)));
figure;
tcOffsetPlot(tcRemoveDC(ratio_all(:,136:140)));

figure;
plot(ratio_all(:, 30));

%subplots of all averages
cycle = 1:epoch;
for iCell = 30;
    figure;
    plot(av(:, iCell),'r-');
    hold on
    for trial = 0:9;
        plot(ratio(cycle+(trial*epoch:trial*epoch), iCell),'k:');
        hold on;
    end
end

%% cell orientation tuning
Oris = [0 pi/6 2*pi/6 3*pi/6 4*pi/6 5*pi/6 6*pi/6 7*pi/6 8*pi/6 9*pi/6 10*pi/6 11*pi/6];
Response_Means = zeros(length(Oris), allcells);
Response_Norms = zeros(length(Oris), allcells);
for iCell = 1:allcells;
    amp0 = sum(ratio(7:11,iCell));
    amp30 = sum(ratio(17:21,iCell));
    amp60 = sum(ratio(27:31,iCell));
    amp90 = sum(ratio(37:41,iCell));
    amp120 = sum(ratio(47:51,iCell));
    amp150 = sum(ratio(57:61,iCell));
    amp180 = sum(ratio(67:71,iCell));
    amp210 = sum(ratio(77:81,iCell));
    amp240 = sum(ratio(87:91,iCell));
    amp270 = sum(ratio(97:101,iCell));
    amp300 = sum(ratio(107:111,iCell));
    amp330 = sum(ratio([1 117 118 119 120],iCell));
    Response = [amp0 amp30 amp60 amp90 amp120 amp150 amp180 amp210 amp240 amp270 amp300 amp330];   
    Response_Means(:,iCell) = Response;
    Response_Norms(:,iCell) = Response./max(Response);
end

%Response threshold map
mask_neuron_R = mask_neuron;
for iCell = 1:allcells;
    mask_neuron_R(find(mask_neuron==iCell)) = (max(Response_Means(:,iCell)))/5;
end
mask_neuron_Rthresh = mask_neuron_R;
mask_neuron_Rthresh(find(mask_neuron_R<5))=0;
mask_neuron_Rthresh = bwlabel(mask_neuron_Rthresh);
figure;
imagesq(mask_neuron_Rthresh);

%polar plots
for iCell = 1:10;
    t = Oris;
    r = Response_Norms(:,iCell)';
    figure;
    polar(t,r, '-o');
end

%Direction maps
mask_neuron_dir = mask_neuron_Rthresh;
for iCell = 1:allcells;
    mask_neuron_dir(find(mask_neuron_Rthresh==iCell)) = find(Response_Norms(:,iCell) == 1);
end
figure;
imagesc(mask_neuron_dir);
title('Direction Map for cells where R>0.5');

%Find magnitude of response for direction opposite to preferred for measure
%of direction selectivity
DirSelectivity = [1, allcells];
for iCell = 1:allcells;
    [I J] = find(Response_Norms(:,iCell) == 1);
    if  I < 6.5;
        DirSelectivity(1,iCell) = Response_Norms(find(Response_Norms(:,iCell) == 1)+6, iCell);
    else
        DirSelectivity(1,iCell) = Response_Norms(find(Response_Norms(:,iCell)==1)-6, iCell);
    end
end

mask_neuron_dirselect = mask_neuron_Rthresh;
for iCell = 1:allcells;
    mask_neuron_dirselect(find(mask_neuron_Rthresh==iCell)) = DirSelectivity(:,iCell);
end
figure;
imagesc(mask_neuron_dirselect);
title('Direction Selectivity Map (Opp Dir/Pref Dir) for cells where R>5');

mask_neuron_dirthresh = mask_neuron_dir;
mask_neuron_dirthresh(find(mask_neuron_dirselect>.6))=0;
figure;
imagesc(mask_neuron_dirthresh);
title('Direction Map for cells where R>5 and DS<0.6');

%Collapse direction map into orientation map
mask_neuron_ori = mask_neuron_CVthresh;
for iCell = 1:allcells
    [I J] = find(Response_Norms(:,iCell) == 1);
    if  I > 6;
        mask_neuron_ori(find(mask_neuron_CVthresh==iCell)) = find(Response_Norms(:,iCell) == 1) - 6;
    else
        mask_neuron_ori(find(mask_neuron_CVthresh==iCell)) = find(Response_Norms(:,iCell) == 1);
    end
end
figure;
imagesc(mask_neuron_ori);
orimap = [1 1 1; 0 1 0; 0 0 1; 1 0 0; 1 1 0];
colormap(orimap);
title('Orientation Map for cells where CV<0.91');

