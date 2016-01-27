mouse = 'AW41';
subN = '641';
date = '151203';
% time = '1027';

% cd('Y:\home\andrew\Behavior\Data')
% load(['data-i' subN '-' date '-' time]);

% cd(['Y:\home\ashley\Data\' subN '\widefield imaging\' date '_' subN '\' experimentname]);
cd(['Y:\home\ashley\Analysis\' mouse '\widefield imaging\' date]);

%% load retinotopy for mask making
% experimentname = 'msFRONT_expt3_dFoverF.tif';
experimentname = 'msFRONT_ret6pos_Substack (7-12).tif';

% cd(['Z:\data\AW36\widefield imaging\' date '_' mouse '\' experimentname]);
data = double(readtiff(experimentname));

%% select ROI
figure; imagesq(data(:,:,4)); colormap gray

bwout = imCellEditInteractive(squeeze(data(:,:,4)));
maskAreas = bwlabel(bwout);

V1poly = impoly;
V1mask = zeros(size(data,1),size(data,2),1);
V1mask(:,:,1) = createMask(V1poly);

%% load experiment data
experiments = [3 4]
plotcolors = {'g';'b'};
figure;
for iexp = 1:length(experiments)
    
clear data
clear V1_TC
experimentname = ['msFRONT_expt' num2str(experiments(iexp)) '_dFoverF.tif'];
data = double(readtiff(experimentname));


%% get timecourse from ROIs
tempTC = stackGetTimeCourses(data,maskAreas);
tempTC = bsxfun(@minus,tempTC,min(tempTC,[],1));
V1_TC = tempTC;
plot(V1_TC)
hold on
end