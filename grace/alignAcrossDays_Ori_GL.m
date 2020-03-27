clear all
clear global
%% 
mouse = 'i1316';
day2 = '200108';
day3 = '200109';
day4 = '200110';
ImgFolder = strvcat('003');
ref_date = '200106';
ref_run = strvcat('003');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
ref_str = catRunName(ref_run, size(ref_run,1));
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P';

%% load data
oriTuning_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_oriTuningAndFits.mat']));
oriTuning_D2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_oriTuningAndFits.mat']));
oriTuning_D3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_oriTuningAndFits.mat']));
oriTuning_D4 = load(fullfile(fnout, [day4 '_' mouse], [day4 '_' mouse '_' run_str], [day4 '_' mouse '_' run_str '_oriTuningAndFits.mat']));

TCs_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_TCs.mat']));
TCs_D2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_TCs.mat']));
TCs_D3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_TCs.mat']));
TCs_D4 = load(fullfile(fnout, [day4 '_' mouse], [day4 '_' mouse '_' run_str], [day4 '_' mouse '_' run_str '_TCs.mat']));

%% define variables
goodfit_D1 = find(oriTuning_D1.fitReliability<22.5);
goodfit_D2 = find(oriTuning_D2.fitReliability<22.5);
goodfit_D3 = find(oriTuning_D3.fitReliability<22.5);
goodfit_D4 = find(oriTuning_D4.fitReliability<22.5);
gdft_D1 = length(goodfit_D1);
gdft_D2 = length(goodfit_D2);
gdft_D3 = length(goodfit_D3);
gdft_D4 = length(goodfit_D4);
goodfit = intersect(goodfit_D1,goodfit_D2);
goodfit2 = intersect(goodfit, goodfit_D3);
goodfit3 = intersect(goodfit2, goodfit_D4);
goodfit3to1 = intersect(goodfit_D1, goodfit_D3);
goodfit4to1 = intersect(goodfit_D1, goodfit_D4);
nCells = size(goodfit,2);

[maxResp_D1 prefOri_D1] = max(squeeze(oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,:)),[],1);
[maxResp_D2 prefOri_D2] = max(squeeze(oriTuning_D2.vonMisesFitAllCellsAllBoots(:,1,:)),[],1);
[maxResp_D3 prefOri_D3] = max(squeeze(oriTuning_D3.vonMisesFitAllCellsAllBoots(:,1,:)),[],1);
[maxResp_D4 prefOri_D4] = max(squeeze(oriTuning_D4.vonMisesFitAllCellsAllBoots(:,1,:)),[],1);

% day 1
load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_input.mat']));
tGratingDir1 = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs1 = unique(tGratingDir1);
nDir1 = length(dirs1);
nOn1 = input.nScansOn;
nOff1 = input.nScansOff;
nFrames1 = nOn1+nOff1;
npSub_tc1 = TCs_D1.npSub_tc;
nCells1 = size(npSub_tc1,2);
nTrials1 = size(tGratingDir1,2);

trial_tc1 = reshape(npSub_tc1,[nFrames1 nTrials1 nCells1]);
trial_f1 = mean(trial_tc1(nOff1/2:nOff1,:,:),1);
trial_dfof1 = (trial_tc1-trial_f1)./trial_f1;

% day 2
load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_input.mat']));
tGratingDir2 = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs2 = unique(tGratingDir1);
nDir2 = length(dirs2);
nOn2 = input.nScansOn;
nOff2 = input.nScansOff;
nFrames2 = nOn2+nOff2;
npSub_tc2 = TCs_D2.npSub_tc;
nCells2 = size(npSub_tc2,2);
nTrials2 = size(tGratingDir2,2);

trial_tc2 = reshape(npSub_tc2,[nFrames2 nTrials2 nCells2]);
trial_f2 = mean(trial_tc2(nOff2/2:nOff2,:,:),1);
trial_dfof2 = (trial_tc2-trial_f2)./trial_f2;

% day 3
load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_input.mat']));
tGratingDir3 = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs3 = unique(tGratingDir1);
nDir3 = length(dirs3);
nOn3 = input.nScansOn;
nOff3 = input.nScansOff;
nFrames3 = nOn3+nOff3;
npSub_tc3 = TCs_D3.npSub_tc;
nCells3 = size(npSub_tc3,2);
nTrials3 = size(tGratingDir3,2);

trial_tc3 = reshape(npSub_tc3,[nFrames3 nTrials3 nCells3]);
trial_f3 = mean(trial_tc3(nOff3/2:nOff3,:,:),1);
trial_dfof3 = (trial_tc3-trial_f3)./trial_f3;

% day 4
load(fullfile(fnout, [day4 '_' mouse], [day4 '_' mouse '_' run_str], [day4 '_' mouse '_' run_str '_input.mat']));
tGratingDir4 = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs4 = unique(tGratingDir1);
nDir4 = length(dirs4);
nOn4 = input.nScansOn;
nOff4 = input.nScansOff;
nFrames4 = nOn4+nOff4;
npSub_tc4 = TCs_D4.npSub_tc;
nCells4 = size(npSub_tc4,2);
nTrials4 = size(tGratingDir4,2);

trial_tc4 = reshape(npSub_tc4,[nFrames4 nTrials4 nCells4]);
trial_f4 = mean(trial_tc4(nOff4/2:nOff4,:,:),1);
trial_dfof4 = (trial_tc4-trial_f4)./trial_f4;

%% original plots 
figure; 
% subplot(2,2,1)
scatter(maxResp_D1(goodfit),maxResp_D2(goodfit))
axis square
xlim([0 max(max(maxResp_D1(goodfit),maxResp_D2(goodfit)))+0.05])
ylim([0 max(max(maxResp_D1(goodfit),maxResp_D2(goodfit)))+0.05])
refline(1,0)
xlabel('Day 1 max dF/F')
ylabel('Day 2 max dF/F')

prefOri_diff = abs(prefOri_D1(goodfit)-prefOri_D2(goodfit));
% prefOri_diff(find(prefOri_diff>90)) = 180-prefOri_diff(find(prefOri_diff>90));
% subplot(2,2,2)
% hist(prefOri_diff)
% xlabel('D1 vs D2 Pref ori')
% ylabel('Number of cells')
% 
% subplot(2,2,3)
% scatter(maxResp_D1(goodfit),prefOri_diff)
% xlabel('Day 1 max dF/F')
% ylabel('D1 vs D2 Pref ori')
% axis square
% 
% maxResp_diff = abs((maxResp_D1(goodfit)-maxResp_D2(goodfit))./(maxResp_D1(goodfit)+maxResp_D2(goodfit)));
% subplot(2,2,4)
% scatter(maxResp_diff,prefOri_diff)
% xlabel('D1 vs D2 max dF/F')
% ylabel('D1 vs D2 Pref ori')
% axis square
% 
% suptitle([mouse ' ' ref_date ' vs ' day2 '- n = ' num2str(length(goodfit)) ' cells'])
% print(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_compAcrossDaysTuning.pdf']),'-dpdf', '-bestfit')
% save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_compAcrossDaysTuning.mat']), 'goodfit', 'prefOri_D1', 'maxResp_D1', 'prefOri_D2', 'maxResp_D2', 'prefOri_diff', 'maxResp_diff');

%% Identify cells that are significantly responsive to at least 1 direction

% Define analysis windows
resp_wind1 = nOff1+1:nOff1+nOn1;
base_wind1 = 1+nOff1-nOn1:nOff1;

resp_wind2 = nOff2+1:nOff2+nOn2;
base_wind2 = 1+nOff2-nOn2:nOff2;

resp_wind3 = nOff3+1:nOff3+nOn3;
base_wind3 = 1+nOff3-nOn3:nOff3;

resp_wind4 = nOff4+1:nOff4+nOn4;
base_wind4 = 1+nOff4-nOn4:nOff4;

% Create two matrices
dfof_resp1 = squeeze(mean(trial_dfof1(resp_wind1,:,:),1));
dfof_base1 = squeeze(mean(trial_dfof1(base_wind1,:,:),1));
dfof_subtract1 = dfof_resp1 - dfof_base1;

dfof_resp2 = squeeze(mean(trial_dfof2(resp_wind2,:,:),1));
dfof_base2 = squeeze(mean(trial_dfof2(base_wind2,:,:),1));
dfof_subtract2 = dfof_resp2 - dfof_base2;

dfof_resp3 = squeeze(mean(trial_dfof3(resp_wind3,:,:),1));
dfof_base3 = squeeze(mean(trial_dfof3(base_wind3,:,:),1));
dfof_subtract3 = dfof_resp3 - dfof_base3;

dfof_resp4 = squeeze(mean(trial_dfof4(resp_wind4,:,:),1));
dfof_base4 = squeeze(mean(trial_dfof4(base_wind4,:,:),1));
dfof_subtract4 = dfof_resp4 - dfof_base4;

% ttest for response significance
% day 1
dfof_dir1 = zeros(nDir1, nCells1, 2);
h1 = zeros(nDir1, nCells1);
p1 = zeros(nDir1, nCells1);
for idir = 1:nDir1
    ind = find(tGratingDir1==dirs1(idir));
    x1 = dfof_base1(ind,:);
    y1 = dfof_resp1(ind,:);
    [h1(idir,:),p1(idir,:)] = ttest(x1,y1,'dim',1,'Alpha',0.05./(nDir1-1),'tail','left');
    dfof_dir1(idir,:,1) = mean(y1-x1,1);
    dfof_dir1(idir,:,2) = std(y1-x1,[],1)./sqrt(length(ind));
end
h1sum = zeros(1, nCells1);
h1over1 = zeros(1, nCells1);
h1sum = sum(h1(:,:));
h1dirsum = sum(h1);
for iCell = 1:nCells1
    h1over1(:,iCell) = h1sum(:,iCell) >= 1;
end
sig_resp_cells1 = sum(h1over1(:,:));

h_1 = h1';
h1sum_dirs = zeros(1,nDir1);
h1sum_dirs = sum(h_1(:,:));

% day 2
dfof_dir1 = zeros(nDir2, nCells2, 2);
h2 = zeros(nDir2, nCells2);
p2 = zeros(nDir2, nCells2);
for idir = 1:nDir2
    ind = find(tGratingDir2==dirs2(idir));
    x2 = dfof_base2(ind,:);
    y2 = dfof_resp2(ind,:);
    [h2(idir,:),p2(idir,:)] = ttest(x2,y2,'dim',1,'Alpha',0.05./(nDir1-1),'tail','left');
    dfof_dir2(idir,:,1) = mean(y2-x2,1);
    dfof_dir2(idir,:,2) = std(y2-x2,[],1)./sqrt(length(ind));
end
h2sum = zeros(1, nCells2);
h2over1 = zeros(1, nCells2);
h2sum = sum(h2(:,:));
h2dirsum = sum(h2);
for iCell = 1:nCells2
    h2over1(:,iCell) = h2sum(:,iCell) >= 1;
end
sig_resp_cells2 = sum(h2over1(:,:));

h_2 = h2';
h2sum_dirs = zeros(1,nDir2);
h2sum_dirs = sum(h_2(:,:));

% day 3
dfof_dir3 = zeros(nDir3, nCells3, 2);
h3 = zeros(nDir3, nCells3);
p3 = zeros(nDir3, nCells3);
for idir = 1:nDir3
    ind = find(tGratingDir3==dirs3(idir));
    x3 = dfof_base3(ind,:);
    y3 = dfof_resp3(ind,:);
    [h3(idir,:),p3(idir,:)] = ttest(x3,y3,'dim',1,'Alpha',0.05./(nDir1-1),'tail','left');
    dfof_dir3(idir,:,1) = mean(y3-x3,1);
    dfof_dir3(idir,:,2) = std(y3-x3,[],1)./sqrt(length(ind));
end
h3sum = zeros(1, nCells3);
h3over1 = zeros(1, nCells3);
h3sum = sum(h3(:,:));
h3dirsum = sum(h3);
for iCell = 1:nCells3
    h3over1(:,iCell) = h3sum(:,iCell) >= 1;
end
sig_resp_cells3 = sum(h3over1(:,:));

h_3 = h3';
h3sum_dirs = zeros(1,nDir3);
h3sum_dirs = sum(h_3(:,:));

% day 4
dfof_dir4 = zeros(nDir4, nCells4, 2);
h4 = zeros(nDir4, nCells4);
p4 = zeros(nDir4, nCells4);
for idir = 1:nDir4
    ind = find(tGratingDir4==dirs4(idir));
    x4 = dfof_base4(ind,:);
    y4 = dfof_resp4(ind,:);
    [h4(idir,:),p4(idir,:)] = ttest(x4,y4,'dim',1,'Alpha',0.05./(nDir1-1),'tail','left');
    dfof_dir4(idir,:,1) = mean(y4-x4,1);
    dfof_dir4(idir,:,2) = std(y4-x4,[],1)./sqrt(length(ind));
end
h4sum = zeros(1, nCells4);
h4over1 = zeros(1, nCells4);
h4sum = sum(h4(:,:));
h4dirsum = sum(h4);
for iCell = 1:nCells4
    h4over1(:,iCell) = h4sum(:,iCell) >= 1;
end
sig_resp_cells4 = sum(h4over1(:,:));

h_4 = h4';
h4sum_dirs = zeros(1,nDir4);
h4sum_dirs = sum(h_4(:,:));

intersect1 = intersect(sig_resp_cells1,sig_resp_cells2);
intersect2 = intersect(intersect1, sig_resp_cells3);
all_sig_cells = intersect(intersect2, sig_resp_cells4);

figure; 
x = [1 2 3 4];
y = [sig_resp_cells1 sig_resp_cells2 sig_resp_cells3 sig_resp_cells4];
cellnames = {'Day 1'; 'Day 2'; 'Day 3'; 'Day 4';};
bar(x,y);
set(gca,'xticklabel',cellnames)
ylim([0 nCells1])
title('Significantly Responsive Cells')
ylabel('Number of cells')
% imagesc(h');
% xlabel('Direction')
% ylabel('Cells')
% figAxForm
if exist(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays'])
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\sig_resp_cells'],'-dpdf') 

else mkdir(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays'])
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\sig_resp_cells'],'-dpdf') 
end

figure;
x = [1 2 3 4];
vals = [sig_resp_cells1 gdft_D1; sig_resp_cells2 gdft_D2; sig_resp_cells3 gdft_D3; sig_resp_cells4 gdft_D4]; 
cellnames = {'Day 1'; 'Day 2'; 'Day 3'; 'Day 4';};
b11 = bar(x,vals);
b11(1).FaceColor = [0.6953 0.1328 0.1328];
b11(2).FaceColor = [0.9102 0.5859 0.4766];
set(gca,'xticklabel',cellnames)
ylim([0 nCells1+2])
% title('Significantly Responsive Cells vs Goodfit Cells')
ylabel('Number of cells')
legend({'Sig Cells', 'Reliably Fit Cells'},'Location','northeast')
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\sig_resp_vs_goodfit'],'-dpdf') 

figure;
b = [h1dirsum; h2dirsum; h3dirsum; h4dirsum]';
b1 = bar(b,'stacked');
b1(1).FaceColor = [0.0977 0.0977 0.5];
b1(2).FaceColor = [0.2334 0.4678 0.700];
b1(3).FaceColor = [0.6558 0.8238 0.8984];
b1(4).FaceColor = [0.9375 0.9708 1.0000];
legend({'Day 1', 'Day 2', 'Day 3', 'Day 4'},'Location','northeast')
xlabel('Cell #')
ylabel('nDirs')
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\nDirs_eachCell'],'-dpdf') 

figure;
b1 = bar(h1dirsum);
b1(1).FaceColor = [0.2334 0.4678 0.700];
xlabel('Cell #')
ylabel('nDirs')
legend({'Day 1', 'Day 2', 'Day 3', 'Day 4'},'Location','northeast')
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\nDirs_eachCellD1'],'-dpdf') 

figure;
y = [h1sum_dirs; h2sum_dirs; h3sum_dirs; h4sum_dirs]';
y1 = bar(y,'stacked');
y1(1).FaceColor = [0.0977 0.0977 0.5];
y1(2).FaceColor = [0.2334 0.4678 0.700];
y1(3).FaceColor = [0.6558 0.8238 0.8984];
y1(4).FaceColor = [0.9375 0.9708 1.0000];
% cellnames = {'0';'22.5';'45';'67.5';'90';'112.5';'135';'157.5';'180';'202.5';'225';'247.5';'270';'292.5';'315';'337.5'};
cellnames = {'0';'45';'90';'135';'180';'225';'270';'315'};
xticks([1:2:16])
set(gca,'xticklabel',cellnames)
xlabel('Direction')
ylabel('Significantly Responsive Cells')
legend({'Day 1', 'Day 2', 'Day 3', 'Day 4'},'Location','northeast')
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\sig_resp_16dir'],'-dpdf') 

%% Preferred Direction Plots
% with goodfit cells
prefOri_D21 = find(((prefOri_D2(goodfit3))<(prefOri_D1(goodfit3)+10)) & ((prefOri_D2(goodfit3))> (prefOri_D1(goodfit3)-10)));
prefOri_D31 = find(((prefOri_D3(goodfit3))<(prefOri_D1(goodfit3)+10)) & ((prefOri_D3(goodfit3))> (prefOri_D1(goodfit3)-10)));
prefOri_D41 = find(((prefOri_D4(goodfit3))<(prefOri_D1(goodfit3)+10)) & ((prefOri_D4(goodfit3))> (prefOri_D1(goodfit3)-10)));

% with all cells
prefOri_D2_1 = find(((prefOri_D2)<(prefOri_D1+10)) & ((prefOri_D2)> (prefOri_D1-10)));
prefOri_D3_1 = find(((prefOri_D3)<(prefOri_D1+10)) & ((prefOri_D3)> (prefOri_D1-10)));
prefOri_D4_1 = find(((prefOri_D4)<(prefOri_D1+10)) & ((prefOri_D4)> (prefOri_D1-10)));

figure;
subplot(1,2,1); 
x = [1 2 3];
y = [size(prefOri_D21,2)/size(goodfit3,2) size(prefOri_D31,2)/size(goodfit3,2) size(prefOri_D41,2)/size(goodfit3,2)]; 
bar(x,y,'facecolor',[0 0.3906 0]);
cellnames = {'Day 2';'Day 3';'Day 4'};
set(gca,'xticklabel',cellnames)
ylabel({'fraction of reliably fit cells';'maintaining their initial tuning'})
axis square
subplot(1,2,2); 
x = [1 2 3];
y = [size(prefOri_D2_1,2)/nCells1 size(prefOri_D3_1,2)/nCells1 size(prefOri_D4_1,2)/nCells1]; 
bar(x,y,'facecolor',[0.0900 0.0977 0.5]);
cellnames = {'Day 2';'Day 3';'Day 4'};
set(gca,'xticklabel',cellnames)
ylabel({'fraction of all cells maintaining'; 'their initial tuning'})
ylim([0 max(y)+0.1])
axis square
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\maintained_tuning'],'-dpdf') 

%% dfof traces
% trace max dfof of each cell over 4 days
figure;
subplot(2,1,1);
nCells = size(goodfit,2);
hold on
x = [1 2 3 4];
xticks([1 2 3 4]);
cellnames = {'Day1'; 'Day2'; 'Day3'; 'Day4';};
for iCell = 1:nCells
    maxResp_all = [maxResp_D1(iCell) maxResp_D2(iCell) maxResp_D3(iCell) maxResp_D4(iCell)];
    plot(x,maxResp_all)
    set(gca,'xticklabel',cellnames)
    ylabel('max dF/F')
    title([mouse ' cells goodfit on days 1 and 2'])
end

subplot(2,1,2);
nCells2 = size(goodfit3,2);
hold on
x = [1 2 3 4];
xticks([1 2 3 4]);
cellnames = {'Day1'; 'Day2'; 'Day3'; 'Day4';};
for iCell = 1:nCells2
    maxResp_all = [maxResp_D1(iCell) maxResp_D2(iCell) maxResp_D3(iCell) maxResp_D4(iCell)];
    plot(x,maxResp_all)
    set(gca,'xticklabel',cellnames)
    ylabel('max dF/F')
    title([mouse ' cells goodfit on all days'])
end

mkdir(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays']))
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_maxDfofAcrossDays.pdf']),'-dpdf', '-bestfit')

% trace avg max dfof over 4 days
avgdfof1 = mean(maxResp_D1(goodfit));
avgdfof2 = mean(maxResp_D2(goodfit));
avgdfof3 = mean(maxResp_D3(goodfit));
avgdfof4 = mean(maxResp_D4(goodfit));
figure;
x = [1 2 3 4];
avgdfof = [avgdfof1 avgdfof2 avgdfof3 avgdfof4];
cellnames = {'Day1'; 'Day2'; 'Day3'; 'Day4';};
bar(x,avgdfof)
set(gca,'xticklabel',cellnames)
ylabel('mean max dF/F')
title(mouse)
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_AvgMaxDfofAcrossDays.pdf']),'-dpdf', '-bestfit')

%% SNR
off_only1 = [ ];
for iTrial = 1:nTrials1
    off_only1 = [off_only1; squeeze(trial_tc1(nOff1/2:nOff1,iTrial,:))];
end
meanf1 = mean(off_only1,1);
stdf1 = std(off_only1,[],1);
SNR1 = meanf1./stdf1;
save(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], [ref_Date '_' mouse '_' run_str '_off_only.mat']), 'off_only1')

off_only2 = [ ];
for iTrial = 1:nTrials2
    off_only2 = [off_only2; squeeze(trial_tc2(nOff2/2:nOff2,iTrial,:))];
end
meanf2 = mean(off_only2,1);
stdf2 = std(off_only2,[],1);
SNR2 = meanf2./stdf2;
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_off_only.mat']), 'off_only2')

off_only3 = [ ];
for iTrial = 1:nTrials3
    off_only3 = [off_only3; squeeze(trial_tc3(nOff3/2:nOff3,iTrial,:))];
end
meanf3 = mean(off_only3,1);
stdf3 = std(off_only3,[],1);
SNR3 = meanf3./stdf3;
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_off_only.mat']), 'off_only3')

off_only4 = [ ];
for iTrial = 1:nTrials4
    off_only4 = [off_only4; squeeze(trial_tc4(nOff4/2:nOff4,iTrial,:))];
end
meanf4 = mean(off_only4,1);
stdf4 = std(off_only4,[],1);
SNR4 = meanf4./stdf4;
save(fullfile(fnout, [day4 '_' mouse], [day4 '_' mouse '_' run_str], [day4 '_' mouse '_' run_str '_off_only.mat']), 'off_only4')

figure;
subplot(1,3,1);
scatter(SNR1,SNR2);
axis square
xlim([0 max(SNR1)+1])
ylim([0 max(SNR2)+1])
refline(1,0)
xlabel('SNR Day 1')
ylabel('SNR Day 2')
R = corrcoef(SNR1,SNR2);
disp(R(1,2));
str = ['    r = ',num2str(R(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(1,3,2)
scatter(SNR1,SNR3);
axis square
xlim([0 max(SNR1)+1])
ylim([0 max(SNR3)+1])
refline(1,0)
xlabel('SNR Day 1')
ylabel('SNR Day 3')
R = corrcoef(SNR1,SNR3);
disp(R(1,2));
str = ['    r = ',num2str(R(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(1,3,3)
scatter(SNR1,SNR4);
axis square
xlim([0 max(SNR1)+1])
ylim([0 max(SNR4)+1])
refline(1,0)
xlabel('SNR Day 1')
ylabel('SNR Day 4')
R = corrcoef(SNR1,SNR4);
disp(R(1,2));
str = ['    r = ',num2str(R(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_SNR1-2.pdf']),'-dpdf', '-bestfit')

avgSNR1 = mean(SNR1);
avgSNR2 = mean(SNR2);
avgSNR3 = mean(SNR3);
avgSNR4 = mean(SNR4);
stdSNR1 = std(SNR1);
stdSNR2 = std(SNR2);
stdSNR3 = std(SNR3);
stdSNR4 = std(SNR4);
figure;
x = [1 2 3 4];
y = [avgSNR1 avgSNR2 avgSNR3 avgSNR4];
err = [stdSNR1/sqrt(nCells1) stdSNR2/sqrt(nCells1) stdSNR3/sqrt(nCells1) stdSNR4/sqrt(nCells1)];
bar(x,y);
hold on
er = errorbar(x,y,err);
er.Color = [0 0 0];
er.LineStyle = 'none';
cellnames = {'Day1'; 'Day2'; 'Day3'; 'Day4'};
set(gca, 'xticklabel', cellnames)
ylabel('average SNR')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_avgSNR.pdf']),'-dpdf', '-bestfit')

days = {ref_date; day2; day3; day4};
sz = size(off_only1);
% off_only = nan(sz(1),sz(2),length(days));
% off_only_mean = nan(sz(1),sz(2),length(days));
off_only_df = nan(sz(1),sz(2),length(days));
for iday = 1:length(days)
    iDay = cell2mat(days(iday));
    off_only = load(fullfile(fnout, [iDay '_' mouse], [iDay '_' mouse '_' run_str], [iDay '_' mouse '_' run_str '_off_only.mat']));
    off_only = struct2array(off_only);
    off_only_mean = mean(off_only,1);
    off_only_prctile = prctile(off_only,10,1);
    sz = size(off_only,1);
    off_only_df(1:sz,:,iday) = (off_only - off_only_prctile)./off_only_mean;
end

[n, n2] = subplotn(length(goodfit3));
start = 1;
figure;
for iC = 1:length(goodfit3)
    iCell = goodfit3(iC);
    subplot(n,n2,start)
    tcOffsetPlot(squeeze(off_only_df(:,iCell,:)));
    xlim([2500 3500])
    title(iCell)
    start = start + 1;
end
[hleg, hobj, hout, mout] = legend({'Day 1', 'Day 2', 'Day 3', 'Day 4'},'Position',[0.6 0.9 0.2 0.1],'Orientation','horizontal');
set(hobj,'linewidth',1.5);
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_TCoffset-mean3.pdf']),'-dpdf', '-bestfit')

off_only_df1 = off_only1 - prctile(off_only1,10,1);
snr1 = mean(off_only_df1,1)./std(off_only_df1,[],1);
off_only_df2 = off_only2 - prctile(off_only2,10,1);
snr2 = mean(off_only_df2,1)./std(off_only_df2,[],1);
off_only_df3 = off_only3 - prctile(off_only3,10,1);
snr3 = mean(off_only_df3,1)./std(off_only_df3,[],1);
off_only_df4 = off_only4 - prctile(off_only4,10,1);
snr4 = mean(off_only_df4,1)./std(off_only_df4,[],1);

[n, n2] = subplotn(length(goodfit3));
start = 1;
figure;
for iC = 1:length(goodfit3)
    iCell = goodfit3(iC);
    subplot(n,n2,start)
    x = [1 2 3 4];
    y = [snr1(:,iCell) snr2(:,iCell) snr3(:,iCell) snr4(:,iCell)];
%     err = [stdSNR1/sqrt(nCells1) stdSNR2/sqrt(nCells1) stdSNR3/sqrt(nCells1) stdSNR4/sqrt(nCells1)];
    bar(x,y);
    cellnames = {'Day1'; 'Day2'; 'Day3'; 'Day4'};
    set(gca, 'xticklabel', cellnames)
    ylabel('SNR')
    title(iCell)
    start = start + 1;
end
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_goodfitSNR-allDays2.pdf']),'-dpdf', '-bestfit')

%% baseline dfof within session
f_base1 = squeeze(mean(trial_tc1(base_wind1,:,:),1));
f_base_1 = mean(f_base1,2);
f_base2 = squeeze(mean(trial_tc2(base_wind2,:,:),1));
f_base_2 = mean(f_base2,2);
f_base3 = squeeze(mean(trial_tc3(base_wind3,:,:),1));
f_base_3 = mean(f_base3,2);
f_base4 = squeeze(mean(trial_tc4(base_wind4,:,:),1));
f_base_4 = mean(f_base4,2);
dfof_base_1 = mean(dfof_base1,2);
dfof_base_2 = mean(dfof_base2,2);
dfof_base_3 = mean(dfof_base3,2);
dfof_base_4 = mean(dfof_base4,2);
figure;
subplot(2,1,1)
xx = [1:nTrials1];
plot(xx,f_base_1')
hold on
plot(xx,f_base_2')
hold on
plot(xx,f_base_3')
hold on
plot(xx,f_base_4')
xlabel('nTrials')
ylabel('avg baseline f')
title(mouse)
legend({'Day 1', 'Day 2', 'Day 3', 'Day 4'},'Location','northwest','Orientation','horizontal')
subplot(2,1,2)
plot(xx,dfof_base_1')
hold on
plot(xx,dfof_base_2')
hold on
plot(xx,dfof_base_3')
hold on
plot(xx,dfof_base_4')
xlabel('nTrials')
ylabel('avg baseline dfof')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_Ftraces.pdf']),'-dpdf', '-bestfit')

%% Cell Maps
maskD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], [ref_date '_' mouse '_' run_str '_mask_cell.mat']));
mask_cell = maskD1.mask_cell;
data_dfof_max = maskD1.data_dfof_max;
reg_shifts = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], [ref_date '_' mouse '_' run_str '_reg_shifts.mat']));
reg = reg_shifts.data_reg_avg;
reg(find(reg>7000)) = 0;
reg = (reg./max(max(abs(reg))));
cell_list = intersect(1:nCells1, unique(mask_cell));
cell_stats = regionprops(mask_cell);

figure;
imagesc(reg); hold on;
bound = cell2mat(bwboundaries(mask_cell(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',0.5); hold on;
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if length(find(mask_cell == iCell))
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    else
        cell_list(iC) = NaN;
    end
    width = 30; height = 30;
    xCenter = cell_stats(iCell).Centroid(1);
    yCenter = cell_stats(iCell).Centroid(2);
    xLeft = xCenter - width/2;
    yBottom = yCenter - height/2;
    rectangle('Position', [xLeft, yBottom, width, height], 'EdgeColor', 'y', 'LineWidth', .8);
end
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_' ref_str '_boxedCellsReg.pdf']),'-dpdf', '-bestfit')

figure;
imagesc(data_dfof_max); hold on;
bound = cell2mat(bwboundaries(mask_cell(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',0.5); hold on;
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if length(find(mask_cell == iCell))
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    else
        cell_list(iC) = NaN;
    end
    width = 30; height = 30;
    xCenter = cell_stats(iCell).Centroid(1);
    yCenter = cell_stats(iCell).Centroid(2);
    xLeft = xCenter - width/2;
    yBottom = yCenter - height/2;
    rectangle('Position', [xLeft, yBottom, width, height], 'EdgeColor', 'y', 'LineWidth', .8);
end
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_' ref_str '_boxedCellsDfof.pdf']),'-dpdf', '-bestfit')

figure;
for iC = 1:length(cell_list)
    subplot(nCells1,1,iC)
    iCell = cell_list(iC);    
    width = 30; height = 30;
    xCenter = cell_stats(iCell).Centroid(2);
    yCenter = cell_stats(iCell).Centroid(1);
    xLeft = round(xCenter - width/2);
    yBottom = round(yCenter - height/2);
    imagesc(reg(xLeft:round(xLeft+width),yBottom:round(height+yBottom)))  
    axis square
    axis off
end
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_' ref_str '_structureMapReg.pdf']),'-dpdf', '-bestfit')

figure;
for iC = 1:length(cell_list)
    subplot(nCells1,1,iC)
    iCell = cell_list(iC);    
    width = 30; height = 30;
    xCenter = cell_stats(iCell).Centroid(2);
    yCenter = cell_stats(iCell).Centroid(1);
    xLeft = round(xCenter - width/2);
    yBottom = round(yCenter - height/2);
    imagesc(data_dfof_max(xLeft:round(xLeft+width),yBottom:round(height+yBottom)))  
    axis square
    axis off
end
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_' ref_str '_structureMapDfof.pdf']),'-dpdf', '-bestfit')

% figure;
% imagesc(reg);
% for iC = 1:length(cell_list)
%     iCell = cell_list(iC);
%     if length(find(mask_cell == iCell))
%         text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
%             'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
%         hold on
%     else
%         cell_list(iC) = NaN;
%     end
% end


%%
transformD2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_transform.mat']));
maskD2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_mask_cell.mat']));
mask_cell2 = maskD2.mask_cell;
reg2ref2 = transformD2.reg2ref;
ref = transformD2.ref;
data_dfof_max_ref = transformD2.data_dfof_max_ref;
reg2ref_dfof2 = transformD2.reg2ref_dfof;


transformD3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_transform.mat']));
maskD3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_mask_cell.mat']));
mask_cell3 = maskD3.mask_cell;
reg2ref3 = transformD3.reg2ref;
reg2ref_dfof3 = transformD3.reg2ref_dfof;

transformD4 = load(fullfile(fnout, [day4 '_' mouse], [day4 '_' mouse '_' run_str], [day4 '_' mouse '_' run_str '_transform.mat']));
maskD4 = load(fullfile(fnout, [day4 '_' mouse], [day4 '_' mouse '_' run_str], [day4 '_' mouse '_' run_str '_mask_cell.mat']));
mask_cell4 = maskD4.mask_cell;
reg2ref4 = transformD4.reg2ref;
reg2ref_dfof4 = transformD4.reg2ref_dfof;

cell_stats = regionprops(mask_cell2);
nCells = size(cell_stats,1);
% i1316
ROI = [75:275; 275:475]; 
% i1314
% ROI = [50:250; 300:500];

% see goodfit cells between day1 and day2 for all days
cell_list = intersect(goodfit, unique(mask_cell2(ROI(1,:),ROI(2,:))));
% figure;
% subplot(2,4,1)
% imagesc(ref);
% for iC = 1:length(cell_list)
%     iCell = cell_list(iC);
%     if length(find(mask_cell2 == iCell)) == length(find(mask_cell2(ROI(1,:),ROI(2,:)) == iCell))
%         text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
%             'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
%         hold on
%     else
%         cell_list(iC) = NaN;
%     end
% end
% xlim([ROI(2,1) ROI(2,end)])
% ylim([ROI(1,1) ROI(1,end)])
% axis square
% axis off
% title('Day 1')
% 
% subplot(2,4,2)
% imagesc(reg2ref2);
% for iC = 1:length(cell_list)
%     iCell = cell_list(iC);
%     if ~isnan(iCell)
%         text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
%             'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
%         hold on
%     end
% end
% xlim([ROI(2,1) ROI(2,end)])
% ylim([ROI(1,1) ROI(1,end)])
% axis square
% axis off
% title('Day 2')
% 
% subplot(2,4,3)
% imagesc(reg2ref3);
% for iC = 1:length(cell_list)
%     iCell = cell_list(iC);
%     if ~isnan(iCell)
%         text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
%             'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
%         hold on
%     end
% end
% xlim([ROI(2,1) ROI(2,end)])
% ylim([ROI(1,1) ROI(1,end)])
% axis square
% axis off
% title('Day 3')
% 
% subplot(2,4,4)
% imagesc(reg2ref4);
% for iC = 1:length(cell_list)
%     iCell = cell_list(iC);
%     if ~isnan(iCell)
%         text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
%             'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
%         hold on
%     end
% end
% xlim([ROI(2,1) ROI(2,end)])
% ylim([ROI(1,1) ROI(1,end)])
% axis square
% axis off
% title('Day 4')
% 
% subplot(2,4,5)
% imagesc(data_dfof_max_ref);
% for iC = 1:length(cell_list)
%     iCell = cell_list(iC);
%     if ~isnan(iCell)
%         text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
%             'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
%         hold on
%     end
% end
% xlim([ROI(2,1) ROI(2,end)])
% ylim([ROI(1,1) ROI(1,end)])
% axis square
% axis off
% 
% subplot(2,4,6)
% imagesc(reg2ref_dfof2);
% for iC = 1:length(cell_list)
%     iCell = cell_list(iC);
%     if ~isnan(iCell)
%         text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
%             'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
%         hold on
%     end
% end
% xlim([ROI(2,1) ROI(2,end)])
% ylim([ROI(1,1) ROI(1,end)])
% axis square
% axis off
% 
% subplot(2,4,7)
% imagesc(reg2ref_dfof3);
% for iC = 1:length(cell_list)
%     iCell = cell_list(iC);
%     if ~isnan(iCell)
%         text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
%             'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
%         hold on
%     end
% end
% xlim([ROI(2,1) ROI(2,end)])
% ylim([ROI(1,1) ROI(1,end)])
% axis square
% axis off
% 
% subplot(2,4,8)
% imagesc(reg2ref_dfof4);
% for iC = 1:length(cell_list)
%     iCell = cell_list(iC);
%     if ~isnan(iCell)
%         text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
%             'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
%         hold on
%     end
% end
% xlim([ROI(2,1) ROI(2,end)])
% ylim([ROI(1,1) ROI(1,end)])
% axis square
% axis off
% % suptitle([mouse ' days' ref_date ' thru' day4], 'color', 'red')
% % print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref '_' mouse '_' ref_str '_cellAlignZoom.pdf']),'-dpdf', '-bestfit')

nC = sum(~isnan(cell_list));
figure; 
start = 1;
for iC = 1:length(cell_list)
    if ~isnan(cell_list(iC))
        iCell = cell_list(iC);
        subplot(6, 2, start)
        plot(0:180, oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D1.avgResponseEaOri(iCell,:), oriTuning_D1.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D1(iCell)) ' deg'])
        subplot(6, 2, start+1)
        plot(0:180, oriTuning_D2.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D2.avgResponseEaOri(iCell,:), oriTuning_D2.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D2(iCell)) ' deg'])
%         subplot(5, 4, start+2)
%         plot(0:180, oriTuning_D3.vonMisesFitAllCellsAllBoots(:,1,iCell))
%         hold on
%         errorbar(0:22.5:180-22.5, oriTuning_D3.avgResponseEaOri(iCell,:), oriTuning_D3.semResponseEaOri(iCell,:))
%         title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D3(iCell)) ' deg'])
%         subplot(5, 4, start+3)
%         plot(0:180, oriTuning_D4.vonMisesFitAllCellsAllBoots(:,1,iCell))
%         hold on
%         errorbar(0:22.5:180-22.5, oriTuning_D4.avgResponseEaOri(iCell,:), oriTuning_D4.semResponseEaOri(iCell,:))
%         title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D4(iCell)) ' deg'])
        start = start+2;
    end
end
suptitle([mouse ' day 1 to day 2 tuning'])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_alignTunin1-2.pdf']),'-dpdf', '-bestfit')

cell_list = intersect(goodfit3to1, unique(mask_cell2(ROI(1,:),ROI(2,:))));
nC = sum(~isnan(cell_list));
figure; 
start = 1;
for iC = 1:length(cell_list)
    if ~isnan(cell_list(iC))
        iCell = cell_list(iC);
        subplot(6, 2, start)
        plot(0:180, oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D1.avgResponseEaOri(iCell,:), oriTuning_D1.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D1(iCell)) ' deg'])
        subplot(6, 2, start+1)
        plot(0:180, oriTuning_D3.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D3.avgResponseEaOri(iCell,:), oriTuning_D3.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D2(iCell)) ' deg'])
        start = start+2;
    end
end
suptitle([mouse ' day 1 to day 3 tuning'])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_alignTunin1-3.pdf']),'-dpdf', '-bestfit')

cell_list = intersect(goodfit4to1, unique(mask_cell2(ROI(1,:),ROI(2,:))));
nC = sum(~isnan(cell_list));
figure; 
start = 1;
for iC = 1:length(cell_list)
    if ~isnan(cell_list(iC))
        iCell = cell_list(iC);
        subplot(6, 2, start)
        plot(0:180, oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D1.avgResponseEaOri(iCell,:), oriTuning_D1.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D1(iCell)) ' deg'])
        subplot(6, 2, start+1)
        plot(0:180, oriTuning_D4.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D4.avgResponseEaOri(iCell,:), oriTuning_D4.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D2(iCell)) ' deg'])
        start = start+2;
    end
end
suptitle([mouse ' day 1 to day 4 tuning'])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_alignTunin1-4.pdf']),'-dpdf', '-bestfit')

%% show cumulative goodfit cells across days

cell_list = intersect(goodfit_D1, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,1)
imagesc(data_dfof_max_ref);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 1')

cell_list = intersect(goodfit, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,2)
imagesc(reg2ref_dfof2);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 2')

cell_list = intersect(goodfit2, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,3)
imagesc(reg2ref_dfof3);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 3')

cell_list = intersect(goodfit3, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,4)
imagesc(reg2ref_dfof4);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 4')

suptitle([mouse ' cumulative goodfit cells across days'])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_goodfitAcrossDays.pdf']),'-dpdf', '-bestfit')

%% how each day individually aligns to day 1
cell_list = intersect(goodfit_D1, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,1)
imagesc(data_dfof_max_ref);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 1')

cell_list = intersect(goodfit, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,2)
imagesc(reg2ref_dfof2);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 2')

cell_list = intersect(goodfit3to1, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,3)
imagesc(reg2ref_dfof3);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 3')

cell_list = intersect(goodfit4to1, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,4)
imagesc(reg2ref_dfof4);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 4')

suptitle([mouse ' days 2-4 goodfit with day 1'])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_goodfitToDay1.pdf']),'-dpdf', '-bestfit')

%% each day's individual goodfit cells 
cell_list = intersect(goodfit_D1, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,1)
imagesc(data_dfof_max_ref);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 1')

cell_list = intersect(goodfit_D2, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,2)
imagesc(reg2ref_dfof2);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 2')

cell_list = intersect(goodfit_D3, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,3)
imagesc(reg2ref_dfof3);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 3')

cell_list = intersect(goodfit_D4, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,4)
imagesc(reg2ref_dfof4);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 4')

% suptitle([mouse ' goodfit cells on each day'])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_goodfitSepDays.pdf']),'-dpdf', '-bestfit')

%% tuning across all days
cell_list = intersect(goodfit3, unique(mask_cell2(ROI(1,:),ROI(2,:))));
nC = sum(~isnan(cell_list));
figure; 
start = 1;
for iC = 1:length(cell_list)
    if ~isnan(cell_list(iC))
        iCell = cell_list(iC);
        subplot(5, 4, start)
        plot(0:180, oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D1.avgResponseEaOri(iCell,:), oriTuning_D1.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D1(iCell)) ' deg'])
        subplot(5, 4, start+1)
        plot(0:180, oriTuning_D2.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D2.avgResponseEaOri(iCell,:), oriTuning_D2.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D2(iCell)) ' deg'])
        subplot(5, 4, start+2)
        plot(0:180, oriTuning_D3.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D3.avgResponseEaOri(iCell,:), oriTuning_D3.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D3(iCell)) ' deg'])
        subplot(5, 4, start+3)
        plot(0:180, oriTuning_D4.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D4.avgResponseEaOri(iCell,:), oriTuning_D4.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D4(iCell)) ' deg'])
        start = start+4;
    end
end
suptitle([mouse ' tuning across all days'])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_alignTuning.pdf']),'-dpdf', '-bestfit')


