%% running and pupil size-tuning experiment analysis
%
%% load csv to define exp
% updated to include ret_str as a header, only use for runningPupil list for now
clear all; clc;
expfile = 'C:\Users\kevin\Documents\Repositories\ImagingCode-Glickfeld-Hull\kevin\Size Tuning Code\runningPupil_experiment_list.txt';
fID = fopen(expfile);
head = textscan(fID,'%s%s%s%s%s%s',1,'delimiter',',');
head = vertcat(head{:});
temp = textscan(fID,'%s%s%s%s%s%s','delimiter',',','HeaderLines',1);
temp = horzcat(temp{:});
expdata = cell2table(temp,'VariableNames',head);
nExp = size(expdata,1);
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);

fprintf(['Size-tuning Running + Pupil analysis - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

%% load each experiment and concatenate data

fprintf('\nBegin loading and concatentating experiment data...\n')
% sizeTuneData
fprintf('Loading sizeTuneData (raw size-tuning data)\n')
sizeTune_all = cell(0);
sizeMean_all = [];
sizeSEM_all = [];
cellDists_all = [];
tc_dfof_all = [];
Ind_struct_all = [];
nCellsExp = zeros(1,nExp);
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = expdata.run_str{i};
    %filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_sizeTuneData.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'sizeTune', 'sizeMean', 'sizeSEM', 'cellDists')
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_Tuning.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'tc_dfof','Ind_struct')
    if size(sizeTune,2) == 6 % if 6 con, only take middle 4 cons (2-5)
        sizeTune_all = cat(3,sizeTune_all,sizeTune(:,2:5,:));
        sizeMean_all = cat(3,sizeMean_all,sizeMean(:,2:5,:));
        sizeSEM_all = cat(3,sizeSEM_all,sizeSEM(:,2:5,:));
    else
        sizeTune_all = cat(3,sizeTune_all,sizeTune);
        sizeMean_all = cat(3,sizeMean_all,sizeMean);
        sizeSEM_all = cat(3,sizeSEM_all,sizeSEM);
        tc_dfof_all = cat(2,tc_dfof_all,tc_dfof);
        Ind_struct_all = cat(1, Ind_struct_all, Ind_struct);
    end
    cellDists_all = [cellDists_all;cellDists];
    nCellsExp(i) = length(cellDists);
    fprintf('done\n')
end

% sizeFitResults_SP
fprintf('Loading sizeFitResults_SP (true fit at all cons)\n')
sizeFits_all = struct([]); % no cells, 4 cons
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = expdata.run_str{i};
    %filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_SP.mat']);
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_SP.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_sizeFitResults_SP.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'sizeFits')
    %sizeFits_all = cat(1,sizeFits_all,sizeFits);
    if size(sizeFits,2) == 6 % if 6 con, only take middle 4 cons (2-5)
        sizeFits_all = cat(1,sizeFits_all,sizeFits(:,2:5));
    else
        sizeFits_all = cat(1,sizeFits_all,sizeFits);
    end
    fprintf('done\n')
end

% lbub_fits
fprintf('Loading good fit inds (ret)\n')
goodfit_ind_all = [];
nCellsExpRet = zeros(1,nExp);
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = expdata.ret_str{i};
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind')
    nCellsExpRet(i) = size(lbub_fits,1);
    tempinds = sum(nCellsExpRet(1:i-1)) + goodfit_ind; % offset by # cells in previous exps
    goodfit_ind_all = [goodfit_ind_all tempinds];
    
    fprintf('done\n')
end

% lbub_fits
fprintf('Loading lbub_fits\n')
lbub_fits_all = [];
goodfit_ind_size_all = [];
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = expdata.run_str{i};
    %filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind_size')
    lbub_fits_all = cat(1,lbub_fits_all,lbub_fits);
    tempinds = sum(nCellsExp(1:i-1)) + goodfit_ind_size; % offset by # cells in previous exps
    goodfit_ind_size_all = [goodfit_ind_size_all tempinds];
    fprintf('done\n')
end

% eye&running trials
fprintf('Loading pupil and running data\n')

filename = fullfile('N:\home\ashley\Manuscripts\Size Tuning\Matlab Analysis', 'eye&RunStruct.mat');
if ~exist(filename, 'file')
    fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
end
load(filename)
for iExp = 1:nExp
    eyeRun(iExp).relEyePosition = eyeRun(iExp).relEyePostion;
end
eyeRun = rmfield(eyeRun,'relEyePostion');
fprintf('done\n')

fprintf(['\nFinished loading all ' num2str(nExp) ' experiments.\n'])

nCellsTot = sum(nCellsExp);
fprintf('%d cells loaded (%d goodfit_size)\n',nCellsTot,length(goodfit_ind_size_all))

expInd = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExp(i))];
end

cons = [0.1 0.2 0.4 0.8]; nCon = length(cons);
szs = 5*1.5.^[0:7]; nSize = length(szs);
szRng = linspace(0,max(szs));

%% examine data
% may need to re-do size fits bootstrapping?
% for now just use bootstrapping data (good fits) from all trials
% and build size tuning curves from running/stationary trials in those cells

% find running/stationary/pupil-in-range trials and build new sizeTune data
run_cutoff = 2; %in cm/s
pup_cutoff = 5; %in deg
for iExp = 1:nExp
    %run_trials = find(eyeRun(iExp).speedCPS > run_cutoff);
    stat_trials = find(eyeRun(iExp).speedCPS <= run_cutoff);
    pup_trials = find(eyeRun(iExp).relEyePosition <= pup_cutoff);
    for iSz = 1:nSize
        Ind_struct_all(iExp,iSz).stat_trials = intersect(Ind_struct_all(iExp,iSz).all_trials,stat_trials);
        Ind_struct_all(iExp,iSz).pup_trials = intersect(Ind_struct_all(iExp,iSz).all_trials,pup_trials);
    end
end
sizeTune_all_stat = cell(size(sizeTune_all));
sizeTune_all_pup = cell(size(sizeTune_all));
for iCell = 1:nCellsTot
    iExp = expInd(iCell);
    for iSz = 1:nSize
        stat_inds = find(ismember(Ind_struct_all(iExp,iSz).all_trials,Ind_struct_all(iExp,iSz).stat_trials));
        pup_inds = find(ismember(Ind_struct_all(iExp,iSz).all_trials,Ind_struct_all(iExp,iSz).pup_trials));
        sizeTune_all_stat{iSz,1,iCell} = sizeTune_all{iSz,1,iCell}(stat_inds);
        sizeTune_all_pup{iSz,1,iCell} = sizeTune_all{iSz,1,iCell}(pup_inds);
    end
end
sizeMean_all_stat = cellfun(@mean,sizeTune_all_stat);
sizeSEM_all_stat = cellfun(@(x) std(x)/sqrt(length(x)),sizeTune_all_stat);
sizeMean_all_pup = cellfun(@mean,sizeTune_all_pup);
sizeSEM_all_pup = cellfun(@(x) std(x)/sqrt(length(x)),sizeTune_all_pup);

%% re-run size tuning fits
szRng = linspace(0,max(szs));
logfit1 = @(coefs,xdata) coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3))))
logfit2 = @(coefs,xdata) coefs(1)./(1+exp(-(coefs(2)+coefs(5))*(xdata-coefs(3)))) - coefs(4)./(1+exp(-coefs(5)*(xdata-(coefs(3)+coefs(6)))))
% store # trials at each size for this con con (same for all cells)
nTr = zeros(nSize,1);
for iSz = 1:nSize
    nTr(iSz) = length(sizeTune_all{iSz,1,1});
end
cd('C:\Users\kevin\Documents\Repositories\ImagingCode-Glickfeld-Hull\kevin\Size Tuning Code')
opts = optimset('Display','off');

% initialize sizeFits struct with last cell, last con
nPts = floor(mean(nTr));
dumdum = zeros(1,nPts);%[0]; % define zero point for dF/F
szs0 = zeros(1,nPts);%[0.1]; % and size
for iSz = 1:nSize
    nPts = nPts + nTr(iSz);
    dumdum = [dumdum sizeTune_all{iSz,1,nCellsTot}'];
    szs0 = [szs0 szs(iSz)*ones(1,nTr(iSz))];
end
% max of each size mean for use in initial guesses
maxMean = max(sizeMean_all(:,1,nCellsTot));
x_max = szs(find(sizeMean_all(:,1,nCellsTot)==maxMean,1));
if isempty(x_max)
    x_max = 15;
end

PLOTIT_FIT = 0;
SAVEALLDATA = 1;
Fit_SizeTuneSmooth_KM % call fit script, returns fit structure s, no plot
%eval(['sizeFits(iCell,iCon)',' = s;']);

f = fieldnames(s)';
f{2,1} = {};
sizeFits_stat=struct(f{:});
sizeFits_stat(nCellsTot,1) = s;
sizeFits_pup=sizeFits_stat;

fprintf('Begin fitting size-tuning curves at all cells, all contrasts...')
for iCell = 1:nCellsTot
    fprintf(['\nCell# ' num2str(iCell) '/' num2str(nCellsTot) ', (RF-stim dist: ' num2str(cellDists_all(iCell)) ' deg)'])
    
    fprintf('.')
    % store # trials at each size for this con con (same for all cells)
    nTr_stat = zeros(nSize,1); nTr_pup=nTr_stat;
    for iSz = 1:nSize
        nTr_stat(iSz) = length(sizeTune_all_stat{iSz,1,iCell});
        nTr_pup(iSz) = length(sizeTune_all_pup{iSz,1,iCell});
    end
    
    % first use stationary trials
    nPts = floor(mean(nTr_stat));
    dumdum = zeros(1,nPts);%[0]; % define zero point for dF/F
    szs0 = zeros(1,nPts);%[0.1]; % and size
    for iSz = 1:nSize
        nPts = nPts + nTr_stat(iSz,1);
        dumdum = [dumdum sizeTune_all_stat{iSz,1,iCell}'];
        szs0 = [szs0 szs(iSz)*ones(1,nTr_stat(iSz,1))];
    end
    % max of each size mean for use in initial guesses
    maxMean = max(sizeMean_all_stat(:,1,iCell));
    x_max = szs(find(sizeMean_all_stat(:,1,iCell)==maxMean,1));
    
    PLOTIT_FIT = 0;
    SAVEALLDATA = 1;
    Fit_SizeTuneSmooth_KM % call fit script, returns fit structure s, no plot
    sizeFits_stat(iCell,1)=s;
    
    % again for pupil trials
    nPts = floor(mean(nTr_pup));
    dumdum = zeros(1,nPts);%[0]; % define zero point for dF/F
    szs0 = zeros(1,nPts);%[0.1]; % and size
    for iSz = 1:nSize
        nPts = nPts + nTr_pup(iSz,1);
        dumdum = [dumdum sizeTune_all_pup{iSz,1,iCell}'];
        szs0 = [szs0 szs(iSz)*ones(1,nTr_pup(iSz,1))];
    end
    % max of each size mean for use in initial guesses
    maxMean = max(sizeMean_all_pup(:,1,iCell));
    x_max = szs(find(sizeMean_all_pup(:,1,iCell)==maxMean,1));
    
    PLOTIT_FIT = 0;
    SAVEALLDATA = 1;
    Fit_SizeTuneSmooth_KM % call fit script, returns fit structure s, no plot
    sizeFits_pup(iCell,1)=s;
end

% Save fit results as sizeFitResults.mat
% save sizeFitResults.mat, with sizeFits struct
%save(filename, 'sizeFits')
fprintf('\nDid not save sizeFitResults_SP.mat\n')


%% plot
% plots:
% Supp 3 - running (set a cutoff for running speed?)
% mean size tuning curve in PM (with distance cutoff at 20), stationary vs all?
% pref size: violin plots for stationary vs all
% ^ for suppind
% Supp 4 - pupil (set a cutoff for pupil position?)
% Comparison of RF size, preferred size and SI on all trials
% vs those curated for being close to the median.

RFcutoff = 20; %20 deg for PM
ind = intersect(goodfit_ind_size_all,find(cellDists_all < RFcutoff));
figure(3);clf 
% first for running
subplot(2,3,1) % mean curves
hold on
%plot(szs,squeeze(sizeMean_all(:,:,ind)),'-','Color',[0 0 0 0.1],'HandleVisibility','off');
errorbar(szs,mean(sizeMean_all(:,:,ind),3),geo_mean(sizeSEM_all(:,:,ind),3))
errorbar(szs,mean(sizeMean_all_stat(:,:,ind),3),geo_mean(sizeSEM_all_stat(:,:,ind),3))
title('Running Analysis')
xlabel('Size (deg)')
ylabel('dF/F')
ylim([0 0.15])
legend('All trials','Stationary','Location','best')
ax = subplot(2,3,2); % SI
% %5 sideways violin
% data_cell{1} = [sizeFits_all(ind).suppInd];
% data_cell{2} = [sizeFits_stat(ind).suppInd];
% y_mean = [mean(data_cell{1}) mean(data_cell{2})];
% y_std = [std(data_cell{1}) std(data_cell{2})];
% for i = 1:2
%     [fi xi] = ksdensity(data_cell{i},'BoundaryCorrection','reflection','Support',[0-eps 1.2+eps]); %,'BoundaryCorrection','reflection'
%     fnorm(:,i) = fi/max(fi)*0.3;
%     xinorm(:,i) = xi;
% end
% colors = get(ax,'ColorOrder');
% for i=1:2
%     hold on
%     h5(i)=fill([xinorm(:,i);flipud(xinorm(:,i))],[fnorm(:,i)+(3-i);flipud((3-i)-fnorm(:,i))],[1 1 1],'EdgeColor','k');
%     p(1)=plot([y_mean(i) y_mean(i)],[interp1(xinorm(:,i),fnorm(:,i)+(3-i),y_mean(i)), interp1(flipud(xinorm(:,i)),flipud((3-i)-fnorm(:,i)),y_mean(i)) ],'k','LineWidth',2);
%     h5(i).FaceColor = colors(i,:);
% end
% axis([0 1 0.5 2.5]);
% legend off
% ax.YTick = [1:2];
% ax.YTickLabel = fliplr(["all","stat"]);
% set(gca,'box','off','TickDir','out')
% ylabel('Trials')
%xlabel('SI')
%histogram([sizeFits_all(ind).suppInd])
%hold on 
%histogram([sizeFits_stat(ind).suppInd])
%legend('All','Stat')
plot([sizeFits_all(ind).suppInd],[sizeFits_stat(ind).suppInd],'.');
hold on
plot([0 1.5],[0 1.5],'-','Color',[0 0 0 0.1])
xlabel('All trials')
ylabel('Stationary trials')
title('SI')
ax = subplot(2,3,3); % prefSize
% %5 sideways violin
% data_cell{1} = [sizeFits_all(ind).prefSize];
% data_cell{2} = [sizeFits_stat(ind).prefSize];
% y_mean = [mean(data_cell{1}) mean(data_cell{2})];
% y_std = [std(data_cell{1}) std(data_cell{2})];
% for i = 1:2
%     [fi xi] = ksdensity(data_cell{i},'BoundaryCorrection','reflection','Support',[0-eps 85+eps]); %,'BoundaryCorrection','reflection'
%     fnorm(:,i) = fi/max(fi)*0.3;
%     xinorm(:,i) = xi;
% end
% colors = get(ax,'ColorOrder');
% for i=1:2
%     hold on
%     h5(i)=fill([xinorm(:,i);flipud(xinorm(:,i))],[fnorm(:,i)+(3-i);flipud((3-i)-fnorm(:,i))],[1 1 1],'EdgeColor','k');
%     p(1)=plot([y_mean(i) y_mean(i)],[interp1(xinorm(:,i),fnorm(:,i)+(3-i),y_mean(i)), interp1(flipud(xinorm(:,i)),flipud((3-i)-fnorm(:,i)),y_mean(i)) ],'k','LineWidth',2);
%     h5(i).FaceColor = colors(i,:);
% end
% axis([0 85 0.5 2.5]);
% legend off
% ax.YTick = [1:2];
% ax.YTickLabel = fliplr(["all","stat"]);
% set(gca,'box','off','TickDir','out')
% ylabel('Trials')
%xlabel('prefSize')
%histogram([sizeFits_all(ind).prefSize])
%hold on 
%histogram([sizeFits_stat(ind).prefSize])
%legend('All','Stat')
plot([sizeFits_all(ind).prefSize],[sizeFits_stat(ind).prefSize],'.');
hold on
plot([0 100],[0 100],'-','Color',[0 0 0 0.1])
xlabel('All trials')
ylabel('Stationary trials')
title('Preferred Size')

% then for pupil
subplot(2,3,4) % mean curves
hold on
%plot(szs,squeeze(sizeMean_all(:,:,ind)),'-','Color',[0 0 0 0.1],'HandleVisibility','off');
errorbar(szs,mean(sizeMean_all(:,:,ind),3),geo_mean(sizeSEM_all(:,:,ind),3))
errorbar(szs,mean(sizeMean_all_pup(:,:,ind),3),geo_mean(sizeSEM_all_pup(:,:,ind),3))
title('Pupil Analysis')
xlabel('Size (deg)')
ylabel('dF/F')
ylim([0 0.15])
legend('All trials','Pupil on (5 deg)','Location','best')
subplot(2,3,5); % SI
plot([sizeFits_all(ind).suppInd],[sizeFits_pup(ind).suppInd],'.');
hold on
plot([0 1.5],[0 1.5],'-','Color',[0 0 0 0.1])
xlabel('All trials')
ylabel('Pupil-on trials')
title('SI')
subplot(2,3,6); % prefSize
plot([sizeFits_all(ind).prefSize],[sizeFits_pup(ind).prefSize],'.');
hold on
plot([0 100],[0 100],'-','Color',[0 0 0 0.1])
xlabel('All trials')
ylabel('Pupil-on trials')
title('Preferred Size')