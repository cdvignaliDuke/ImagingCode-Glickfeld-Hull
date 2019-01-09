%% compare size-tuning across visual areas
% script that loads given experiments and sorts by visual area to compare
% size tuning parameters across V1 + HVAs

% given csv .txt file for experiment list
% includes: date, mouse, run_str(size-tuning exp), HVA, indicator

% data to load:
% ?ret results - lbub_fits (w/ RF location), goodfit cells (maybe?)
% ?size tuning run input - stim location (maybe?)
% _sizeTuneData.mat - sizeTune, sizeMean, sizeSEM, cellDists
% _sizeFitResults_SP.mat - sizeFits struct (all cons)
% _Fit_struct.mat - Fit_struct (highest con, with shuffles)
% _lbub_fits - lbub_fits, goodfit_ind_size

%% load csv to define exp

clear all; clc;
expfile = '\\CRASH.dhe.duke.edu\data\home\kevin\Code\Ai9x_experiment_list.txt';
fID = fopen(expfile);
head = textscan(fID,'%s%s%s%s%s',1,'delimiter',',');
head = vertcat(head{:});
temp = textscan(fID,'%s%s%s%s%s','delimiter',',','HeaderLines',1);
temp = horzcat(temp{:});
expdata = cell2table(temp,'VariableNames',head);
nExp = size(expdata,1);
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);

fprintf(['Size-tuning visual-area comparison analysis - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

%% load each experiment and concatenate data

fprintf('\nBegin loading and concatentating experiment data...\n')
% sizeTuneData
fprintf('Loading sizeTuneData (raw size-tuning data)\n')
sizeTune_all = cell(0);
sizeMean_all = [];
sizeSEM_all = [];
cellDists_all = [];
nCellsExp = zeros(1,nExp);
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = expdata.run_str{i};
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_sizeTuneData.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'sizeTune', 'sizeMean', 'sizeSEM', 'cellDists')
    sizeTune_all = cat(3,sizeTune_all,sizeTune);
    sizeMean_all = cat(3,sizeMean_all,sizeMean);
    sizeSEM_all = cat(3,sizeSEM_all,sizeSEM);
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
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_SP.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_sizeFitResults_SP.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'sizeFits')
    sizeFits_all = cat(1,sizeFits_all,sizeFits);
    fprintf('done\n')
end

% % Fit_struct - need this?
% fprintf('Loading Fit_struct (highest con, with shuffling)\n')
% Fit_struct_all = struct;
% for i=1:nExp
%     fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
%     date = expdata.date{i};
%     mouse = expdata.mouse{i};
%     run_str = expdata.run_str{i};
%     filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Fit_struct.mat']);
%     if ~exist(filename, 'file')
%         fprintf([[date '_' mouse '_' run_str '_Fit_struct.mat'] ' not found! Please remove from list\n'])
%     end
%     load(filename, 'Fit_struct')
%     Fit_struct_all = cat(1,Fit_struct_all,Fit_struct);
%     fprintf('done\n')
% end

% lbub_fits
fprintf('Loading lbub_fits\n')
lbub_fits_all = [];
goodfit_ind_size_all = [];
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = expdata.run_str{i};
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind_size')
    lbub_fits_all = cat(1,lbub_fits_all,lbub_fits);
    tempinds = sum(nCellsExp(1:i-1)) + goodfit_ind_size; % offset by # cells in previous exps
    goodfit_ind_size_all = [goodfit_ind_size_all tempinds];
    fprintf('done\n')
end
fprintf(['\nFinished loading all ' num2str(nExp) ' experiments.\n'])

%% continue setting up data
% define model1, model2 cells (0,1,2)
% then extract each visual area, using RF-dist <10 cutoff

expInd = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExp(i))];
end

% is model1 + is model2
ism1 = find(~lbub_fits_all(goodfit_ind_size_all,8,4));
ism2 = find(lbub_fits_all(goodfit_ind_size_all,8,4));

% extract each visual area (V1, LM, AL, PM)
% take raw data (sizeTune, sizeMean, sizeSEM), sizeFits, lbub_fits, ism1/2?
% find indices with: goodfit_ind_size, expdata.area, RFcutoff
RFcutoff_V1 = find(cellDists_all<10);
RFcutoff_LMAL = find(cellDists_all<15);
RFcutoff_PM = find(cellDists_all<20);

% area experiment indices
V1exp = find(cellfun(@(x) strcmp(x,'V1'), expdata.area, 'UniformOutput', 1));
LMexp = find(cellfun(@(x) strcmp(x,'LM'), expdata.area, 'UniformOutput', 1));
ALexp = find(cellfun(@(x) strcmp(x,'AL'), expdata.area, 'UniformOutput', 1));
PMexp = find(cellfun(@(x) strcmp(x,'PM'), expdata.area, 'UniformOutput', 1));
% cell indices by area (goodfit_size+RF<10)
V1ind = intersect(intersect(find(ismember(expInd,V1exp)),goodfit_ind_size_all),RFcutoff_V1);
LMind = intersect(intersect(find(ismember(expInd,LMexp)),goodfit_ind_size_all),RFcutoff_LMAL);
ALind = intersect(intersect(find(ismember(expInd,ALexp)),goodfit_ind_size_all),RFcutoff_LMAL);
PMind = intersect(intersect(find(ismember(expInd,PMexp)),goodfit_ind_size_all),RFcutoff_PM);

% now extract sizeTune,sizeMean,sizeSEM,sizeFits,lbub_fits : from each area
%v1
sizeTune_V1 = sizeTune_all{:,:,V1ind}; % (size,con,cell)
sizeMean_V1 = sizeMean_all(:,:,V1ind);
sizeSEM_V1 = sizeSEM_all(:,:,V1ind);
sizeFits_V1 = sizeFits_all(V1ind,:); %cell,con
lbub_fits_V1 = lbub_fits_all(V1ind,:,:); %cell,par,val (low up mean true stdev)
ism1_V1 = reshape(~[sizeFits_V1.Ftest],size(sizeFits_V1));
ism2_V1 = reshape([sizeFits_V1.Ftest],size(sizeFits_V1));
%ism1_V1 = find(~lbub_fits_V1(:,8,4));
%ism2_V1 = find(lbub_fits_V1(:,8,4));
%lm
sizeTune_LM = sizeTune_all{:,:,LMind}; % (size,con,cell)
sizeMean_LM = sizeMean_all(:,:,LMind);
sizeSEM_LM = sizeSEM_all(:,:,LMind);
sizeFits_LM = sizeFits_all(LMind,:); %cell,con
lbub_fits_LM = lbub_fits_all(LMind,:,:); %cell,par,val (low up mean true stdev)
ism1_LM = reshape(~[sizeFits_LM.Ftest],size(sizeFits_LM));
ism2_LM = reshape([sizeFits_LM.Ftest],size(sizeFits_LM));
%al
sizeTune_AL = sizeTune_all{:,:,ALind}; % (size,con,cell)
sizeMean_AL = sizeMean_all(:,:,ALind);
sizeSEM_AL = sizeSEM_all(:,:,ALind);
sizeFits_AL = sizeFits_all(ALind,:); %cell,con
lbub_fits_AL = lbub_fits_all(ALind,:,:); %cell,par,val (low up mean true stdev)
ism1_AL = reshape(~[sizeFits_AL.Ftest],size(sizeFits_AL));
ism2_AL = reshape([sizeFits_AL.Ftest],size(sizeFits_AL));
%pm
sizeTune_PM = sizeTune_all{:,:,PMind}; % (size,con,cell)
sizeMean_PM = sizeMean_all(:,:,PMind);
sizeSEM_PM = sizeSEM_all(:,:,PMind);
sizeFits_PM = sizeFits_all(PMind,:); %cell,con
lbub_fits_PM = lbub_fits_all(PMind,:,:); %cell,par,val (low up mean true stdev)
ism1_PM = reshape(~[sizeFits_PM.Ftest],size(sizeFits_PM));
ism2_PM = reshape([sizeFits_PM.Ftest],size(sizeFits_PM));

nExp_area = [length(V1exp) length(LMexp) length(ALexp) length(PMexp)];
nCells_area = [length(V1ind) length(LMind) length(ALind) length(PMind)];

%% present each exp to examine cells+fits and choose example cells
% show each cell rawdata w/ model1+2 overlaid

% (technically should load this from experiment file)
cons = [0.1 0.2 0.4 0.8]; nCon = length(cons);
szs = 5*1.5.^[0:7]; nSize = length(szs);
szRng = linspace(0,max(szs));

for iExp = 1:nExp
    % select only cells in exp i with goodfits and RF<10
    
    switch expdata.area{iExp}
        case 'V1'
            RFcutoff = find(cellDists_all<10);
        case {'LM', 'AL'}
            RFcutoff = find(cellDists_all<15);
        case 'PM'
            RFcutoff = find(cellDists_all<20);
    end
    
    inds = intersect(intersect(find(ismember(expInd,iExp)),goodfit_ind_size_all),RFcutoff);
    n = length(inds);
    
    figure;
    start = 1;
    ifig = 1;
    for i=1:n
        iCell = inds(i);
        s = sizeFits_all(iCell,nCon);
        if start ==37
            suptitle(['Exp: ' num2str(iExp)])
            set(gcf, 'Position', [0 0 800 1000]);
            fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P\HVAcomparison', ['exp' num2str(iExp) '_SizeTuneFits' num2str(ifig) '.pdf']);
            print(fn_out,'-dpdf')
            figure;
            ifig = 1+ifig;
            start = 1;
        end
        h = subplot(6,6,start);
        errorbar([0 szs],[0 sizeMean_all(:,nCon,iCell)'],[0 sizeSEM_all(:,nCon,iCell)'])
        hold on
        plot(s.szs0,s.data,'.k')
        plot(szRng,s.fitout1,'-b')
        plot(szRng,s.fitout2,'-r')
        hold off
        ylim([min([-0.5*s.maxResp1 min(s.data)]) 1.2*max([s.maxResp2 max(s.data)])])
        if s.Fstr1
            title(['#' num2str(iCell), ' m1 ' num2str(cellDists_all(iCell),3)]);
        else
            title(['#' num2str(iCell), ' m2 ' num2str(cellDists_all(iCell),3)]);
        end
        
        start = start+1;
    end
    suptitle(['Exp: ' num2str(iExp)])
    set(gcf, 'Position', [0 0 800 1000]);
    fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P\HVAcomparison', ['exp_' num2str(iExp) '_SizeTuneFits' num2str(ifig) '.pdf']);
    print(fn_out,'-dpdf')
end

%% now look at each area
areas_c = categorical({'V1','LM','AL','PM'},{'V1','LM','AL','PM'});

% first examine cell counts, proportion of each model
figure(1);clf;
modelcounts = [sum(~[sizeFits_V1.Ftest]) sum([sizeFits_V1.Ftest]); ...
    sum(~[sizeFits_LM.Ftest]) sum([sizeFits_LM.Ftest]); ...
    sum(~[sizeFits_AL.Ftest]) sum([sizeFits_AL.Ftest]); ...
    sum(~[sizeFits_PM.Ftest]) sum([sizeFits_PM.Ftest])];
modelcounts_norm = modelcounts./sum(modelcounts,2);
bar(areas_c,modelcounts_norm,'stacked')
suptitle('Relative proportion of each suppression model across areas')
text(1,0.65,['m1:' num2str(modelcounts_norm(1,1)*100,3) '%'],'HorizontalAlignment','center')
text(1,0.6,['m2:' num2str(modelcounts_norm(1,2)*100,3) '%'],'HorizontalAlignment','center')
text(1,0.5,['n_{cell}=' num2str(nCells_area(1))],'HorizontalAlignment','center')
text(1,0.45,['n_{exp}=' num2str(nExp_area(1))],'HorizontalAlignment','center')
text(2,0.65,['m1:' num2str(modelcounts_norm(2,1)*100,3) '%'],'HorizontalAlignment','center')
text(2,0.6,['m2:' num2str(modelcounts_norm(2,2)*100,3) '%'],'HorizontalAlignment','center')
text(2,0.5,['n_{cell}=' num2str(nCells_area(2))],'HorizontalAlignment','center')
text(2,0.45,['n_{exp}=' num2str(nExp_area(2))],'HorizontalAlignment','center')
text(3,0.65,['m1:' num2str(modelcounts_norm(3,1)*100,3) '%'],'HorizontalAlignment','center')
text(3,0.6,['m2:' num2str(modelcounts_norm(3,2)*100,3) '%'],'HorizontalAlignment','center')
text(3,0.5,['n_{cell}=' num2str(nCells_area(3))],'HorizontalAlignment','center')
text(3,0.45,['n_{exp}=' num2str(nExp_area(3))],'HorizontalAlignment','center')
text(4,0.65,['m1:' num2str(modelcounts_norm(4,1)*100,3) '%'],'HorizontalAlignment','center')
text(4,0.6,['m2:' num2str(modelcounts_norm(4,2)*100,3) '%'],'HorizontalAlignment','center')
text(4,0.5,['n_{cell}=' num2str(nCells_area(4))],'HorizontalAlignment','center')
text(4,0.45,['n_{exp}=' num2str(nExp_area(4))],'HorizontalAlignment','center')
ylabel('Frac. [cells x cons]')
legend('model1','model2','Location','ne')

% now examine counts vs contrast for each area
cons_c = categorical({'0.1' '0.2' '0.4' '0.8'},{'0.1' '0.2' '0.4' '0.8'});

figure(2)
% v1
modelcounts_V1 = [sum(ism1_V1); sum(ism2_V1)]';
modelcounts_V1_norm = modelcounts_V1./sum(modelcounts_V1,2);
subplot(2,2,1)
bar(cons_c,modelcounts_V1_norm,'stacked')
title(['V1 (n=' num2str(nCells_area(1)) ')'])
xlabel('Contrast')
ylabel('Frac. cells')
legend('model1','model2','Location','ne')
% lm
modelcounts_LM = [sum(ism1_LM); sum(ism2_LM)]';
modelcounts_LM_norm = modelcounts_LM./sum(modelcounts_LM,2);
subplot(2,2,2)
bar(cons_c,modelcounts_LM_norm,'stacked')
title(['LM (n=' num2str(nCells_area(2)) ')'])
xlabel('Contrast')
ylabel('Frac. cells')
legend('model1','model2','Location','ne')
% al
modelcounts_AL = [sum(ism1_AL); sum(ism2_AL)]';
modelcounts_AL_norm = modelcounts_AL./sum(modelcounts_AL,2);
subplot(2,2,3)
bar(cons_c,modelcounts_AL_norm,'stacked')
title(['AL (n=' num2str(nCells_area(3)) ')'])
xlabel('Contrast')
ylabel('Frac. cells')
legend('model1','model2','Location','ne')
% pm
modelcounts_PM = [sum(ism1_PM); sum(ism2_PM)]';
modelcounts_PM_norm = modelcounts_PM./sum(modelcounts_PM,2);
subplot(2,2,4)
bar(cons_c,modelcounts_PM_norm,'stacked')
title(['PM (n=' num2str(nCells_area(4)) ')'])
xlabel('Contrast')
ylabel('Frac. cells')
legend('model1','model2','Location','ne')
%suptitle('Relative proportion of each suppression model by contrast')


%% now look at prefSize trends
%v1
prefSize_V1 = reshape([sizeFits_V1.prefSize],size(sizeFits_V1));
%prefMean_V1_all = mean(prefSize_V1,1);
%prefSEM_V1_all = std(prefSize_V1,[],1)./sqrt(nCells_area(1));
for i = 1:nCon
    temp1 = prefSize_V1(find(ism1_V1(:,i)),i);
    temp2 = prefSize_V1(find(ism2_V1(:,i)),i);
    prefMean_V1_m1(i) = mean(temp1,1);
    prefSEM_V1_m1(i) = std(temp1,[],1)./sqrt(length(temp1));
    prefMean_V1_m2(i) = mean(temp2,1);
    prefSEM_V1_m2(i) = std(temp2,[],1)./sqrt(length(temp2));
    n_V1(i,:) = [length(temp1) length(temp2)];
end
%lm
prefSize_LM = reshape([sizeFits_LM.prefSize],size(sizeFits_LM));
for i = 1:nCon
    temp1 = prefSize_LM(find(ism1_LM(:,i)),i);
    temp2 = prefSize_LM(find(ism2_LM(:,i)),i);
    prefMean_LM_m1(i) = mean(temp1,1);
    prefSEM_LM_m1(i) = std(temp1,[],1)./sqrt(length(temp1));
    prefMean_LM_m2(i) = mean(temp2,1);
    prefSEM_LM_m2(i) = std(temp2,[],1)./sqrt(length(temp2));
    n_LM(i,:) = [length(temp1) length(temp2)];
end
%al
prefSize_AL = reshape([sizeFits_AL.prefSize],size(sizeFits_AL));
for i = 1:nCon
    temp1 = prefSize_AL(find(ism1_AL(:,i)),i);
    temp2 = prefSize_AL(find(ism2_AL(:,i)),i);
    prefMean_AL_m1(i) = mean(temp1,1);
    prefSEM_AL_m1(i) = std(temp1,[],1)./sqrt(length(temp1));
    prefMean_AL_m2(i) = mean(temp2,1);
    prefSEM_AL_m2(i) = std(temp2,[],1)./sqrt(length(temp2));
    n_AL(i,:) = [length(temp1) length(temp2)];
end
%pm
prefSize_PM = reshape([sizeFits_PM.prefSize],size(sizeFits_PM));
for i = 1:nCon
    temp1 = prefSize_PM(find(ism1_PM(:,i)),i);
    temp2 = prefSize_PM(find(ism2_PM(:,i)),i);
    prefMean_PM_m1(i) = mean(temp1,1);
    prefSEM_PM_m1(i) = std(temp1,[],1)./sqrt(length(temp1));
    prefMean_PM_m2(i) = mean(temp2,1);
    prefSEM_PM_m2(i) = std(temp2,[],1)./sqrt(length(temp2));
    n_PM(i,:) = [length(temp1) length(temp2)];
end

% do a subplot for each area, show mean prefSize vs con with errorbars, overlay model1/2
figure(3);clf;
%v1
ax1=subplot(2,2,1);
%errorbar(cons-0.001,prefMean_V1_all,prefSEM_V1_all,'k')
hold on
errorbar(cons,prefMean_V1_m1,prefSEM_V1_m1,'s')
errorbar(cons,prefMean_V1_m2,prefSEM_V1_m2,'^')
hold off
title('V1')
xlabel('Contrast')
ylabel('PrefSize (deg)')
legend('model1','model2','Location','ne')
%lm
ax2=subplot(2,2,2);
hold on
errorbar(cons,prefMean_LM_m1,prefSEM_LM_m1,'s')
errorbar(cons,prefMean_LM_m2,prefSEM_LM_m2,'^')
hold off
title('LM')
xlabel('Contrast')
ylabel('PrefSize (deg)')
legend('model1','model2','Location','ne')
%al
ax3=subplot(2,2,3);
hold on
errorbar(cons,prefMean_AL_m1,prefSEM_AL_m1,'s')
errorbar(cons,prefMean_AL_m2,prefSEM_AL_m2,'^')
hold off
title('AL')
xlabel('Contrast')
ylabel('PrefSize (deg)')
legend('model1','model2','Location','ne')
%pm
ax4=subplot(2,2,4);
hold on
errorbar(cons,prefMean_PM_m1,prefSEM_PM_m1,'s')
errorbar(cons,prefMean_PM_m2,prefSEM_PM_m2,'^')
hold off
title('PM')
xlabel('Contrast')
ylabel('PrefSize (deg)')
legend('model1','model2','Location','ne')
linkaxes([ax1 ax2 ax3 ax4])
xlim([0 1])
ylim([0 60])
suptitle('Preferred Size vs Contrast across Visual Areas')


%% repeat for suppInd
% currently exluding SI values >2 and <-1 to prevent major data skewing
SI_up = 1;
SI_low = 0;
%v1
suppInd_V1 = reshape([sizeFits_V1.suppInd],size(sizeFits_V1));
%suppMean_V1_all = mean(suppInd_V1,1);
%suppSEM_V1_all = std(suppInd_V1,[],1)./sqrt(nCells_area(1));
for i = 1:nCon
    temp1 = suppInd_V1(find(ism1_V1(:,i)),i);
    temp2 = suppInd_V1(find(ism2_V1(:,i)),i);
    %temp2(temp2>SI_up)=NaN; temp2(temp2<SI_low)=NaN;
    temp2(temp2>SI_up)=SI_up; temp2(temp2<SI_low)=NaN;
    suppMean_V1_m1(i) = mean(temp1,1);
    suppSEM_V1_m1(i) = std(temp1,[],1)./sqrt(length(temp1));
    suppMean_V1_m2(i) = nanmean(temp2,1);
    suppSEM_V1_m2(i) = nanstd(temp2,[],1)./sqrt(sum(~isnan(temp2)));
    n_V1(i,:) = [length(temp1) sum(~isnan(temp2))];
    suppMean_V1_all(i) = nanmean([temp1; temp2],1);
    suppSEM_V1_all(i) = nanstd([temp1; temp2],[],1)./sqrt(sum(~isnan([temp1; temp2])));
end
%lm
suppInd_LM = reshape([sizeFits_LM.suppInd],size(sizeFits_LM));
for i = 1:nCon
    temp1 = suppInd_LM(find(ism1_LM(:,i)),i);
    temp2 = suppInd_LM(find(ism2_LM(:,i)),i);
    %temp2(temp2>SI_up)=NaN; temp2(temp2<SI_low)=NaN;
    temp2(temp2>SI_up)=SI_up; temp2(temp2<SI_low)=NaN;
    suppMean_LM_m1(i) = mean(temp1,1);
    suppSEM_LM_m1(i) = std(temp1,[],1)./sqrt(length(temp1));
    suppMean_LM_m2(i) = nanmean(temp2,1);
    suppSEM_LM_m2(i) = nanstd(temp2,[],1)./sqrt(sum(~isnan(temp2)));
    n_LM(i,:) = [length(temp1) sum(~isnan(temp2))];
    suppMean_LM_all(i) = nanmean([temp1; temp2],1);
    suppSEM_LM_all(i) = nanstd([temp1; temp2],[],1)./sqrt(sum(~isnan([temp1; temp2])));
end
%al
suppInd_AL = reshape([sizeFits_AL.suppInd],size(sizeFits_AL));
for i = 1:nCon
    temp1 = suppInd_AL(find(ism1_AL(:,i)),i);
    temp2 = suppInd_AL(find(ism2_AL(:,i)),i);
    %temp2(temp2>SI_up)=NaN; temp2(temp2<SI_low)=NaN;
    temp2(temp2>SI_up)=SI_up; temp2(temp2<SI_low)=NaN;
    suppMean_AL_m1(i) = mean(temp1,1);
    suppSEM_AL_m1(i) = std(temp1,[],1)./sqrt(length(temp1));
    suppMean_AL_m2(i) = nanmean(temp2,1);
    suppSEM_AL_m2(i) = nanstd(temp2,[],1)./sqrt(sum(~isnan(temp2)));
    n_AL(i,:) = [length(temp1) sum(~isnan(temp2))];
    suppMean_AL_all(i) = nanmean([temp1; temp2],1);
    suppSEM_AL_all(i) = nanstd([temp1; temp2],[],1)./sqrt(sum(~isnan([temp1; temp2])));
end
%pm
suppInd_PM = reshape([sizeFits_PM.suppInd],size(sizeFits_PM));
for i = 1:nCon
    temp1 = suppInd_PM(find(ism1_PM(:,i)),i);
    temp2 = suppInd_PM(find(ism2_PM(:,i)),i);
    %temp2(temp2>SI_up)=NaN; temp2(temp2<SI_low)=NaN;
    temp2(temp2>SI_up)=SI_up; temp2(temp2<SI_low)=NaN;
    suppMean_PM_m1(i) = mean(temp1,1);
    suppSEM_PM_m1(i) = std(temp1,[],1)./sqrt(length(temp1));
    suppMean_PM_m2(i) = nanmean(temp2,1);
    suppSEM_PM_m2(i) = nanstd(temp2,[],1)./sqrt(sum(~isnan(temp2)));
    n_PM(i,:) = [length(temp1) sum(~isnan(temp2))];
    suppMean_PM_all(i) = nanmean([temp1; temp2],1);
    suppSEM_PM_all(i) = nanstd([temp1; temp2],[],1)./sqrt(sum(~isnan([temp1; temp2])));
end

% do a subplot for each area, show mean suppInd vs con with errorbars, overlay model1/2
figure(4);clf;
suptitle('Suppression Index vs Contrast across Visual Areas')
%v1
ax1=subplot(2,2,1);
%errorbar(cons-0.001,suppMean_V1_all,suppSEM_V1_all,'k')
hold on
errorbar(cons,suppMean_V1_m1,suppSEM_V1_m1,'s')
errorbar(cons,suppMean_V1_m2,suppSEM_V1_m2,'^')
errorbar(cons,suppMean_V1_all,suppSEM_V1_all,'kx')
hold off
title(['V1 (n=' num2str(max(sum(n_V1'))) ')'])
xlabel('Contrast')
ylabel('SI')
legend('model1','model2','all','Location','best')
%lm
ax2=subplot(2,2,2);
hold on
errorbar(cons,suppMean_LM_m1,suppSEM_LM_m1,'s')
errorbar(cons,suppMean_LM_m2,suppSEM_LM_m2,'^')
errorbar(cons,suppMean_LM_all,suppSEM_LM_all,'kx')
hold off
title(['LM (n=' num2str(max(sum(n_LM'))) ')'])
xlabel('Contrast')
ylabel('SI')
legend('model1','model2','all','Location','best')
%al
ax3=subplot(2,2,3);
hold on
errorbar(cons,suppMean_AL_m1,suppSEM_AL_m1,'s')
errorbar(cons,suppMean_AL_m2,suppSEM_AL_m2,'^')
errorbar(cons,suppMean_AL_all,suppSEM_AL_all,'kx')
hold off
title(['AL (n=' num2str(max(sum(n_AL'))) ')'])
xlabel('Contrast')
ylabel('SI')
legend('model1','model2','all','Location','best')
%pm
ax4=subplot(2,2,4);
hold on
errorbar(cons,suppMean_PM_m1,suppSEM_PM_m1,'s')
errorbar(cons,suppMean_PM_m2,suppSEM_PM_m2,'^')
errorbar(cons,suppMean_PM_all,suppSEM_PM_all,'kx')
hold off
title(['PM (n=' num2str(max(sum(n_PM'))) ')'])
xlabel('Contrast')
ylabel('SI')
legend('model1','model2','all','Location','best')
linkaxes([ax1 ax2 ax3 ax4])
xlim([0 1])
ylim([0 1])


%% examine by GECI

% indicator experiments + cell indices
slowExp = find(cellfun(@(x) strcmp(x,'s'), expdata.GECI, 'UniformOutput', 1));
fastExp = find(cellfun(@(x) strcmp(x,'f'), expdata.GECI, 'UniformOutput', 1));
slowInd = find(ismember(expInd,slowExp));
fastInd = find(ismember(expInd,fastExp));

% area x indicator indices
V1ind_s = intersect(V1ind, slowInd);
V1ind_f = intersect(V1ind, fastInd);
LMind_s = intersect(LMind, slowInd);
LMind_f = intersect(LMind, fastInd);
ALind_s = intersect(ALind, slowInd);
ALind_f = intersect(ALind, fastInd);
PMind_s = intersect(PMind, slowInd);
PMind_f = intersect(PMind, fastInd);

% extract for GCaMP6s
%v1_s
sizeTune_V1_s = sizeTune_all{:,:,V1ind_s}; % (size,con,cell)
sizeMean_V1_s = sizeMean_all(:,:,V1ind_s);
sizeSEM_V1_s = sizeSEM_all(:,:,V1ind_s);
sizeFits_V1_s = sizeFits_all(V1ind_s,:); %cell,con
lbub_fits_V1_s = lbub_fits_all(V1ind_s,:,:); %cell,par,val (low up mean true stdev)
ism1_V1_s = reshape(~[sizeFits_V1_s.Ftest],size(sizeFits_V1_s));
ism2_V1_s = reshape([sizeFits_V1_s.Ftest],size(sizeFits_V1_s));
%ism1_V1 = find(~lbub_fits_V1(:,8,4));
%ism2_V1 = find(lbub_fits_V1(:,8,4));
%lm_s
sizeTune_LM_s = sizeTune_all{:,:,LMind_s}; % (size,con,cell)
sizeMean_LM_s = sizeMean_all(:,:,LMind_s);
sizeSEM_LM_s = sizeSEM_all(:,:,LMind_s);
sizeFits_LM_s = sizeFits_all(LMind_s,:); %cell,con
lbub_fits_LM_s = lbub_fits_all(LMind_s,:,:); %cell,par,val (low up mean true stdev)
ism1_LM_s = reshape(~[sizeFits_LM_s.Ftest],size(sizeFits_LM_s));
ism2_LM_s = reshape([sizeFits_LM_s.Ftest],size(sizeFits_LM_s));
%al_s
sizeTune_AL_s = sizeTune_all{:,:,ALind_s}; % (size,con,cell)
sizeMean_AL_s = sizeMean_all(:,:,ALind_s);
sizeSEM_AL_s = sizeSEM_all(:,:,ALind_s);
sizeFits_AL_s = sizeFits_all(ALind_s,:); %cell,con
lbub_fits_AL_s = lbub_fits_all(ALind_s,:,:); %cell,par,val (low up mean true stdev)
ism1_AL_s = reshape(~[sizeFits_AL_s.Ftest],size(sizeFits_AL_s));
ism2_AL_s = reshape([sizeFits_AL_s.Ftest],size(sizeFits_AL_s));
%pm_s
sizeTune_PM_s = sizeTune_all{:,:,PMind_s}; % (size,con,cell)
sizeMean_PM_s = sizeMean_all(:,:,PMind_s);
sizeSEM_PM_s = sizeSEM_all(:,:,PMind_s);
sizeFits_PM_s = sizeFits_all(PMind_s,:); %cell,con
lbub_fits_PM_s = lbub_fits_all(PMind_s,:,:); %cell,par,val (low up mean true stdev)
ism1_PM_s = reshape(~[sizeFits_PM_s.Ftest],size(sizeFits_PM_s));
ism2_PM_s = reshape([sizeFits_PM_s.Ftest],size(sizeFits_PM_s));
% again for GCaMP6f
%v1_f
sizeTune_V1_f = sizeTune_all{:,:,V1ind_f}; % (size,con,cell)
sizeMean_V1_f = sizeMean_all(:,:,V1ind_f);
sizeSEM_V1_f = sizeSEM_all(:,:,V1ind_f);
sizeFits_V1_f = sizeFits_all(V1ind_f,:); %cell,con
lbub_fits_V1_f = lbub_fits_all(V1ind_f,:,:); %cell,par,val (low up mean true stdev)
ism1_V1_f = reshape(~[sizeFits_V1_f.Ftest],size(sizeFits_V1_f));
ism2_V1_f = reshape([sizeFits_V1_f.Ftest],size(sizeFits_V1_f));
%ism1_V1 = find(~lbub_fits_V1(:,8,4));
%ism2_V1 = find(lbub_fits_V1(:,8,4));
%lm_f
sizeTune_LM_f = sizeTune_all{:,:,LMind_f}; % (size,con,cell)
sizeMean_LM_f = sizeMean_all(:,:,LMind_f);
sizeSEM_LM_f = sizeSEM_all(:,:,LMind_f);
sizeFits_LM_f = sizeFits_all(LMind_f,:); %cell,con
lbub_fits_LM_f = lbub_fits_all(LMind_f,:,:); %cell,par,val (low up mean true stdev)
ism1_LM_f = reshape(~[sizeFits_LM_f.Ftest],size(sizeFits_LM_f));
ism2_LM_f = reshape([sizeFits_LM_f.Ftest],size(sizeFits_LM_f));
%al_f
sizeTune_AL_f = sizeTune_all{:,:,ALind_f}; % (size,con,cell)
sizeMean_AL_f = sizeMean_all(:,:,ALind_f);
sizeSEM_AL_f = sizeSEM_all(:,:,ALind_f);
sizeFits_AL_f = sizeFits_all(ALind_f,:); %cell,con
lbub_fits_AL_f = lbub_fits_all(ALind_f,:,:); %cell,par,val (low up mean true stdev)
ism1_AL_f = reshape(~[sizeFits_AL_f.Ftest],size(sizeFits_AL_f));
ism2_AL_f = reshape([sizeFits_AL_f.Ftest],size(sizeFits_AL_f));
%pm_f
sizeTune_PM_f = sizeTune_all{:,:,PMind_f}; % (size,con,cell)
sizeMean_PM_f = sizeMean_all(:,:,PMind_f);
sizeSEM_PM_f = sizeSEM_all(:,:,PMind_f);
sizeFits_PM_f = sizeFits_all(PMind_f,:); %cell,con
lbub_fits_PM_f = lbub_fits_all(PMind_f,:,:); %cell,par,val (low up mean true stdev)
ism1_PM_f = reshape(~[sizeFits_PM_f.Ftest],size(sizeFits_PM_f));
ism2_PM_f = reshape([sizeFits_PM_f.Ftest],size(sizeFits_PM_f));

% total exps and cells in each area/GECI
nExp_area_s = [length(intersect(slowExp,V1exp)) length(intersect(slowExp,LMexp)) length(intersect(slowExp,ALexp)) length(intersect(slowExp,PMexp))];
nCells_area_s = [length(V1ind_s) length(LMind_s) length(ALind_s) length(PMind_s)];
nExp_area_f = [length(intersect(fastExp,V1exp)) length(intersect(fastExp,LMexp)) length(intersect(fastExp,ALexp)) length(intersect(fastExp,PMexp))];
nCells_area_f = [length(V1ind_f) length(LMind_f) length(ALind_f) length(PMind_f)];

%% look at size tuning curves
figure(7);clf;
for q=1:4 % quadrants = (V1,LM,AL,PM)
    switch q
        case 1
            area = 'V1';
            offset = [1 2];
            sizeMean_f = sizeMean_V1_f;
            sizeSEM_f = sizeSEM_V1_f;
            ism1_f = find(ism1_V1_f(:,nCon));
            ism2_f = find(ism2_V1_f(:,nCon));
            sizeMean_s = sizeMean_V1_s;
            sizeSEM_s = sizeSEM_V1_s;
            ism1_s = find(ism1_V1_s(:,nCon));
            ism2_s = find(ism2_V1_s(:,nCon));
        case 2
            area = 'LM';
            offset = [3 4];
            sizeMean_f = sizeMean_LM_f;
            sizeSEM_f = sizeSEM_LM_f;
            ism1_f = find(ism1_LM_f(:,nCon));
            ism2_f = find(ism2_LM_f(:,nCon));
            sizeMean_s = sizeMean_LM_s;
            sizeSEM_s = sizeSEM_LM_s;
            ism1_s = find(ism1_LM_s(:,nCon));
            ism2_s = find(ism2_LM_s(:,nCon));
        case 3
            area = 'AL';
            offset = [5 6];
            sizeMean_f = sizeMean_AL_f;
            sizeSEM_f = sizeSEM_AL_f;
            ism1_f = find(ism1_AL_f(:,nCon));
            ism2_f = find(ism2_AL_f(:,nCon));
            sizeMean_s = sizeMean_AL_s;
            sizeSEM_s = sizeSEM_AL_s;
            ism1_s = find(ism1_AL_s(:,nCon));
            ism2_s = find(ism2_AL_s(:,nCon));
        case 4
            area = 'PM';
            offset = [7 8];
            sizeMean_f = sizeMean_PM_f;
            sizeSEM_f = sizeSEM_PM_f;
            ism1_f = find(ism1_PM_f(:,nCon));
            ism2_f = find(ism2_PM_f(:,nCon));
            sizeMean_s = sizeMean_PM_s;
            sizeSEM_s = sizeSEM_PM_s;
            ism1_s = find(ism1_PM_s(:,nCon));
            ism2_s = find(ism2_PM_s(:,nCon));
    end
    
    norm_f = max(mean(sizeMean_f(:,nCon,:),3));
    norm_s = max(mean(sizeMean_s(:,nCon,:),3));
    %     %norm_f = max(mean(sizeMean_f(:,nCon,ism1_f),3));
    %     %norm_s = max(mean(sizeMean_s(:,nCon,ism1_s),3));
    %     subplot(2,4,offset(1))
    %     hold on
    %     errorbar(szs,mean(sizeMean_f(:,nCon,ism1_f),3)/norm_f,geomean(sizeSEM_f(:,nCon,ism1_f),3)/norm_f,'r')
    %     errorbar(szs,mean(sizeMean_s(:,nCon,ism1_s),3)/norm_s,geomean(sizeSEM_s(:,nCon,ism1_s),3)/norm_s,'b')
    %     title([area ' model1'])
    %     xlabel('Size (deg)')
    %     ylabel('dF/F (normalized)')
    %     legend(['fast (n=' num2str(length(ism1_f)) ')'], ['slow (n=' num2str(length(ism1_s)) ')'])
    %
    %     %norm_f = max(mean(sizeMean_f(:,nCon,ism2_f),3));
    %     %norm_s = max(mean(sizeMean_s(:,nCon,ism2_s),3));
    %     subplot(2,4,offset(2))
    %     hold on
    %     errorbar(szs,mean(sizeMean_f(:,nCon,ism2_f),3)/norm_f,geomean(sizeSEM_f(:,nCon,ism2_f),3)/norm_f,'r')
    %     errorbar(szs,mean(sizeMean_s(:,nCon,ism2_s),3)/norm_s,geomean(sizeSEM_s(:,nCon,ism2_s),3)/norm_s,'b')
    %     title([area ' model2'])
    %     xlabel('Size (deg)')
    %     ylabel('dF/F')
    %     legend(['fast (n=' num2str(length(ism2_f)) ')'], ['slow (n=' num2str(length(ism2_s)) ')'])
    %
    subplot(2,2,q)
    hold on
    errorbar(szs,mean(sizeMean_f(:,nCon,:),3)/norm_f,geomean(sizeSEM_f(:,nCon,:),3)/norm_f,'r')
    errorbar(szs,mean(sizeMean_s(:,nCon,:),3)/norm_s,geomean(sizeSEM_s(:,nCon,:),3)/norm_s,'b')
    title([area ' mean size-tuning curve'])
    xlabel('Size (deg)')
    ylabel('dF/F')
    legend(['fast (n=' num2str(size(sizeMean_f,3)) ')'], ['slow (n=' num2str(size(sizeMean_s,3)) ')'],'location','best')
end

% look at proportion of model1 vs model2
% % first examine cell counts, proportion of each model
% figure(1);clf;
% modelcounts = [sum(~[sizeFits_V1.Ftest]) sum([sizeFits_V1.Ftest]); ...
%     sum(~[sizeFits_LM.Ftest]) sum([sizeFits_LM.Ftest]); ...
%     sum(~[sizeFits_AL.Ftest]) sum([sizeFits_AL.Ftest]); ...
%     sum(~[sizeFits_PM.Ftest]) sum([sizeFits_PM.Ftest])];
% modelcounts_norm = modelcounts./sum(modelcounts,2);
% bar(areas_c,modelcounts_norm,'stacked')
% suptitle('Relative proportion of each suppression model across areas')
% text(1,0.65,['m1:' num2str(modelcounts_norm(1,1)*100,3) '%'],'HorizontalAlignment','center')
% text(1,0.6,['m2:' num2str(modelcounts_norm(1,2)*100,3) '%'],'HorizontalAlignment','center')
% text(1,0.5,['n_{cell}=' num2str(nCells_area(1))],'HorizontalAlignment','center')
% text(1,0.45,['n_{exp}=' num2str(nExp_area(1))],'HorizontalAlignment','center')
% text(2,0.65,['m1:' num2str(modelcounts_norm(2,1)*100,3) '%'],'HorizontalAlignment','center')
% text(2,0.6,['m2:' num2str(modelcounts_norm(2,2)*100,3) '%'],'HorizontalAlignment','center')
% text(2,0.5,['n_{cell}=' num2str(nCells_area(2))],'HorizontalAlignment','center')
% text(2,0.45,['n_{exp}=' num2str(nExp_area(2))],'HorizontalAlignment','center')
% text(3,0.65,['m1:' num2str(modelcounts_norm(3,1)*100,3) '%'],'HorizontalAlignment','center')
% text(3,0.6,['m2:' num2str(modelcounts_norm(3,2)*100,3) '%'],'HorizontalAlignment','center')
% text(3,0.5,['n_{cell}=' num2str(nCells_area(3))],'HorizontalAlignment','center')
% text(3,0.45,['n_{exp}=' num2str(nExp_area(3))],'HorizontalAlignment','center')
% text(4,0.65,['m1:' num2str(modelcounts_norm(4,1)*100,3) '%'],'HorizontalAlignment','center')
% text(4,0.6,['m2:' num2str(modelcounts_norm(4,2)*100,3) '%'],'HorizontalAlignment','center')
% text(4,0.5,['n_{cell}=' num2str(nCells_area(4))],'HorizontalAlignment','center')
% text(4,0.45,['n_{exp}=' num2str(nExp_area(4))],'HorizontalAlignment','center')
% ylabel('Frac. [cells x cons]')
% legend('model1','model2','Location','ne')

% now examine counts vs contrast in each area, comparing fast vs slow GECI
cons_c = categorical({'0.1' '0.2' '0.4' '0.8'},{'0.1' '0.2' '0.4' '0.8'});

figure(8);clf;
colormap([1 0 0;0 0 1])
% v1
modelcounts_V1 = [sum(ism2_V1_f)/nCells_area_f(1); sum(ism2_V1_s)/nCells_area_s(1)]';
subplot(2,2,1)
bar(cons_c,modelcounts_V1,'grouped')
title(['V1 (n_f=' num2str(nCells_area_f(1)) ', n_s=' num2str(nCells_area_s(1)) ')'])
xlabel('Contrast')
ylabel('Frac. model2 cells')
legend('f','s','Location','se')
% lm
modelcounts_LM = [sum(ism2_LM_f)/nCells_area_f(2); sum(ism2_LM_s)/nCells_area_s(2)]';
subplot(2,2,2)
bar(cons_c,modelcounts_LM,'grouped')
title(['LM (n_f=' num2str(nCells_area_f(2)) ', n_s=' num2str(nCells_area_s(2)) ')'])
xlabel('Contrast')
ylabel('Frac. model2 cells')
legend('f','s','Location','se')
% al
modelcounts_AL = [sum(ism2_AL_f)/nCells_area_f(3); sum(ism2_AL_s)/nCells_area_s(3)]';
subplot(2,2,3)
bar(cons_c,modelcounts_AL,'grouped')
title(['AL (n_f=' num2str(nCells_area_f(3)) ', n_s=' num2str(nCells_area_s(3)) ')'])
xlabel('Contrast')
ylabel('Frac. model2 cells')
legend('f','s','Location','se')
% pm
modelcounts_PM = [sum(ism2_PM_f)/nCells_area_f(4); sum(ism2_PM_s)/nCells_area_s(4)]';
subplot(2,2,4)
bar(cons_c,modelcounts_PM,'grouped')
title(['PM (n_f=' num2str(nCells_area_f(4)) ', n_s=' num2str(nCells_area_s(4)) ')'])
xlabel('Contrast')
ylabel('Frac. model2 cells')
legend('f','s','Location','nw')
ylim([0 1])


%% now look at prefSize trends
prefSize_V1_f = reshape([sizeFits_V1_f.prefSize],size(sizeFits_V1_f));
prefSize_LM_f = reshape([sizeFits_LM_f.prefSize],size(sizeFits_LM_f));
prefSize_AL_f = reshape([sizeFits_AL_f.prefSize],size(sizeFits_AL_f));
prefSize_PM_f = reshape([sizeFits_PM_f.prefSize],size(sizeFits_PM_f));
prefSize_V1_s = reshape([sizeFits_V1_s.prefSize],size(sizeFits_V1_s));
prefSize_LM_s = reshape([sizeFits_LM_s.prefSize],size(sizeFits_LM_s));
prefSize_AL_s = reshape([sizeFits_AL_s.prefSize],size(sizeFits_AL_s));
prefSize_PM_s = reshape([sizeFits_PM_s.prefSize],size(sizeFits_PM_s));
%prefMean_V1_all = mean(prefSize_V1,1);
%prefSEM_V1_all = std(prefSize_V1,[],1)./sqrt(nCells_area(1));
for i = 1:nCon
    %v1
    prefMean_V1_f(i) = mean(prefSize_V1_f(:,i),1);
    prefSEM_V1_f(i) = std(prefSize_V1_f(:,i),[],1)./sqrt(nCells_area_f(1));
    prefMean_V1_s(i) = mean(prefSize_V1_s(:,i),1);
    prefSEM_V1_s(i) = std(prefSize_V1_s(:,i),[],1)./sqrt(nCells_area_s(1));
    %lm
    prefMean_LM_f(i) = mean(prefSize_LM_f(:,i),1);
    prefSEM_LM_f(i) = std(prefSize_LM_f(:,i),[],1)./sqrt(nCells_area_f(2));
    prefMean_LM_s(i) = mean(prefSize_LM_s(:,i),1);
    prefSEM_LM_s(i) = std(prefSize_LM_s(:,i),[],1)./sqrt(nCells_area_s(2));
    %al
    prefMean_AL_f(i) = mean(prefSize_AL_f(:,i),1);
    prefSEM_AL_f(i) = std(prefSize_AL_f(:,i),[],1)./sqrt(nCells_area_f(3));
    prefMean_AL_s(i) = mean(prefSize_AL_s(:,i),1);
    prefSEM_AL_s(i) = std(prefSize_AL_s(:,i),[],1)./sqrt(nCells_area_s(3));
    %pm
    prefMean_PM_f(i) = mean(prefSize_PM_f(:,i),1);
    prefSEM_PM_f(i) = std(prefSize_PM_f(:,i),[],1)./sqrt(nCells_area_f(4));
    prefMean_PM_s(i) = mean(prefSize_PM_s(:,i),1);
    prefSEM_PM_s(i) = std(prefSize_PM_s(:,i),[],1)./sqrt(nCells_area_s(4));
end

% do a subplot for each area, show mean prefSize vs con with errorbars, overlay model1/2
figure(3);clf;
colormap([1 0 0;0 0 1])
%suptitle('Preferred Size vs Contrast across Visual Areas')
%v1
ax1=subplot(2,2,1);
%errorbar(cons-0.001,prefMean_V1_all,prefSEM_V1_all,'k')
hold on
errorbar(cons,prefMean_V1_f,prefSEM_V1_f,'s-')
errorbar(cons,prefMean_V1_s,prefSEM_V1_s,'^-')
hold off
title(['V1 (n_f=' num2str(nCells_area_f(1)) ', n_s=' num2str(nCells_area_s(1)) ')'])
xlabel('Contrast')
ylabel('PrefSize (deg)')
legend('f','s','Location','nw')
%lm
ax2=subplot(2,2,2);
hold on
errorbar(cons,prefMean_LM_f,prefSEM_LM_f,'s-')
errorbar(cons,prefMean_LM_s,prefSEM_LM_s,'^-')
hold off
title(['LM (n_f=' num2str(nCells_area_f(2)) ', n_s=' num2str(nCells_area_s(2)) ')'])
xlabel('Contrast')
ylabel('PrefSize (deg)')
legend('f','s','Location','nw')
%al
ax3=subplot(2,2,3);
hold on
errorbar(cons,prefMean_AL_f,prefSEM_AL_f,'s-')
errorbar(cons,prefMean_AL_s,prefSEM_AL_s,'^-')
hold off
title(['AL (n_f=' num2str(nCells_area_f(3)) ', n_s=' num2str(nCells_area_s(3)) ')'])
xlabel('Contrast')
ylabel('PrefSize (deg)')
legend('f','s','Location','nw')
%pm
ax4=subplot(2,2,4);
hold on
errorbar(cons,prefMean_PM_f,prefSEM_PM_f,'s-')
errorbar(cons,prefMean_PM_s,prefSEM_PM_s,'^-')
hold off
title(['PM (n_f=' num2str(nCells_area_f(4)) ', n_s=' num2str(nCells_area_s(4)) ')'])
xlabel('Contrast')
ylabel('PrefSize (deg)')
legend('f','s','Location','nw')
linkaxes([ax1 ax2 ax3 ax4])
xlim([0 1])
ylim([0 70])

%% repeat for suppInd
% currently exluding SI values >2 and <-1 to prevent major data skewing
SI_up = 1;
SI_low = 0;

suppInd_V1_f = reshape([sizeFits_V1_f.suppInd],size(sizeFits_V1_f));
suppInd_LM_f = reshape([sizeFits_LM_f.suppInd],size(sizeFits_LM_f));
suppInd_AL_f = reshape([sizeFits_AL_f.suppInd],size(sizeFits_AL_f));
suppInd_PM_f = reshape([sizeFits_PM_f.suppInd],size(sizeFits_PM_f));
suppInd_V1_s = reshape([sizeFits_V1_s.suppInd],size(sizeFits_V1_s));
suppInd_LM_s = reshape([sizeFits_LM_s.suppInd],size(sizeFits_LM_s));
suppInd_AL_s = reshape([sizeFits_AL_s.suppInd],size(sizeFits_AL_s));
suppInd_PM_s = reshape([sizeFits_PM_s.suppInd],size(sizeFits_PM_s));
%suppMean_V1_all = mean(suppInd_V1,1);
%suppSEM_V1_all = std(suppInd_V1,[],1)./sqrt(nCells_area(1));
for i = 1:nCon
    temp1 = suppInd_V1_f(:,i);
    temp2 = suppInd_V1_s(:,i);
    temp1(temp1>SI_up)=SI_up; temp1(temp1<SI_low)=SI_low;
    temp2(temp2>SI_up)=SI_up; temp2(temp2<SI_low)=SI_low;
    suppMean_V1_f(i) = nanmean(temp1,1);
    suppSEM_V1_f(i) = nanstd(temp1,[],1)./sqrt(sum(~isnan(temp1)));
    suppMean_V1_s(i) = nanmean(temp2,1);
    suppSEM_V1_s(i) = nanstd(temp2,[],1)./sqrt(sum(~isnan(temp2)));
    %n_V1(i,:) = [length(temp1) length(temp2)];
    %lm
    temp1 = suppInd_LM_f(:,i);
    temp2 = suppInd_LM_s(:,i);
    temp1(temp1>SI_up)=SI_up; temp1(temp1<SI_low)=SI_low;
    temp2(temp2>SI_up)=SI_up; temp2(temp2<SI_low)=SI_low;
    suppMean_LM_f(i) = nanmean(temp1,1);
    suppSEM_LM_f(i) = nanstd(temp1,[],1)./sqrt(sum(~isnan(temp1)));
    suppMean_LM_s(i) = nanmean(temp2,1);
    suppSEM_LM_s(i) = nanstd(temp2,[],1)./sqrt(sum(~isnan(temp2)));
    %al
    temp1 = suppInd_AL_f(:,i);
    temp2 = suppInd_AL_s(:,i);
    temp1(temp1>SI_up)=SI_up; temp1(temp1<SI_low)=SI_low;
    temp2(temp2>SI_up)=SI_up; temp2(temp2<SI_low)=SI_low;
    suppMean_AL_f(i) = nanmean(temp1,1);
    suppSEM_AL_f(i) = nanstd(temp1,[],1)./sqrt(sum(~isnan(temp1)));
    suppMean_AL_s(i) = nanmean(temp2,1);
    suppSEM_AL_s(i) = nanstd(temp2,[],1)./sqrt(sum(~isnan(temp2)));
    %pm
    temp1 = suppInd_PM_f(:,i);
    temp2 = suppInd_PM_s(:,i);
    temp1(temp1>SI_up)=SI_up; temp1(temp1<SI_low)=SI_low;
    temp2(temp2>SI_up)=SI_up; temp2(temp2<SI_low)=SI_low;
    suppMean_PM_f(i) = nanmean(temp1,1);
    suppSEM_PM_f(i) = nanstd(temp1,[],1)./sqrt(sum(~isnan(temp1)));
    suppMean_PM_s(i) = nanmean(temp2,1);
    suppSEM_PM_s(i) = nanstd(temp2,[],1)./sqrt(sum(~isnan(temp2)));
end

% do a subplot for each area, show mean suppInd vs con with errorbars, overlay model1/2
figure(4);clf;
colormap([1 0 0;0 0 1])
%suptitle('Suppression Index vs Contrast across Visual Areas')
%v1
ax1=subplot(2,2,1);
bar(cons_c,[suppMean_V1_f; suppMean_V1_s;]','grouped')
title(['V1 (n_f=' num2str(nCells_area_f(1)) ', n_s=' num2str(nCells_area_s(1)) ')'])
xlabel('Contrast')
ylabel('SI')
legend('f','s','Location','se')
%lm
ax2=subplot(2,2,2);
bar(cons_c,[suppMean_LM_f; suppMean_LM_s;]','grouped')
title(['LM (n_f=' num2str(nCells_area_f(2)) ', n_s=' num2str(nCells_area_s(2)) ')'])
xlabel('Contrast')
ylabel('SI')
legend('f','s','Location','se')
%al
ax3=subplot(2,2,3);
bar(cons_c,[suppMean_AL_f; suppMean_AL_s;]','grouped')
title(['AL (n_f=' num2str(nCells_area_f(3)) ', n_s=' num2str(nCells_area_s(3)) ')'])
xlabel('Contrast')
ylabel('SI')
legend('f','s','Location','se')
%pm
ax4=subplot(2,2,4);
bar(cons_c,[suppMean_PM_f; suppMean_PM_s;]','grouped')
title(['PM (n_f=' num2str(nCells_area_f(4)) ', n_s=' num2str(nCells_area_s(4)) ')'])
xlabel('Contrast')
ylabel('SI')
legend('f','s','Location','nw')
linkaxes([ax1 ax2 ax3 ax4])
xlim([0 1])
ylim([0 1])

%% extract unique cells from full set
% plots figures for model fraction, average curves, prefSize, suppInd
% using the experiment conditions specified below
conds = ["V1" "f";
    "V1" "s";
    "LM" "f";
    "LM" "s";
    "AL" "f";
    "AL" "s";
    "PM" "f";
    "PM" "s";];

expInd = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExp(i))];
end

close all
choosefig = [5:7];
% choose figs: 1=modelcounts; 2=averagecurves; 3=prefSize; 4=suppInd; 5=conresp; 6=ex.cells; 7=medianfits;
legStrs = strings(1,length(conds));

nExp_area = zeros(1,length(conds));
nCells_area = nExp_area;
modelcountsi = nExp_area; prefMeani=nExp_area; prefSEMi=prefMeani;
szs = 5*1.5.^(0:7); nSz=length(szs);
cons = 0.1*2.^(0:3); nCon=length(cons);
c_areas = categorical({'V1' 'LM' 'AL' 'PM' 'all'},{'V1' 'LM' 'AL' 'PM' 'all'});
for i = 1:length(conds)
    % select exps matching cond [area, GECI]
    expArea = find(cellfun(@(x) strcmp(x,conds(i,1)), expdata.area, 'UniformOutput', 1));
    expGECI = find(cellfun(@(x) strcmp(x,conds(i,2)), expdata.GECI, 'UniformOutput', 1));
    expIndi = intersect(expArea,expGECI);
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_all);
    
    % cutoff by cellDist
    switch conds(i,1)
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
        case {'LM','AL'}
            cutoff = 10; %alm cutoff at 15
        case 'PM'
            cutoff = 10; %pm cutoff at 20
    end
    ind = intersect(ind,find(cellDists_all<cutoff));
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    nExp_area(i) = nExpi;
    nCells_area(i) = nCellsi;
    sizeTune = sizeTune_all{:,:,ind}; % (size,con,cell)
    sizeMean = sizeMean_all(:,:,ind);
    sizeSEM = sizeSEM_all(:,:,ind);
    sizeFits = sizeFits_all(ind,:); %cell,con
    lbub_fits = lbub_fits_all(ind,:,:); %cell,par,val (low up mean true stdev)
    ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
    ism2 = reshape([sizeFits.Ftest],size(sizeFits));
    
    cons_c = categorical({'0.1' '0.2' '0.4' '0.8'});
    
    modelcounts = sum(ism2)/nCellsi;
    modelcountsi(i) = modelcounts(nCon);
    if sum(choosefig==1)
        figure(1);if i==1;clf;end %figure 1 = proportions of model2 by con
        subplot(2,4,i)
        bar(cons_c,modelcounts,'grouped')
        title({sprintf('Area:%s, GECI:%s',conds(i,1),conds(i,2));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Contrast')
        ylabel('Frac. model2 cells')
    end
    
    if sum(choosefig==2)
        figure(2);if i==1;clf;end %figure 2 = average size tuning curves (normalized)
        subplot(4,4,2*(i-1)+1)
        sizeMean_norm = sizeMean*0; sizeSEM_norm = sizeSEM*0;
        for iCell = 1:nCellsi
            dum = sizeMean(:,:,iCell); % take all sizeMean values for cell
            norm = max(dum(:)); % take max of all dF/F's including all cons
            sizeMean_norm(:,:,iCell) = sizeMean(:,:,iCell)/norm; % normalize by this max for the individual cell
            sizeSEM_norm(:,:,iCell) = sizeSEM(:,:,iCell)/norm;
        end
        for iCon = 1:nCon
            errorbar(szs,mean(sizeMean_norm(:,iCon,find(ism1(:,iCon))),3),geomean(sizeSEM_norm(:,iCon,find(ism1(:,iCon))),3))
            hold on
        end
        title({sprintf('Model1: Area:%s, GECI:%s',conds(i,1),conds(i,2));['(n=' num2str(mean(sum(ism1))) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Size (deg)')
        ylabel('dF/F (norm)')
        ylim([0 1.2])
        subplot(4,4,2*(i-1)+2)
        for iCon = 1:nCon
            errorbar(szs,mean(sizeMean_norm(:,iCon,find(ism2(:,iCon))),3),geomean(sizeSEM_norm(:,iCon,find(ism2(:,iCon))),3))
            hold on
        end
        title({sprintf('Model2: Area:%s, GECI:%s',conds(i,1),conds(i,2));['(n=' num2str(mean(sum(ism2))) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Size (deg)')
        ylabel('dF/F (norm)')
        ylim([0 1.2])
        if i==8;legend(num2str(cons'));end
    end
    
    prefSize = reshape([sizeFits.prefSize],size(sizeFits));
    prefMean1=zeros(1,nCon);prefSEM1=prefMean1;prefMean2=prefMean1;prefSEM2=prefMean1;prefMeanAll=prefMean1;prefSEMAll=prefMean1;
    for iCon=1:nCon
        prefMean1(iCon) = mean(prefSize(find(ism1(:,iCon)),iCon));
        prefSEM1(iCon) = std(prefSize(find(ism1(:,iCon)),iCon))./sqrt(sum(ism1(:,iCon)));
        prefMean2(iCon) = mean(prefSize(find(ism2(:,iCon)),iCon));
        prefSEM2(iCon) = std(prefSize(find(ism2(:,iCon)),iCon))./sqrt(sum(ism2(:,iCon)));
        prefMeanAll(iCon) = mean(prefSize(:,iCon));
        prefSEMAll(iCon) = std(prefSize(:,iCon))./sqrt(nCellsi);
    end
    prefMeani(i) = prefMeanAll(nCon);
    prefSEMi(i) = prefSEMAll(nCon);
    if sum(choosefig==3)
        figure(3);if i==1;clf;end %figure 3 = prefSize vs con
        subplot(2,4,i)
        errorbar(cons,prefMean1,prefSEM1,'s-');
        hold on
        errorbar(cons,prefMean2,prefSEM2,'^-');
        errorbar(cons,prefMeanAll,prefSEMAll,'kx-');
        hold off
        title({sprintf('Area:%s, GECI:%s',conds(i,1),conds(i,2));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Contrast')
        ylabel('PrefSize')
        ylim([0 max(szs)+1])
        if i==8;legend('m1','m2','all');end
    end
    
    suppInd = reshape([sizeFits.suppInd],size(sizeFits));
    suppInd(suppInd<0)=0;suppInd(suppInd>1)=1;
    suppMean2=zeros(1,nCon);suppSEM2=suppMean2;suppMeanAll=suppMean2;suppSEMAll=suppMean2;
    for iCon=1:nCon
        suppMean2(iCon) = mean(suppInd(find(ism2(:,iCon)),iCon));
        suppSEM2(iCon) = std(suppInd(find(ism2(:,iCon)),iCon))./sqrt(sum(ism2(:,iCon)));
        suppMeanAll(iCon) = mean(suppInd(:,iCon));
        suppSEMAll(iCon) = std(suppInd(:,iCon))./sqrt(nCellsi);
    end
    suppMeani(i) = suppMeanAll(nCon);
    suppSEMi(i) = suppSEMAll(nCon);
    if sum(choosefig==4)
        figure(4);if i==1;clf;end %figure 4 = suppInd vs con
        subplot(2,4,i)
        errorbar(cons,suppMean2,suppSEM2,'^-');
        hold on
        errorbar(cons,suppMeanAll,suppSEMAll,'kx-');
        hold off
        title({sprintf('Area:%s, GECI:%s',conds(i,1),conds(i,2));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Contrast')
        ylabel('Supp Ind')
        ylim([0 1])
        if i==8;legend('m2 only','all');end
    end
    
    
end

nCells_area = reshape(nCells_area,2,4); %V1 LM AL PM, 1st row f, 2nd row s
if sum(choosefig==5) % model fractions at highest con, f vs s
    modelcounts2 = zeros(2,5);
    modelcounts2(:,1:4) = reshape(modelcountsi,2,4);
    modelcounts2(:,5) = sum(nCells_area.*modelcounts2(:,1:4),2)./sum(nCells_area,2);
    figure(5);clf %figure 5 = proportions of model2 at con=0.8, f vs s
    colormap([1 0 0;0 0 1])
    bar(c_areas,modelcounts2','grouped')
    title('Model fractions at con=0.8')
    xlabel('Area')
    ylabel('Frac. model2 cells')
    legend('f','s','location','se')
end

if sum(choosefig==6 )% prefSize at highest con, f vs s
    prefMean2 = reshape(prefMeani,2,4);
    prefSEM2 = reshape(prefSEMi,2,4);
    figure(6);clf %figure 5 = proportions of model2 at con=0.8, f vs s
    colormap([1 0 0;0 0 1])
    errorbar([1:4]-0.05,prefMean2(1,:),prefSEM2(1,:),'ro')
    hold on
    errorbar([1:4]+0.05,prefMean2(2,:),prefSEM2(2,:),'bo')
    title('Preferred Size at con=0.8')
    xlabel('Area')
    ylabel('Pref. Size')
    legend('f','s','location','se')
    xlim([0 5])
    ylim([0 40])
    set(gca,'XTick',1:4,'XTickLabel',{'V1','LM','AL','PM'})
end

if sum(choosefig==7) % SI at highest con, f vs s
    suppMean2 = reshape(suppMeani,2,4);
    suppSEM2 = reshape(suppSEMi,2,4);
    figure(7);clf %figure 5 = proportions of model2 at con=0.8, f vs s
    colormap([1 0 0;0 0 1])
    errorbar([1:4]-0.05,suppMean2(1,:),suppSEM2(1,:),'ro')
    hold on
    errorbar([1:4]+0.05,suppMean2(2,:),suppSEM2(2,:),'bo')
    title('Suppression Index at con=0.8')
    xlabel('Area')
    ylabel('SI')
    legend('f','s','location','se')
    xlim([0 5])
    ylim([0 1])
    set(gca,'XTick',1:4,'XTickLabel',{'V1','LM','AL','PM'})
end