%% look at ret in all these experiments

%% load csv to define exp

clear all; clc;
%expfile = '\\CRASH.dhe.duke.edu\data\home\kevin\Code\Ai9x_experiment_list.txt';
expfile = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Code\Ai9x_experiment_list.txt';
fID = fopen(expfile);
head = textscan(fID,'%s%s%s%s%s',1,'delimiter',',');
head = vertcat(head{:});
temp = textscan(fID,'%s%s%s%s%s','delimiter',',','HeaderLines',1);
temp = horzcat(temp{:});
expdata = cell2table(temp,'VariableNames',head);
nExp = size(expdata,1);
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);

fprintf(['Receptive-field-size visual-area comparison analysis - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

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
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_sizeTuneData.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'sizeTune', 'sizeMean', 'sizeSEM', 'cellDists')
    if size(sizeTune,2) == 6 % if 6 con, only take middle 4 cons (2-5)
        sizeTune_all = cat(3,sizeTune_all,sizeTune(:,2:5,:));
        sizeMean_all = cat(3,sizeMean_all,sizeMean(:,2:5,:));
        sizeSEM_all = cat(3,sizeSEM_all,sizeSEM(:,2:5,:));
    else
        sizeTune_all = cat(3,sizeTune_all,sizeTune);
        sizeMean_all = cat(3,sizeMean_all,sizeMean);
        sizeSEM_all = cat(3,sizeSEM_all,sizeSEM);
    end
    cellDists_all = [cellDists_all;cellDists];
    nCellsExp(i) = length(cellDists);
    fprintf('done\n')
end

% lbub_fits
fprintf('Loading lbub_fits\n')
lbub_fits_all = [];
goodfit_ind_all = [];
stimAz = zeros(1,nExp); stimEl = stimAz;
cellAz_all = []; cellAz_norm_all = [];
cellEl_all = []; cellEl_norm_all = [];
nCellsExpRet = zeros(1,nExp);
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = 'runs-002'; %expdata.run_str{i};
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind')
    lbub_fits_all = cat(1,lbub_fits_all,lbub_fits);
    nCellsExpRet(i) = size(lbub_fits,1);
    tempinds = sum(nCellsExpRet(1:i-1)) + goodfit_ind; % offset by # cells in previous exps
    goodfit_ind_all = [goodfit_ind_all tempinds];
    
    run_str = expdata.run_str{i};
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_input.mat'] ' not found! Please remove from list\n'])
    end
    load(filename)
    
    stimAz(i) = double(input.gratingAzimuthDeg);
    stimEl(i) = double(input.gratingElevationDeg);
    cellAz_all = cat(1,cellAz_all,lbub_fits(goodfit_ind,4,4));
    cellEl_all = cat(1,cellEl_all,lbub_fits(goodfit_ind,5,4) - stimEl(i));
    cellAz_norm_all = cat(1,cellAz_norm_all,lbub_fits(goodfit_ind,4,4) - stimAz(i));
    cellEl_norm_all = cat(1,cellEl_norm_all,lbub_fits(goodfit_ind,5,4) - stimEl(i));
    fprintf('done\n')
end

% goodfit_ind_size
fprintf('Loading goodfit_ind_size\n')
goodfit_ind_size_all = [];
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = expdata.run_str{i};
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'goodfit_ind_size')
    tempinds = sum(nCellsExp(1:i-1)) + goodfit_ind_size; % offset by # cells in previous exps
    goodfit_ind_size_all = [goodfit_ind_size_all tempinds];
    fprintf('done\n')
end

fprintf('Loading sizeFitResults_SP (true fit at all cons)\n')
sizeFits_all = struct([]); % no cells, 4 cons
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = expdata.run_str{i};
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

fprintf(['\nFinished loading all ' num2str(nExp) ' experiments.\n'])

%% pref size vs azimuth analysis
Az_all = cellAz_all(goodfit_ind_size_all);
Aznorm_all = cellAz_norm_all(goodfit_ind_size_all);
prefSize_all = reshape([sizeFits_all(goodfit_ind_size_all,4).prefSize],size(goodfit_ind_size_all));

figure(1);clf;
subplot(2,2,1)
plot(cellDists_all(goodfit_ind_size_all),prefSize_all,'.')
title('prefSize vs cellDist')
xlabel('Cell Dist (deg)')
ylabel('Pref size (deg)')
hold on
x = cellDists_all(goodfit_ind_size_all);
xsamp = 0.1:0.1:40;
p = polyfit(x,prefSize_all',1)
pSfit = polyval(p,xsamp);
plot(xsamp,pSfit,'-r')
SSresid = sum([prefSize_all' - polyval(p,x)].^2);
SStot = (length(prefSize_all)-1)*var(prefSize_all);
Rsq = 1-SSresid/SStot
text(30, 50,['slope=' num2str(p(1),3)])
text(30, 40,['R^2=' num2str(Rsq,3)])
xlim([0 50])
subplot(2,2,2)
plot(Az_all,prefSize_all,'.')
title('prefSize vs Az (raw)')
xlabel('Pref size (deg)')
ylabel('Azimuth (deg)')
subplot(2,2,3)
plot(abs(Aznorm_all),cellDists_all(goodfit_ind_size_all),'.')
title('cellDist vs |Az (norm)|')
xlabel('Cell Dist (deg)')
ylabel('|Azimuth| (deg)')
subplot(2,2,4)
plot(Aznorm_all,prefSize_all,'.')
title('prefSize vs |Az (norm)|')
xlabel('Pref Size (deg)')
ylabel('|Azimuth| (deg)')
hold on
x = Aznorm_all;
%xsamp = 0.1:0.1:30;
xsamp = -30:0.1:30;
p = polyfit(x,prefSize_all',2)
pSfit = polyval(p,xsamp);
plot(xsamp,pSfit,'-r')
SSresid = sum([prefSize_all' - polyval(p,x)].^2);
SStot = (length(prefSize_all)-1)*var(prefSize_all);
Rsq = 1-SSresid/SStot
text(20, 40,['R^2=' num2str(Rsq,3)])

%% RF size vs RF-stim distance
% for each area plot all cells:
% histogram of RF-stim distance (with vertical line at cutoff)
% histogram of RF size
% plot of RF-size vs RF-stim distance (with maesure of correlation)

% extract RF size (lbub_fits_all -> only goodfit_ind -> only goodfit_ind_size)
%sigmax = lbub_fits_all(goodfit_ind_all,2,4);
%sigmay = lbub_fits_all(goodfit_ind_all,3,4);
sigmax = lbub_fits_all(:,2,4);
sigmay = lbub_fits_all(:,3,4);
RFsize_all = 2*sqrt(2*log(2))*geo_mean([sigmax sigmay],2);

fprintf('Examine cells from each area:\n')
areas = ["V1","LM","AL","PM"];

expInd = [];
areaInd = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExpRet(i))];
    switch cell2mat(expdata.area(i))
        case 'V1'
            areaInd = [areaInd repmat(1,1,nCellsExpRet(i))];
        case 'LM'
            areaInd = [areaInd repmat(2,1,nCellsExpRet(i))];
        case 'AL'
            areaInd = [areaInd repmat(3,1,nCellsExpRet(i))];
        case 'PM'
            areaInd = [areaInd repmat(4,1,nCellsExpRet(i))];
    end
end

x=[];y_geo=[];

for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good size fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_all);%(goodfit_ind_size_all));
    % or comment out and just take all cells with good ret
    %ind = find(ismember(expInd,expIndi));
    
    % no cutoff for now
    % cutoff by cellDist
    % try looking with different cutoffs
    switch areas(i)
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
            RFsize_V1 = RFsize_all(ind);
        case 'LM'
            cutoff = 15; %alm cutoff at 15
            RFsize_LM = RFsize_all(ind);
        case 'AL'
            cutoff = 15; %alm cutoff at 15
            RFsize_AL = RFsize_all(ind);
        case 'PM'
            cutoff = 20; %pm cutoff at 20
            RFsize_PM = RFsize_all(ind);
    end
    %ind = intersect(ind,find(cellDists_all<cutoff));
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    cellDists = cellDists_all(ismember(goodfit_ind_all,ind));
    RFsize = RFsize_all(ind);
    
%     sizeFits = sizeFits_all(ind,:); %cell,con
%     ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
%     ism2 = reshape([sizeFits.Ftest],size(sizeFits));
    
    cons_c = categorical({'0.1' '0.2' '0.4' '0.8'});
    
    RFrads_geo{i} = RFsize;
    x = [x; i*ones(size(RFsize))];
    y_geo = [y_geo; RFsize];
    
    % rf-stim distance
    figure(2);if i==1;clf;end
    subplot(2,2,i)
    histogram(cellDists,[0:1:40])
    hold on
    line([cutoff cutoff],[0 1000],'color','red','LineStyle','--')
    title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    xlabel('RF-stim distance')
    ylabel('num cells')
    ylim([0 nCellsi/6])
    text(2,nCellsi/7,['n_{cut}=' num2str(sum(cellDists<cutoff))])
    %if i==4;legend('m1','m2');end
    
    % rf size
    figure(3);if i==1;clf;end
    subplot(2,2,i)
    histogram(RFsize,[0:1:40])
    title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    xlabel('RF size (FWHM, deg)')
    ylabel('num cells')
    
    % rf size vs rf-stim dist
    figure(4);if i==1;clf;end
    subplot(2,2,i)
    plot(cellDists,RFsize,'.')
    title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    %cc = corrcoef(cellDists,RFsize);
    %text(1,2,['R^2=' num2str(cc(2))])
    xlabel('RF-stim dist')
    ylabel('RF size (FWHM, deg)')
    
    %     ind = intersect(ind,find(cellDists_all<cutoff));
    %     cellAz_norm = cellAz_norm_all(ind);
    %     cellEl_norm = cellEl_norm_all(ind);
    % plot locations (normalized)
    %     figure(5);if i==1;clf;end
    %     subplot(2,2,i)
    %     plot(cellAz_norm,cellEl_norm,'o')
    %     title({sprintf('Area:%s',areas(i));['(n=' num2str(length(ind)) ', n_{exp}=' num2str(nExpi) ')']})
    %     xlabel('Az rel to stim (deg)')
    %     ylabel('El rel to stim (deg)')
end
figure(1);clf;
boxplot(y_geo,x)
hold on
y_mean = [mean(RFrads_geo{1}) mean(RFrads_geo{2}) mean(RFrads_geo{3}) mean(RFrads_geo{4})];
y_std = [std(RFrads_geo{1}) std(RFrads_geo{2}) std(RFrads_geo{3}) std(RFrads_geo{4})];
plot(1:4,y_mean,'x')
hold off
set(gca,'XTick',1:4,'XTickLabel',{['V1 (n=' num2str(length(RFrads_geo{1})) ')'],['LM (n=' num2str(length(RFrads_geo{2})) ')'],['AL (n=' num2str(length(RFrads_geo{3})) ')'],['PM (n=' num2str(length(RFrads_geo{4})) ')']})
set(gcf,'Color','w')
%set(gca,'XAxisLocation','bottom','YAxisLocation','left','TickDir','out')
set(gca,'box','off','TickDir','out')
%title('RF radius')
xlabel('Area')
ylabel('FWHM (deg)')
hold on
plot([3 4],2*[28 28],'-k', 'LineWidth',2)
plot([3.4 3.5 3.6],2*[29 29 29],'*k')
plot([1 2],2*[28 28],'-k', 'LineWidth',2)
plot([1.4 1.5 1.6],2*[29 29 29],'*k')
plot([1 3],2*[31 31],'-k', 'LineWidth',2)
plot([1.9 2.0 2.1],2*[32 32 32],'*k')
plot([2 4],2*[34 34],'-k', 'LineWidth',2)
plot([2.9 3.0 3.1],2*[35 35 35],'*k')
plot([1 4],2*[37 37],'-k', 'LineWidth',2)
plot([2.4 2.5 2.6],2*[38 38 38],'*k')
ylim([0 2*40])
%%
[p,~,stat_geo] = anova1(y_geo,x,'Display','off');
[results, means] = multcompare(stat_geo,'CType','hsd');
[p,tbl,stat_geo]=kruskalwallis(y_geo,x)
[results, ~] = multcompare(stat_geo,'CType','hsd')
%mean and sd, no boxplot
figure(1);clf;
y_mean = [mean(RFrads_geo{1}) mean(RFrads_geo{2}) mean(RFrads_geo{3}) mean(RFrads_geo{4})];
y_std = [std(RFrads_geo{1}) std(RFrads_geo{2}) std(RFrads_geo{3}) std(RFrads_geo{4})];
set(gcf,'Color','w')
errorbar(1:4,means(:,1),means(:,2),'.k');
set(gca,'box','off','TickDir','out')
set(gca,'XTick',1:4,'XTickLabel',areas,'TickLength',[0.015 0.015])
ylabel('HWHM (deg)')
xlim([0.5 4.5])
ylim([0 20])
% 
% hold on
% plot([3 4],[28 28],'-k', 'LineWidth',2)
% plot([3.4 3.5 3.6],[29 29 29],'*k')
% plot([1 2],[28 28],'-k', 'LineWidth',2)
% plot([1.4 1.5 1.6],[29 29 29],'*k')
% plot([1 3],[31 31],'-k', 'LineWidth',2)
% plot([1.9 2.0 2.1],[32 32 32],'*k')
% plot([2 4],[34 34],'-k', 'LineWidth',2)
% plot([2.9 3.0 3.1],[35 35 35],'*k')
% plot([1 4],[37 37],'-k', 'LineWidth',2)
% plot([2.4 2.5 2.6],[38 38 38],'*k')
% ylim([0 40])

%% try 4 new plots
% average(current), dotplot, CDF, violin plot, KDE

figure(2);clf;set(gcf,'Color','w');
% % 1. average
% subplot(1,2,1)
% errorbar(1:4,means(:,1),y_std,'ok');
% set(gca,'box','off','TickDir','out')
% set(gca,'XTick',1:4,'XTickLabel',areas,'TickLength',[0.015 0.015])
% ylabel('HWHM (deg)')
% xlim([0.5 4.5])
% ylim([0 20])
%5 sideways violin
for i = 1:length(areas)
    [fi xi] = ksdensity(RFrads_geo{i});
    fnorm(:,i) = fi/max(fi)*0.3;
    xinorm(:,i) = xi;
end
ax = subplot(1,1,1);
colors = get(ax,'ColorOrder');
for i=1:length(areas)
    hold on
    h5(i)=fill([xinorm(:,i);flipud(xinorm(:,i))],[fnorm(:,i)+(5-i);flipud((5-i)-fnorm(:,i))],[1 1 1],'EdgeColor','k');
    p(1)=plot([means(i,1) y_mean(i)],[interp1(xinorm(:,i),fnorm(:,i)+(5-i),means(i,1)), interp1(flipud(xinorm(:,i)),flipud((5-i)-fnorm(:,i)),means(i,1)) ],'k','LineWidth',2);
    h5(i).FaceColor = colors(i,:);
end
axis([0 40 0.5 length(areas)+0.5]);
legend off
ax.YTick = [1:4];
ax.YTickLabel = fliplr(areas);
set(gca,'box','off','TickDir','out')
ylabel('Area')
xlabel('RF size (deg)')
%filename = 'N:\home\kevin\ppts\_paper figs\plots\RF_testplots.pdf';
%set(gcf,'PaperSize',[12 3])
%print(filename,'-dpdf','-fillpage')

%% eccentricity analysis
% for each experiment, take mean Az and mean El of cells
% also take stim Az and stim El and thirdly normalize by cell RFs by stim
% plots:
% distribution of eccentricity by area (overlaid histograms or violin?)
% 4 plots for each area: scatter of individual exp average RF center, and line to stim, and mean
% ?
areas = ["V1","LM","AL","PM"];

expInd = []; areaInd =[];
expArea = zeros(1,nExp);
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExpRet(i))];
    switch cell2mat(expdata.area(i))
        case 'V1'
            areaInd = [areaInd repmat(1,1,nCellsExpRet(i))];
        case 'LM'
            areaInd = [areaInd repmat(2,1,nCellsExpRet(i))];
        case 'AL'
            areaInd = [areaInd repmat(3,1,nCellsExpRet(i))];
        case 'PM'
            areaInd = [areaInd repmat(4,1,nCellsExpRet(i))];
    end
    expArea(i) = find(strcmp(expdata.area{i},areas));
end

expAz = zeros(1,nExp);
expEl = expAz;
expAz_norm = expAz;
expEl_norm = expAz;
for iExp = 1:nExp
    ind = find(ismember(goodfit_ind_all,find(expInd==iExp)));
    expAz(iExp) = mean(cellAz_all(ind));
    expEl(iExp) = mean(cellEl_all(ind));
    expAz_norm(iExp) = mean(cellAz_norm_all(ind));
    expEl_norm(iExp) = mean(cellEl_norm_all(ind));
end
ecc = sqrt(expAz.^2 + expEl.^2);
ecc_norm = sqrt(expAz_norm.^2 + expEl_norm.^2);
ecc_stim = sqrt(stimAz.^2 + stimEl.^2);

figure(1);clf
for i=1:length(areas)
    ecc_cell{i} = ecc(expArea==i)
    subplot(1,3,1)
    hold on
    histogram(ecc(expArea==i))
    hold off
    ecc_norm_cell{i} = ecc(expArea==i)
    subplot(1,3,2)
    hold on
    histogram(ecc_norm(expArea==i))
    hold off
    ecc_stim_cell{i} = ecc(expArea==i)
    subplot(1,3,3)
    hold on
    histogram(ecc_stim(expArea==i))
    hold off
end

figure(2);clf
for iF=1:3
    switch iF
        case 1
            data = ecc;
            data_cell = ecc_cell;
            xlab = 'Mean Exp. Eccentricity (deg)';
        case 2
            data = ecc_norm;
            data_cell = ecc_norm_cell;
            xlab = 'Exp Eccent. norm to Stim (deg)';
        case 3
            data = ecc_stim;
            data_cell = ecc_stim_cell;
            xlab = 'Stimulus Eccentricity (deg)';
    end
    [p,~,stat_geo] = anova1(data,expArea,'off');
    [results, means] = multcompare(stat_geo,'CType','hsd','Display','off')

    for i = 1:length(areas)
        [fi xi] = ksdensity(data_cell{i});
        fnorm(:,i) = fi/max(fi)*0.3;
        xinorm(:,i) = xi;
    end
    figure(2);
    ax = subplot(1,3,iF);
    colors = get(ax,'ColorOrder');
    for i=1:length(areas)
        hold on
        h5(i)=fill([xinorm(:,i);flipud(xinorm(:,i))],[fnorm(:,i)+(5-i);flipud((5-i)-fnorm(:,i))],[1 1 1],'EdgeColor','k');
        p(1)=plot([means(i,1) means(i,1)],[interp1(xinorm(:,i),fnorm(:,i)+(5-i),means(i,1)), interp1(flipud(xinorm(:,i)),flipud((5-i)-fnorm(:,i)),means(i,1)) ],'k','LineWidth',2);
        h5(i).FaceColor = colors(i,:);
    end
    axis([-5 25 0.5 length(areas)+0.5]);
    legend off
    ax.YTick = [1:4];
    ax.YTickLabel = fliplr(areas);
    set(gca,'box','off','TickDir','out')
    ylabel('Area')
    xlabel(xlab)
end

%% scatter of exp eccentricities by area
figure(3);clf;
linkaxes([subplot(2,2,1) subplot(2,2,2) subplot(2,2,3) subplot(2,2,4)],'xy')
xlim([-10 20])
ylim([-10 10])
for iExp=1:nExp
    subplot(2,2,expArea(iExp))
    hold on
    plot([expAz(iExp) stimAz(iExp)],[expEl(iExp) stimEl(iExp)],'-','Color',colors(expArea(iExp),:))
    plot(expAz(iExp),expEl(iExp),'o','Color',colors(expArea(iExp),:))
    plot(stimAz(iExp),stimEl(iExp),'s','Color',colors(expArea(iExp),:))
    hold off
end
for i=1:length(areas)
    subplot(2,2,i)
    title(areas(i))
    hold on
    plot(mean(expAz(expArea==i)),mean(expEl(expArea==i)),'*','Color',colors(i,:))
end

% same plot but normalized
figure(4);clf;
linkaxes([subplot(2,2,1) subplot(2,2,2) subplot(2,2,3) subplot(2,2,4)],'xy')
xlim([-20 20])
ylim([-20 20])
for iExp=1:nExp
    subplot(2,2,expArea(iExp))
    hold on
    plot(expAz_norm(iExp),expEl_norm(iExp),'o','Color',colors(expArea(iExp),:))
    hold off
end
for i=1:length(areas)
    subplot(2,2,i)
    title(areas(i))
    hold on
    plot(mean(expAz_norm(expArea==i)),mean(expEl_norm(expArea==i)),'*','Color',colors(i,:))
    fprintf('\nMean eccentricity in %s: %.3f\n',areas(i),mean(ecc(expArea==i)))
    fprintf('Mean norm eccentricity in %s: %.3f\n',areas(i),mean(ecc_norm(expArea==i)))
    fprintf('Mean stim eccentricity in %s: %.3f\n',areas(i),mean(ecc_stim(expArea==i)))
end

%% supp fig
% A: RF centers across areas
% do we need a scatter of the positions? => 4 subplots or one overlay (too cluttered)
% or we can do violins/histograms
% B: RF center vs RF Size
% C: RF center vs SI

figure(17);clf
sp = [1 4 7 10];
linkaxes([subplot(4,3,1) subplot(4,3,4) subplot(4,3,7) subplot(4,3,10)],'xy')
xlim([-20 20])
ylim([-20 20])
for iExp=1:nExp
    subplot(4,3,sp(expArea(iExp)))
    hold on
    plot(expAz_norm(iExp),expEl_norm(iExp),'.','Color',colors(expArea(iExp),:))
    hold off
end
for i=1:length(areas)
    subplot(4,3,sp(i))
    title(areas(i))
    hold on
    plot(mean(expAz_norm(expArea==i)),mean(expEl_norm(expArea==i)),'o','Color',colors(i,:))
    xlabel('Norm Az (deg)')
    ylabel('Norm El (deg)')
    axis square
end

% B: eccentricity vs RF size
% by experiments or by individual cells? assuming cells
% for now show all cells in one plot, just color by area
% too crowded, will need to do 4 subplots to see each area clearly
sp = [2 5 8 11];
cellEcc_all = sqrt(cellAz_all.^2 + cellEl_all.^2);
cellEcc_norm_all = sqrt(cellAz_norm_all.^2 + cellEl_norm_all.^2);
%cellEcc_norm_all = cellEcc_all;
fprintf('\nEcc vs RF size slope analysis:\n')
for i = 1:length(areas)
    subplot(4,3,sp(i))
    ind = find(ismember(goodfit_ind_all,find(ismember(expInd,find(expArea==i)))));
    ind = intersect(ind, goodfit_ind_size_all);
    hold on
    plot(cellEcc_norm_all(ind),RFsize_all(goodfit_ind_all(ind)),'.','Color',colors(i,:));
    title([char(areas(i)) ' (n=' num2str(length(ind)) ')']) 
    xlabel('RF Eccentricity (deg)')
    ylabel('RF size (FWHM, deg)')
    x = cellEcc_norm_all(ind);
    y = RFsize_all(goodfit_ind_all(ind));
    if i == 1
        xsamp = 0.1:0.1:20;
    else
        xsamp = 0.1:0.1:40;
    end
    [p, S] = polyfit(x,y,1);
    ci = polyparci(p,S);
    pSfit = polyval(p,xsamp);
    plot(xsamp,pSfit,'-r')
    resid = y - polyval(p,x);
    df = length(y)-2;
    SSresid = sum(resid.^2);
    SStot = (length(y)-1)*var(y);
    Rsq = 1-SSresid/SStot;
    %text(5, 25,['slope=' num2str(p(1),3)])
    %text(5, 20,['R^2=' num2str(Rsq,3)])
    sb1 = sqrt(SSresid./df)./sqrt(sum((x-mean(x)).^2));
    t = p(1)/sb1;
    pval = 1-tcdf(t,df);
    fprintf('Area %s: slope=%0.3f, p=%f, ci = (%0.3f, %0.3f)\n',areas(i),p(1),pval,ci(1,1),ci(2,1))
end

% C: eccentricity vs SI
% by experiments or by individual cells? assuming cells
% for now show all cells in one plot, just color by area
sp = [3 6 9 12];
nCon = 4;
fprintf('\nEcc vs SI slope analysis:\n')
for i = 1:length(areas)
    subplot(4,3,sp(i))
    ind = find(ismember(goodfit_ind_all,find(ismember(expInd,find(expArea==i)))));
    ind = intersect(ind, goodfit_ind_size_all);
    SI = [sizeFits_all(ind,nCon).suppInd];
    hold on
    plot(cellEcc_norm_all(ind),SI,'.','Color',colors(i,:));
    xlabel('RF Eccentricity (deg)')
    ylabel('Supp Index')
    x = cellEcc_norm_all(ind);
    x = x(SI~=0);
    y = SI(SI~=0)'; %(SI~=0)
    if i == 1
        xsamp = 0.1:0.1:20;
    else
        xsamp = 0.1:0.1:40;
    end
    [p, S] = polyfit(x,y,1);
    ci = polyparci(p,S);
    pSfit = polyval(p,xsamp);
    plot(xsamp,pSfit,'-r') 
    resid = y - polyval(p,x);
    df = length(y)-2;
    SSresid = sum(resid.^2);
    SStot = (length(y)-1)*var(y);
    Rsq = 1-SSresid/SStot;
    %text(5, 1.4,['slope=' num2str(p(1),3)])
    %text(5, 1.1,['R^2=' num2str(Rsq,3)])
    sb1 = sqrt(SSresid./df)./sqrt(sum((x-mean(x)).^2));
    t = p(1)/sb1;
    pval = tcdf(t,df);
    fprintf('Area %s: slope=%0.3f, p=%f, ci = (%0.3f, %0.3f)\n',areas(i),p(1),pval,ci(1,1),ci(2,1))
end
ind(cellEcc_norm_all(ind)>10) = [];
mean(RFsize_all(goodfit_ind_all(ind)))
SI = [sizeFits_all(ind,nCon).suppInd];
mean(SI)
mean(SI(SI~=0))

%%
function CI = polyparci(PolyPrms,PolyS,alpha)
% POLYPARCI takes PolyPrms, the 1xN vector of parameter estimates from polyfit,
%   PolyS, the structure returned by polyfit, and alpha, and returns the
%   ?alpha? confidence intervals for the parameters in CI.  The default
%   for alpha is 0.95.  
% 
%   RETURNS: CI (2xN matrix of confidence intervals whose columns match the 
%            1xN vector of coefficients returned by polyfit)
% 
% CALLING POLYPARCI: 
% [p,S,mu] = polyfit(x,y,n);    % Fit the data (polyparci does not use ?mu?)
% ci = polyparci(p,S);          % NOTE: specifying a value for ?alpha? is 
%                                 optional
% 
% Star Strider ? 2012 07 05;    
% Update         2014 02 10 (Compatible with 2013b changes. Corrected typo. 
%                   Inverted ?CI? matrix to 2xN be compatible with 1xN row 
%                   vector of coefficients ?polyfit? produces.)
% 
% Check for a specified value of alpha: 
if nargin < 3
    alpha = 0.95;
end
% Check for out-of-range values for alpha and substitute if necessary: 
if alpha < 1.0E-010
    alpha = 1.0E-010;
elseif alpha > (1 - 1.0E-010)
    alpha = 1 - 1.0E-010;
end
% Calculate the covariance matrix of the parameters (COVB) and the standard
%   errors (SE) from the ?S? structure.  (See the ?polyfit? documentation.)   
COVB = (PolyS.R'*PolyS.R)\eye(size(PolyS.R)) * PolyS.normr^2/PolyS.df;
SE = sqrt(diag(COVB));                              % Standard Errors
[PrmSizR,PrmSizC] = size(PolyPrms);                 % Convert parameter vector to column if necessary
if PrmSizR < PrmSizC
    PolyPrms = PolyPrms';
end
% Calculate t-statistic and inverse t-statistic:
%       ?tstat? is an implicit function that calculates the probability
%       from the cumulative t-distribution ?t_cdf? for various values of
%       ?tval? supplied to it by the ?fzero? call until it converges to the
%       required value of alpha.  The ?t_cdf? function calculates the
%       cumulative t-distribution.  The two lines that calculate ?tstat?
%       and ?T? together calculate the inverse t-distribution, returning it
%       as ?T?.
tstat = @(tval) (alpha - t_cdf(tval,PolyS.df) );    % Function to calculate t-statistic for p = ?alpha? and v = ?PolyS.df?
[T,fval] = fzero(tstat, 1);                         % Calculate t-statistic for p = ?alpha? and v = ?PolyS.df?
T = abs(T);                                         % Critical +ve value from t-distribution
ts = T * [-1  1];                                   % Create vector of ? t-statistic values
CI  = bsxfun(@plus,bsxfun(@times,SE,ts),PolyPrms)'; % Confidence Intervals Matrix (2xN)
% CALCULATE THE CUMULATIVE T-DISTRIBUTION: 
    function PT = t_cdf(t,v)
        % t_cdf(t,v) calculates the cumulative t-distribution probability
        %   given the t-statistic ?t? and degrees-of-freedom ?v?.  The
        %   routine to calculate the inverse t-distribution uses this,
        %   ?tstat?, and the ?fzero? call.  Compared to the Statistics
        %   Toolbox function ?tcdf? and ?tinv?, ?t_cdf? and ?T? have
        %   relative errors of about 1E-12.
        
        IBx = v./(t.^2 + v);                      % ?x? for IxZW ? NOTE: function of ?t?-statistic
        IBZ = v/2;                                % ?Z? for IxZW
        IBW = 0.5;                                % ?W? for IxZW
        Ixzw = betainc(IBx, IBZ, IBW);            % Incomplete beta function
        PT = 1-0.5*Ixzw;                          % Cumulative t-distribution        
    end
end
