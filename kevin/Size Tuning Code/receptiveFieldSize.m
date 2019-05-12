%% look at ret in all these experiments

%% load csv to define exp

clear all; clc;
%expfile = '\\CRASH.dhe.duke.edu\data\home\kevin\Code\Ai9x_experiment_list.txt';
expfile = 'C:\Users\kevin\Documents\Repositories\ImagingCode-Glickfeld-Hull\kevin\Size Tuning Code\Ai9x_experiment_list.txt';
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
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
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
cellAz_norm_all = []; cellAz_all = [];
cellEl_norm_all = [];
nCellsExpRet = zeros(1,nExp);
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = 'runs-002'; %expdata.run_str{i};
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind')
    lbub_fits_all = cat(1,lbub_fits_all,lbub_fits);
    nCellsExpRet(i) = size(lbub_fits,1);
    tempinds = sum(nCellsExpRet(1:i-1)) + goodfit_ind; % offset by # cells in previous exps
    goodfit_ind_all = [goodfit_ind_all tempinds];
    
    run_str = expdata.run_str{i};
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_input.mat'] ' not found! Please remove from list\n'])
    end
    load(filename)
    
    stimAz = double(input.gratingAzimuthDeg);
    cellAz_all = cat(1,cellAz_all,lbub_fits(goodfit_ind,4,4));
    cellAz_norm_all = cat(1,cellAz_norm_all,lbub_fits(goodfit_ind,4,4) - stimAz);
    %cellEl_norm_all = cat(1,cellEl_norm_all,lbub_fits(goodfit_ind,5,4) - stimPosExp(i,1));
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
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
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
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_SP.mat']);
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
subplot(2,2,2)
plot(Az_all,prefSize_all,'.')
title('prefSize vs Az (raw)')
subplot(2,2,3)
plot(abs(Aznorm_all),cellDists_all(goodfit_ind_size_all),'.')
title('cellDist vs |Az (norm)|')
subplot(2,2,4)
x = Aznorm_all;
%xsamp = 0.1:0.1:30;
xsamp = -30:0.1:30;
plot(x,prefSize_all,'.')
title('prefSize vs |Az (norm)|')
hold on
p = polyfit(x,prefSize_all',2)
pSfit = polyval(p,xsamp);
plot(xsamp,pSfit,'-r')
SSresid = sum([prefSize_all' - polyval(p,x)].^2);
SStot = (length(prefSize_all)-1)*var(prefSize_all);
Rsq = 1-SSresid/SStot
text(20, 40,['R^2=' num2str(Rsq)])
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
RFsize_all = sqrt(2*log(2))*geomean([sigmax sigmay],2);

fprintf('Examine cells from each area:\n')
areas = ["V1","LM","AL","PM"];

expInd = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExpRet(i))];
end

x=[];y_geo=[];

for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good size fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_all(goodfit_ind_size_all));
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
    %cellDists = cellDists_all(ind);
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
    %histogram(cellDists,[0:1:40])
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
    xlabel('RF size')
    ylabel('num cells')
    
    % rf size vs rf-stim dist
    figure(4);if i==1;clf;end
    subplot(2,2,i)
    %plot(cellDists,RFsize,'.')
    title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    %cc = corrcoef(cellDists,RFsize);
    %text(1,2,['R^2=' num2str(cc(2))])
    xlabel('RF-stim dist')
    ylabel('RF size')
    
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
ylabel('HWHM (deg)')
hold on
plot([3 4],[28 28],'-k', 'LineWidth',2)
plot([3.4 3.5 3.6],[29 29 29],'*k')
plot([1 2],[28 28],'-k', 'LineWidth',2)
plot([1.4 1.5 1.6],[29 29 29],'*k')
plot([1 3],[31 31],'-k', 'LineWidth',2)
plot([1.9 2.0 2.1],[32 32 32],'*k')
plot([2 4],[34 34],'-k', 'LineWidth',2)
plot([2.9 3.0 3.1],[35 35 35],'*k')
plot([1 4],[37 37],'-k', 'LineWidth',2)
plot([2.4 2.5 2.6],[38 38 38],'*k')
ylim([0 40])
%%
[p,~,stat_geo] = anova1(y_geo,x)
[results, means] = multcompare(stat_geo,'CType','hsd')
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
% 1. average
subplot(1,2,1)
errorbar(1:4,means(:,1),y_std,'ok');
set(gca,'box','off','TickDir','out')
set(gca,'XTick',1:4,'XTickLabel',areas,'TickLength',[0.015 0.015])
ylabel('HWHM (deg)')
xlim([0.5 4.5])
ylim([0 20])
%5 sideways violin
for i = 1:length(areas)
    [fi xi] = ksdensity(RFrads_geo{i});
    fnorm(:,i) = fi/max(fi)*0.3;
    xinorm(:,i) = xi;
end
ax = subplot(1,2,2);
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