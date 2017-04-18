ds = '_V1';
%% load data
close all
eval(['awFSAVdatasets' ds])
rc = behavConstsAV;
titleStr = ds;
if strcmp(titleStr, '')
    titleStr = 'V1_100ms';
else
    titleStr = titleStr(2:end);
end

if strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' ds];
else
    dataGroup = [];
end
str = unique({expt.SubNum});
mouse_str = ['i' strjoin(str,'_i')];
load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' ds '.mat']));
load(fullfile(rc.caOutputDir,dataGroup,[titleStr '_' mouse_str '_modCells.mat']));

fnout = fullfile(rc.caOutputDir,dataGroup,'tarAlgin_');
%% 
for C = 1:7
    if C ~= 7
        doCellInd = 1;
    else
        doCellInd = 0;
    end
        
titleStr = ds;
if strcmp(titleStr, '')
    titleStr = 'V1_100ms';
else
    titleStr = titleStr(2:end);
end
if doCellInd
    cell_ind_name = {'0'; '45'; '90'; '135'; 'tar inc'; 'base resp'};
    titleStr = [titleStr '-' cell_ind_name{C}];
    C_ind_used = {2; 3; 4; 5; 10; [8 9 12]};
    C_ind = C_ind_used{C};
end
%% experiment parameters
vis = 1; % visual trials
cAlign = 3; % invalid trial data
minTrLMs = 2000; % cutoff time for short and long trials
oneS_fr = expt(1).frame_rate;
% outcomes
hit_val = 1;
miss_val = 2;
hit_inv = 3;
miss_inv = 4;

% catch data size
ncatch = sum(nCatchExpt(expt,mouse) > 0);

% analysis windows
pre_win = mouse(1).expt(2).win(1).frames;
trans_win = mouse(1).expt(2).win(2).frames;
pre_event_frames = mouse(1).expt(2).pre_event_frames;


%% get just the data needed
nMice = size(mouse,2);

invTars = cell(1,ncatch); % set of matching targets used

% responses for each direction
hitVal = cell(1,ncatch); 
missVal = cell(1,ncatch);
hitInv = cell(1,ncatch);
missInv = cell(1,ncatch);
allVal = cell(1,ncatch); 
allInv = cell(1,ncatch); 

hitVal_trL = cell(1,ncatch); 
missVal_trL = cell(1,ncatch); 
hitInv_trL = cell(1,ncatch); 
missInv_trL = cell(1,ncatch);
allVal_trL = cell(1,ncatch); 
allInv_trL = cell(1,ncatch); 
cell_ind = cell(1,ncatch); 

auROC_expt = cell(1,ncatch); 
auROC_tar_expt = cell(1,ncatch); 
rst_expt = cell(1,ncatch); 
rst_tar_expt = cell(1,ncatch); 

for im = 1:nMice
    if im == 1
        i = 1;
    end
    for iexp = 1:size(mouse(im).expt,2)
       if mouse(im).expt(iexp).info.isCatch;
           %exp info
           cycTimeMs = mouse(im).expt(iexp).info.cyc_time_ms;
           invTars{i} = round(mouse(im).expt(iexp).info.cDirs);
           
           %temp data
           d = mouse(im).expt(iexp).align(cAlign).av(vis);
           
           % cell type indexes
           if doCellInd
               % cell index
               c = mouse(im).expt(iexp).cells;
               cell_ind = c(C_ind).ind;
               %data to analyze
               hitVal{i} = cellfun(@(x) x(:,cell_ind,:), d.outcome(hit_val).stimResp, 'unif',0);
               missVal{i} = cellfun(@(x) x(:,cell_ind,:), d.outcome(miss_val).stimResp, 'unif',0);
               hitInv{i} = cellfun(@(x) x(:,cell_ind,:), d.outcome(hit_inv).stimResp, 'unif',0);
               missInv{i} = cellfun(@(x) x(:,cell_ind,:), d.outcome(miss_inv).stimResp, 'unif',0);
           else
               %data to analyze
               hitVal{i} = d.outcome(hit_val).stimResp;
               missVal{i} = d.outcome(miss_val).stimResp;
               hitInv{i} = d.outcome(hit_inv).stimResp;
               missInv{i} = d.outcome(miss_inv).stimResp;
           end
           
           % auROC cell groups
           if ~doCellInd
               mi = msModCells(im).expt(iexp).av(vis).outcome(1).mi(1).comp(1).auc(3);
               auROC_expt{i} = mi.value;
               rst_expt{i} = mi.rst;
               
               tar_ind = ismember(mi.dirs,invTars{i});
               auROC_tar_expt{i} = mi.value_dirs(tar_ind);
               rst_tar_expt{i} = mi.rst_dirs_logical(tar_ind);
           end
               
           
           hitVal_trL{i} = cellfun(@(x) x.*cycTimeMs, d.outcome(hit_val).tcyc, 'unif',0);
           missVal_trL{i} = cellfun(@(x) x.*cycTimeMs, d.outcome(miss_val).tcyc, 'unif',0);
           hitInv_trL{i} = cellfun(@(x) x.*cycTimeMs, d.outcome(hit_inv).tcyc, 'unif',0);
           missInv_trL{i} = cellfun(@(x) x.*cycTimeMs, d.outcome(miss_inv).tcyc, 'unif',0);
           
           for itar = 1:length(invTars{i})
               allVal{i}{itar} = cat(3,hitVal{i}{itar},missVal{i}{itar});
               allInv{i}{itar} = cat(3,hitInv{i}{itar},missInv{i}{itar});
               allVal_trL{i}{itar} = cat(2,hitVal_trL{i}{itar},missVal_trL{i}{itar});
               allInv_trL{i}{itar} = cat(2,hitInv_trL{i}{itar},missInv_trL{i}{itar});
           end
           
           i = i+1;
       end        
    end
end

%% match targets across valid and invalid trials for short trials
    hitVal_short = cell(1,ncatch);
    hitInv_short = cell(1,ncatch);
    missVal_short = cell(1,ncatch);
    missInv_short = cell(1,ncatch);
    allVal_short = cell(1,ncatch);
    allInv_short = cell(1,ncatch);
    
    hitVal_long = cell(1,ncatch);
    hitInv_long = cell(1,ncatch);
    missVal_long = cell(1,ncatch);
    missInv_long = cell(1,ncatch);
    allVal_long = cell(1,ncatch);
    allInv_long = cell(1,ncatch);
    
    hitVal_all = cell(1,ncatch);
    hitInv_all = cell(1,ncatch);
    missVal_all = cell(1,ncatch);
    missInv_all = cell(1,ncatch);
    allVal_all = cell(1,ncatch);
    allInv_all = cell(1,ncatch);
    
    hitAuROC_all = [];
    hitRst_all = [];
    hitAuROC_long = [];
    hitRst_long = [];
    hitAuROC_short = [];
    hitRst_short = [];
    
    missAuROC_all = [];
    missRst_all = [];
    missAuROC_long = [];
    missRst_long = [];
    missAuROC_short = [];
    missRst_short = [];
    
    allAuROC_all = [];
    allRst_all = [];
    allAuROC_long = [];
    allRst_long = [];
    allAuROC_short = [];
    allRst_short = [];
    
    ncells_expt = cell(1,18);
    [ncells_expt{:}] = deal(zeros(ncatch,1));
    ntrials_expt = cell(1,18);
    [ntrials_expt{:}] = deal(zeros(ncatch,1));
    cond_str = {'hit_short_val';'hit_long_val';'hit_all_val';'miss_short_val';'miss_long_val';'miss_all_val';'all_short_val';'all_long_val';'all_all_val';...
                'hit_short_inv';'hit_long_inv';'hit_all_inv';'miss_short_inv';'miss_long_inv';'miss_all_inv';'all_short_inv';'all_long_inv';'all_all_inv'};

for i = 1:ncatch
%     get short trials  
    [hitVal_short_ind, hitInv_short_ind] = matchTrialLengthIndex(hitVal_trL{i}, hitInv_trL{i}, invTars{i}, minTrLMs,'short');
    [missVal_short_ind, missInv_short_ind] = matchTrialLengthIndex(missVal_trL{i}, missInv_trL{i}, invTars{i}, minTrLMs,'short');
    [allVal_short_ind, allInv_short_ind] = matchTrialLengthIndex(allVal_trL{i}, allInv_trL{i}, invTars{i}, minTrLMs,'short');
      
%     get long trials  
    [hitVal_long_ind, hitInv_long_ind] = matchTrialLengthIndex(hitVal_trL{i}, hitInv_trL{i}, invTars{i}, minTrLMs,'long');
    [missVal_long_ind, missInv_long_ind] = matchTrialLengthIndex(missVal_trL{i}, missInv_trL{i}, invTars{i}, minTrLMs,'long');
    [allVal_long_ind, allInv_long_ind] = matchTrialLengthIndex(allVal_trL{i}, allInv_trL{i}, invTars{i}, minTrLMs,'long');
    
%     all trials
    [hitVal_all_ind, hitInv_all_ind] =  matchTrialLengthIndex(hitVal_trL{i}, hitInv_trL{i}, invTars{i}, minTrLMs,'all');    
    [missVal_all_ind, missInv_all_ind] = matchTrialLengthIndex(missVal_trL{i}, missInv_trL{i}, invTars{i}, minTrLMs,'all');
    [allVal_all_ind, allInv_all_ind] = matchTrialLengthIndex(allVal_trL{i}, allInv_trL{i}, invTars{i}, minTrLMs,'all');
    
%     auROC
    if ~doCellInd
        [hitAuROC_all, hitRst_all] = catIfIndexNotEmpty(hitAuROC_all,auROC_expt{i},hitRst_all,rst_expt{i},hitInv_all_ind,2);
        [hitAuROC_short, hitRst_short] = catIfIndexNotEmpty(hitAuROC_short,auROC_expt{i},hitRst_short,rst_expt{i},hitInv_short_ind,2);
        [hitAuROC_long, hitRst_long] = catIfIndexNotEmpty(hitAuROC_long,auROC_expt{i},hitRst_long,rst_expt{i},hitInv_long_ind,2);
        
        [missAuROC_all, missRst_all] = catIfIndexNotEmpty(missAuROC_all,auROC_expt{i},missRst_all,rst_expt{i},missInv_all_ind,2);
        [missAuROC_short, missRst_short] = catIfIndexNotEmpty(missAuROC_short,auROC_expt{i},missRst_short,rst_expt{i},missInv_short_ind,2);
        [missAuROC_long, missRst_long] = catIfIndexNotEmpty(missAuROC_long,auROC_expt{i},missRst_long,rst_expt{i},missInv_long_ind,2);
        
        [allAuROC_all, allRst_all] = catIfIndexNotEmpty(allAuROC_all,auROC_expt{i},allRst_all,rst_expt{i},allInv_all_ind,2);
        [allAuROC_short, allRst_short] = catIfIndexNotEmpty(allAuROC_short,auROC_expt{i},allRst_short,rst_expt{i},allInv_short_ind,2);
        [allAuROC_long, allRst_long] = catIfIndexNotEmpty(allAuROC_long,auROC_expt{i},allRst_long,rst_expt{i},allInv_long_ind,2);
    end
        
%     get mean across match target trials
    hitVal_short{i} = nanmean(catIndexedTrialsCell(hitVal{i},hitVal_short_ind),3);
    hitInv_short{i} = nanmean(catIndexedTrialsCell(hitInv{i},hitInv_short_ind),3);
    missVal_short{i} = nanmean(catIndexedTrialsCell(missVal{i},missVal_short_ind),3);
    missInv_short{i} = nanmean(catIndexedTrialsCell(missInv{i},missInv_short_ind),3);
    allVal_short{i} = nanmean(catIndexedTrialsCell(allVal{i},allVal_short_ind),3);
    allInv_short{i} = nanmean(catIndexedTrialsCell(allInv{i},allInv_short_ind),3);
    
    hitVal_long{i} = nanmean(catIndexedTrialsCell(hitVal{i},hitVal_long_ind),3);
    hitInv_long{i} = nanmean(catIndexedTrialsCell(hitInv{i},hitInv_long_ind),3);
    missVal_long{i} = nanmean(catIndexedTrialsCell(missVal{i},missVal_long_ind),3);
    missInv_long{i} = nanmean(catIndexedTrialsCell(missInv{i},missInv_long_ind),3);
    allVal_long{i} = nanmean(catIndexedTrialsCell(allVal{i},allVal_long_ind),3);
    allInv_long{i} = nanmean(catIndexedTrialsCell(allInv{i},allInv_long_ind),3);
    
    hitVal_all{i} = nanmean(catIndexedTrialsCell(hitVal{i},hitVal_all_ind),3);
    hitInv_all{i} = nanmean(catIndexedTrialsCell(hitInv{i},hitInv_all_ind),3);
    missVal_all{i} = nanmean(catIndexedTrialsCell(missVal{i},missVal_all_ind),3);
    missInv_all{i} = nanmean(catIndexedTrialsCell(missInv{i},missInv_all_ind),3);
    allVal_all{i} = nanmean(catIndexedTrialsCell(allVal{i},allVal_all_ind),3);
    allInv_all{i} = nanmean(catIndexedTrialsCell(allInv{i},allInv_all_ind),3);
    
    
    % number of trials and cells per condition
    nCellsTrialsCounts
end  
% concatentate all experiments, subtract baseline;

[hitVal_short, hitInv_short] = catExptAndNorm2BL(hitVal_short, hitInv_short, pre_win);

[missVal_short, missInv_short] = catExptAndNorm2BL(missVal_short, missInv_short, pre_win);

[allVal_short, allInv_short] = catExptAndNorm2BL(allVal_short, allInv_short, pre_win);

[hitVal_long, hitInv_long] = catExptAndNorm2BL(hitVal_long, hitInv_long, pre_win);

[missVal_long, missInv_long] = catExptAndNorm2BL(missVal_long, missInv_long, pre_win);

[allVal_long, allInv_long] = catExptAndNorm2BL(allVal_long, allInv_long, pre_win);

[hitVal_all, hitInv_all] = catExptAndNorm2BL(hitVal_all, hitInv_all, pre_win);

[missVal_all, missInv_all] = catExptAndNorm2BL(missVal_all, missInv_all, pre_win);

[allVal_all, allInv_all] = catExptAndNorm2BL(allVal_all, allInv_all, pre_win);

%% plot n trials per cell per condition
% ncond = length(cond_str);
% 
% nTrPerCell = cell(1,ncond);
% for icond = 1:ncond
%     ntc = [];
%    for icatch = 1:ncatch
%        nc = ncells_expt{icond}(icatch);
%        nt = ntrials_expt{icond}(icatch);
%       if nc ~= 0
%          ntt = ones(nc,1)*nt;
%          ntc = cat(1,ntc,ntt);
%       end
%    end
%    nTrPerCell{icond} = ntc;
% end
% 
% edges = round(linspace(1,50,50));
% for icond = 1:ncond/2;
%    val_ind = icond;
%    inv_ind = icond+(ncond/2);
%    figure;
%    suptitle(cond_str{val_ind}(1:end-4))
%    subplot(1,2,1)
%    h = histogram(nTrPerCell{val_ind},edges);
%    figXAxis(h, 'n trials',[])
%    figYAxis(h,'n cells',[]);
%    title('valid')
%    figAxForm(h.Parent)
%    subplot(1,2,2)
%    h = histogram(nTrPerCell{inv_ind},edges);
% %    legend({'val';'inv'})
%    figXAxis(h, 'n trials',[])
%    figYAxis(h,'n cells',[]);
%    title('invalid')
%    figAxForm(h.Parent)
% end
%% sorting
[hitVal_short_sort, sort_ind] = sortCells4Heatmap(hitVal_short,trans_win);
[hitInv_short_sort, sort_ind] = sortCells4Heatmap(hitInv_short,trans_win);

if ~doCellInd
hitRst_short_sort = logical(hitRst_short(fliplr(sort_ind)));
hitAuROC_short_sort = hitAuROC_short(fliplr(sort_ind));
end


[missVal_short_sort, sort_ind] = sortCells4Heatmap(missVal_short,trans_win);
[missInv_short_sort, sort_ind] = sortCells4Heatmap(missInv_short,trans_win);
if ~doCellInd
missRst_short_sort = logical(missRst_short(fliplr(sort_ind)));
missAuROC_short_sort = missAuROC_short(fliplr(sort_ind));
end

[allVal_short_sort, sort_ind] = sortCells4Heatmap(allVal_short,trans_win);
[allInv_short_sort, sort_ind] = sortCells4Heatmap(allInv_short,trans_win);
if ~doCellInd
allRst_short_sort = logical(allRst_short(fliplr(sort_ind)));
allAuROC_short_sort = allAuROC_short(fliplr(sort_ind));
end

[hitVal_long_sort, sort_ind] = sortCells4Heatmap(hitVal_long,trans_win);
[hitInv_long_sort, sort_ind] = sortCells4Heatmap(hitInv_long,trans_win);
if ~doCellInd
hitRst_long_sort = logical(hitRst_long(fliplr(sort_ind)));
hitAuROC_long_sort = hitAuROC_long(fliplr(sort_ind));
end

[missVal_long_sort, sort_ind] = sortCells4Heatmap(missVal_long,trans_win);
[missInv_long_sort, sort_ind] = sortCells4Heatmap(missInv_long,trans_win);
if ~doCellInd
missRst_long_sort = logical(missRst_long(fliplr(sort_ind)));
missAuROC_long_sort = missAuROC_long(fliplr(sort_ind));
end

[allVal_long_sort, sort_ind] = sortCells4Heatmap(allVal_long,trans_win);
[allInv_long_sort, sort_ind] = sortCells4Heatmap(allInv_long,trans_win);
if ~doCellInd
allRst_long_sort = logical(allRst_long(fliplr(sort_ind)));
allAuROC_long_sort = allAuROC_long(fliplr(sort_ind));
end

[hitVal_all_sort, sort_ind] = sortCells4Heatmap(hitVal_all,trans_win);
[hitInv_all_sort, sort_ind] = sortCells4Heatmap(hitInv_all,trans_win);
if ~doCellInd
hitRst_all_sort = logical(hitRst_all(fliplr(sort_ind)));
hitAuROC_all_sort = hitAuROC_all(fliplr(sort_ind));
end

[missVal_all_sort, sort_ind] = sortCells4Heatmap(missVal_all,trans_win);
[missInv_all_sort, sort_ind] = sortCells4Heatmap(missInv_all,trans_win);
if ~doCellInd
missRst_all_sort = logical(missRst_all(fliplr(sort_ind)));
missAuROC_all_sort = missAuROC_all(fliplr(sort_ind));
end

[allVal_all_sort, sort_ind] = sortCells4Heatmap(allVal_all,trans_win);
[allInv_all_sort, sort_ind] = sortCells4Heatmap(allInv_all,trans_win);
if ~doCellInd
allRst_all_sort = logical(allRst_all(fliplr(sort_ind)));
allAuROC_all_sort = allAuROC_all(fliplr(sort_ind));
end
%% variable name matrices

dn_short_val = {'hitVal_short_sort';'missVal_short_sort';'allVal_short_sort'};
dn_short_inv = {'hitInv_short_sort';'missInv_short_sort';'allInv_short_sort'};
dn_long_val = {'hitVal_long_sort';'missVal_long_sort';'allVal_long_sort'};
dn_long_inv = {'hitInv_long_sort';'missInv_long_sort';'allInv_long_sort'};
dn_all_val = {'hitVal_all_sort';'missVal_all_sort';'allVal_all_sort'};
dn_all_inv = {'hitInv_all_sort';'missInv_all_sort';'allInv_all_sort'};

sp_title = {'hits'; 'misses'; 'all'};
trL_title = {'short';'long';'all'};

%% plot params
ms250_fr  = ceil(oneS_fr/4);
bl_fr = ms250_fr;
stim_fr = 2*ms250_fr;
trL_ind = pre_event_frames - bl_fr: pre_event_frames + stim_fr;

tr_tick_fr = (0:bl_fr:bl_fr+stim_fr)+1;
tr_tick_s = chop(-0.25:0.25:0.5,2);

ttMs = chop(linspace(-0.25,0.5,length(trL_ind)),2);
x_axis_lim = [ttMs(1) ttMs(end)];
y_axis_lim = [-0.02 0.1];
resp_axis_lim = [-0.1 0.1];
%% heatmaps

cb_max = 0.1;

short_hm_fig = figure; setFigParams4Print('portrait')
suptitle({titleStr;'short trials, matched cells, sorted by avg resp to val target'})
colormap(brewermap([],'*RdBu'))
long_hm_fig = figure; setFigParams4Print('portrait')
suptitle({titleStr;'long trials, matched cells, sorted by avg resp to val target'})
colormap(brewermap([],'*RdBu'))
all_hm_fig = figure; setFigParams4Print('portrait')
suptitle({titleStr;'all trials, matched cells, sorted by avg resp to val target'})
colormap(brewermap([],'*RdBu'))
sub_hm_fig = figure; setFigParams4Print('portrait')
suptitle({titleStr;'all trials, sorted by avg resp to val target'})
colormap(brewermap([],'*RdBu'))

for iplot = 1:3
    figure(short_hm_fig);
    p1 = 1 + ((iplot-1)*2);
    p2 = 2 + ((iplot-1)*2);

    subplot(3,2, p1)
    d = eval(dn_short_val{iplot});
    if ~isempty(d)
        cell_ind = ~isnan(mean(d,2));
        f = imagesc(d(cell_ind,trL_ind));
        figXAxis(f.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
        figAxForm(f.Parent);
        colorbar
        caxis([-cb_max cb_max])
        title(['valid ' sp_title{iplot}])
        ylabel('n cells')
    end
    
    subplot(3,2, p2)
    d = eval(dn_short_inv{iplot});
    if ~isempty(d)
    % % % cell_ind = ~isnan(mean(d,2));
    f = imagesc(d(cell_ind,trL_ind));
    figXAxis(f.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
    figAxForm(f.Parent);
    colorbar
    caxis([-cb_max cb_max])
    title(['invalid ' sp_title{iplot}])
    ylabel('n cells')
    end
end
print([fnout titleStr '_hm_short'],'-dpdf','-fillpage')

for iplot = 1:3
    figure(long_hm_fig);
    p1 = 1 + ((iplot-1)*2);
    p2 = 2 + ((iplot-1)*2);

    subplot(3,2, p1)
    d = eval(dn_long_val{iplot});
    if ~isempty(d)
    cell_ind = ~isnan(mean(d,2));
    f = imagesc(d(cell_ind,trL_ind));
    figXAxis(f.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
    figAxForm(f.Parent);
    colorbar
    caxis([-cb_max cb_max])
    title(['valid ' sp_title{iplot}])
    ylabel('n cells')
    end

    subplot(3,2, p2)
    d = eval(dn_long_inv{iplot});
    if ~isempty(d)
    % % % cell_ind = ~isnan(mean(d,2));
    f = imagesc(d(cell_ind,trL_ind));
    figXAxis(f.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
    figAxForm(f.Parent);
    colorbar
    caxis([-cb_max cb_max])
    title(['invalid ' sp_title{iplot}])
    ylabel('n cells')
    end
end
print([fnout titleStr '_hm_long'],'-dpdf','-fillpage')

for iplot = 1:3
    figure(all_hm_fig);
    p1 = 1 + ((iplot-1)*2);
    p2 = 2 + ((iplot-1)*2);

    subplot(3,2, p1)
    d = eval(dn_all_val{iplot});
    if ~isempty(d)
    cell_ind = ~isnan(mean(d,2));
    f = imagesc(d(cell_ind,trL_ind));
    figXAxis(f.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
    figAxForm(f.Parent);
    colorbar
    caxis([-cb_max cb_max])
    title(['valid ' sp_title{iplot}])
    ylabel('n cells')
    end

    subplot(3,2, p2)
    d = eval(dn_all_inv{iplot});
    if ~isempty(d)
    % % % cell_ind = ~isnan(mean(d,2));
    f = imagesc(d(cell_ind,trL_ind));
    figXAxis(f.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
    figAxForm(f.Parent);
    colorbar
    caxis([-cb_max cb_max])
    title(['invalid ' sp_title{iplot}])
    ylabel('n cells')
    end
end
print([fnout titleStr '_hm_all'],'-dpdf','-fillpage')

    figure(sub_hm_fig);

    subplot(1,3, 1)
    d1 = eval(dn_all_val{1});
    f = imagesc(d1(:,trL_ind));
    figXAxis(f.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
    figAxForm(f.Parent);
    colorbar
    caxis([-cb_max cb_max])
    title(['valid ' sp_title{iplot}])
    ylabel('n cells')

    subplot(1,3, 2)
    d2 = eval(dn_all_inv{1});
    f = imagesc(d2(:,trL_ind));
    figXAxis(f.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
    figAxForm(f.Parent);
    colorbar
    caxis([-cb_max cb_max])
    title(['invalid ' sp_title{iplot}])
    ylabel('n cells')
    
    subplot(1,3, 3)
    d3 = d1-d2;
    f = imagesc(d3(:,trL_ind));
    figXAxis(f.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
    figAxForm(f.Parent);
    colorbar
    caxis([-cb_max cb_max])
    title('val-inv')
    ylabel('n cells')

print([fnout titleStr '_hm_sub'],'-dpdf','-fillpage')

%% mean timecourse

tc_trL_fig = figure; setFigParams4Print('portrait')
suptitle({titleStr;'mean tc across matched cells, sorted by trial length'})
colormap(brewermap([],'*RdBu'))
tc_all_fig = figure; setFigParams4Print('portrait')
suptitle({titleStr;'mean tc across matched cells'})
colormap(brewermap([],'*RdBu'))

for iplot = 1:3
    figure(tc_trL_fig);
    p1 = 1 + ((iplot-1)*2);
    p2 = 2 + ((iplot-1)*2);

    sp = subplot(3,2, p1);
    d = eval(dn_short_val{iplot});
    if ~isempty(d)
    cell_ind = ~isnan(mean(d,2));
    tc = mean(d(cell_ind,trL_ind),1);
    tc_err = ste(d(cell_ind,trL_ind),1);
    shadedErrorBar(ttMs,tc,tc_err,'k');
    hold on
    d = eval(dn_short_inv{iplot});
    tc = mean(d(cell_ind,trL_ind),1);
    tc_err = ste(d(cell_ind,trL_ind),1);
    shadedErrorBar(ttMs,tc,tc_err,'c');
    figXAxis(sp,'time (s)',x_axis_lim,tr_tick_s,tr_tick_s);
    figYAxis(sp,'dF/F',y_axis_lim);
    figAxForm(sp);
    hold on
    vline(0,'k--')
    title(['short ' sp_title{iplot}])
    end


    sp = subplot(3,2, p2);
    d = eval(dn_long_val{iplot});
    if ~isempty(d)
    cell_ind = ~isnan(mean(d,2));
    tc = mean(d(cell_ind,trL_ind),1);
    tc_err = ste(d(cell_ind,trL_ind),1);
    shadedErrorBar(ttMs,tc,tc_err,'k');
    hold on
    d = eval(dn_long_inv{iplot});
    tc = mean(d(cell_ind,trL_ind),1);
    tc_err = ste(d(cell_ind,trL_ind),1);
    shadedErrorBar(ttMs,tc,tc_err,'c');

    figXAxis(sp,'time (s)',x_axis_lim,tr_tick_s,tr_tick_s);
    figYAxis(sp,'dF/F',y_axis_lim);
    figAxForm(sp);
    hold on
    vline(0,'k--')
    title(['long ' sp_title{iplot}])
    end
end
print([fnout titleStr '_tc_trL'],'-dpdf','-fillpage')

for iplot = 1:3
    figure(tc_all_fig);
%     p1 = 1 + ((iplot-1)*2);
%     p2 = 2 + ((iplot-1)*2);

    sp = subplot(1,3, iplot);
    d = eval(dn_all_val{iplot});
    if ~isempty(d)
    cell_ind = ~isnan(mean(d,2));
    tc = mean(d(cell_ind,trL_ind),1);
    tc_err = ste(d(cell_ind,trL_ind),1);
    shadedErrorBar(ttMs,tc,tc_err,'k');
    hold on
    d = eval(dn_all_inv{iplot});
    tc = mean(d(cell_ind,trL_ind),1);
    tc_err = ste(d(cell_ind,trL_ind),1);
    shadedErrorBar(ttMs,tc,tc_err,'c');


    figXAxis(sp,'time (s)',x_axis_lim,tr_tick_s,tr_tick_s);
    figYAxis(sp,'dF/F',y_axis_lim);
    figAxForm(sp);
    hold on
    vline(0,'k--')
    title([sp_title{iplot}])
    end
end
print([fnout titleStr '_tc_allTr'],'-dpdf','-fillpage')

%% scatter mean respose all cells

scat_trL_fig = figure; setFigParams4Print('portrait')
suptitle({titleStr;'mean resp across matched cells, sorted by trial length'})
colormap(brewermap([],'*RdBu'))
scat_all_fig = figure; setFigParams4Print('portrait')
suptitle({titleStr;'mean resp across matched cells'})
colormap(brewermap([],'*RdBu'))

for iplot = 1:3
    figure(scat_trL_fig);
    p1 = 1 + ((iplot-1)*2);
    p2 = 2 + ((iplot-1)*2);

    sp = subplot(3,2, p1);
    d = eval(dn_short_val{iplot});
    if ~isempty(d)
%     cell_ind = ~isnan(mean(d,2));
    r_val = mean(d(:,trans_win),2);
    d = eval(dn_short_inv{iplot});
    r_inv = mean(d(:,trans_win),2);
    s = scatter(r_val,r_inv,100,'k.');
    hold on
    vl = plot(-20:1:20,-20:1:20,'k--');
    figXAxis(s,'val dF/F',resp_axis_lim);
    figYAxis(s,{trL_title{iplot};'inv dF/F'},resp_axis_lim);
    figAxForm(vl.Parent);   
    title(['short ' sp_title{iplot}])
    end
    
    sp = subplot(3,2, p2);
    d = eval(dn_long_val{iplot});
    if ~isempty(d)
%     cell_ind = ~isnan(mean(d,2));
    r_val = mean(d(:,trans_win),2);
    d = eval(dn_long_inv{iplot});
    r_inv = mean(d(:,trans_win),2);
    s = scatter(r_val,r_inv,100,'k.');
    hold on
    vl = plot(-20:1:20,-20:1:20,'k--');
    figXAxis(s,'val dF/F',resp_axis_lim);
    figYAxis(s,{trL_title{iplot};'inv dF/F'},resp_axis_lim);
    figAxForm(vl.Parent);   
    title(['long ' sp_title{iplot}])
    end
end
print([fnout titleStr '_scat_trL'],'-dpdf','-fillpage')

for iplot = 1:3
    figure(scat_all_fig);
%     p1 = 1 + ((iplot-1)*2);
%     p2 = 2 + ((iplot-1)*2);

    sp = subplot(1,3, iplot);
    d = eval(dn_all_val{iplot});
    if ~isempty(d)
    r_val = mean(d(:,trans_win),2);
    d = eval(dn_all_inv{iplot});
    r_inv = mean(d(:,trans_win),2);
    s = scatter(r_val,r_inv,100,'k.');
    hold on
    vl = plot(-20:1:20,-20:1:20,'k--');
    figXAxis(s,'val dF/F',resp_axis_lim);
    figYAxis(s,{trL_title{iplot};'inv dF/F'},resp_axis_lim);
    figAxForm(vl.Parent);   
    title([sp_title{iplot}])
    end
end
print([fnout titleStr '_scat_allTr'],'-dpdf','-fillpage')

%% auROC group comparison
if ~doCellInd
    auroc_plots 
end

end