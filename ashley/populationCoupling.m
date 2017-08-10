clear all
close all

ds = 'awData_audMod_V1';
rc = behavConstsAV;
eval(ds)
fnout = 'Z:\Analysis\Expt Summaries\population coupling testing';
%%
nexp = size(expt,2);
timebinS1 = 0.12;
timebinS2 = 0.48;
%%
pCpl = cell(1,nexp);
pCpl_FR = cell(1,nexp);
pCpl_1st = cell(1,nexp);
pCpl_2nd = cell(1,nexp);
pCorr_1 = cell(1,nexp);
nCorrMat = cell(1,nexp);
nFRMat = cell(1,nexp);
for iexp = 1:nexp;
        
    subNum = expt(iexp).SubNum;
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    frRateHz = expt(iexp).frame_rate;
    disp([mouse ' ' expDate])
    
    data = [];
    
    for irun = 1:expt(iexp).nrun
        runFolder = expt(iexp).runs(irun,:);
        runTime = expt(iexp).time_mat(irun,:);
        fName = [runFolder '_000_000'];
        
        % load cell timecourses and mworks file
        fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate,runFolder);
        load(fullfile(fn,'timecourses.mat'))
        % combine data
        data = cat(1,data,data_tc_subnp);        
    end
    f = mean(data,1);
%     [n, bin_edge] = histcounts(f);
%     [~, max_n] = max(n);
%     f = bin_edge(max_n);
    dff = bsxfun(@rdivide, bsxfun(@minus,data,f), f);
    
    
    nfr = size(dff,1);
    dff_1st = dff(1:floor(nfr/2),:);
    dff_2nd = dff(floor(nfr/2)+1:end,:);
    
    offset_sample = linspaceNDim(1:1000:nfr,500:1000:nfr,500);
    offset_sample = int64(offset_sample(:));
    dff_samp1 = dff(offset_sample,:);
    dff_samp2 = dff(setdiff(1:nfr,offset_sample),:);
    
    pCpl_1st{iexp} = popCouple(dff_1st,frRateHz);
    pCpl_2nd{iexp} = popCouple(dff_2nd,frRateHz);  
    pCpl_samp1{iexp} = popCouple(dff_samp1,frRateHz);
    pCpl_samp2{iexp} = popCouple(dff_samp2,frRateHz);    
    [pCpl{iexp},pCpl_FR{iexp}] = popCouple(dff,frRateHz);
    
    pCorr_1{iexp} = popCorr(dff,frRateHz,timebinS1);
    pCorr_2{iexp} = popCorr(dff,frRateHz,timebinS2);
    
    [nCorrMat{iexp},nFRMat{iexp}] = neurCorrMat(dff,frRateHz,timebinS1);
end
%%
pCpl_all = cell2mat(cellfun(@zscore, pCpl, 'unif', 0));
% pCpl_FR_all = cell2mat(cellfun(@zscore, pCpl_FR, 'unif', 0));
pCorr_1_all = cell2mat(pCorr_1);


% pCpl_all = cell2mat(pCpl);
pCpl_FR_all = cell2mat(pCpl_FR);

nbins = 10;
fr_bins = linspace(min(pCpl_FR_all),max(pCpl_FR_all),nbins);
[n,~,bin_ind] = histcounts(pCpl_FR_all,fr_bins);
binnedFR = splitapply(@mean, pCpl_FR_all, findgroups(bin_ind));
binnedCorr = splitapply(@mean, pCorr_1_all, findgroups(bin_ind));
binnedCpl = splitapply(@mean, pCpl_all, findgroups(bin_ind));
    
figure;
subplot 131
h = scatter(pCorr_1_all,pCpl_all);
h.MarkerFaceColor = 'k';
h.MarkerEdgeColor = [1,1,1];
figXAxis([],'Correlation to summed population',[])
figYAxis([],'Population coupling',[])
figAxForm([])
subplot 132
h = scatter(pCpl_FR_all,pCorr_1_all);
h.MarkerFaceColor = [0.5 0.5 0.5];
h.MarkerEdgeColor = [1,1,1];
hold on
h = scatter(binnedFR,binnedCorr);
h.MarkerFaceColor = 'k';
figXAxis([],'Firing rate',[])
figYAxis([],'Correlation to summed population',[])
figAxForm([])
subplot 133
h = scatter(pCpl_FR_all,pCpl_all);
h.MarkerFaceColor = [0.5 0.5 0.5];
h.MarkerEdgeColor = [1,1,1];
hold on
h = scatter(binnedFR,binnedCpl);
h.MarkerFaceColor = 'k';
figXAxis([],'Firing rate',[])
figYAxis([],'Population coupling',[])
figAxForm([])
print(fullfile(fnout,'pcpl_vsCorr'),'-dpdf','-fillpage')

%% compare pearson correlations across down-sampling rates
pCorr_2_all = cell2mat(pCorr_2);
figure;
h = scatter(pCorr_1_all,pCorr_2_all);
h.MarkerFaceColor = 'k';
h.MarkerEdgeColor = [1,1,1];
figXAxis([],[num2str(timebinS1) ' S bins'],[])
figYAxis([],[num2str(timebinS2) ' S bins'],[])
figAxForm([])
title('neuron correlation with population')


%%
cpl_lim = [-8 10];
pCpl_1st_all = cell2mat(cellfun(@zscore, pCpl_1st, 'unif', 0));
pCpl_2nd_all = cell2mat(cellfun(@zscore, pCpl_2nd, 'unif', 0));
pCpl_samp1_all = cell2mat(cellfun(@zscore, pCpl_samp1, 'unif', 0));
pCpl_samp2_all = cell2mat(cellfun(@zscore, pCpl_samp2, 'unif', 0));

figure;
subplot 121
h = scatter(pCpl_1st_all,pCpl_2nd_all);
h.MarkerFaceColor = 'k';
h.MarkerEdgeColor = [1,1,1];
hold on
plot(cpl_lim,cpl_lim,'k--')
figXAxis([],'Population coupling 1st half',cpl_lim)
figYAxis([],'Population coupling 2nd half',cpl_lim)
figAxForm([])
title('compare beginning and end')
subplot 122
h = scatter(pCpl_samp1_all,pCpl_samp2_all);
h.MarkerFaceColor = 'k';
h.MarkerEdgeColor = [1,1,1];
hold on
plot(cpl_lim,cpl_lim,'k--')
figXAxis([],'Population coupling dff set 1',cpl_lim)
figYAxis([],'Population coupling dff set 2',cpl_lim)
figAxForm([])
title('interleaved sampling')

print(fullfile(fnout,'pcpl_vsSelf'),'-dpdf','-fillpage')

%%
corr_lim = [0 1];
resp_lim = [];
nc_expt = cellfun(@length, pCpl,'unif',0);

pair_ind = cellfun(@(x) logical(tril(ones(x),1)),nc_expt, 'unif', 0);
nCorrPair = cellfun(@(x,y) x(y), nCorrMat, pair_ind, 'unif', 0);
nFRPair = cellfun(@(x,y) x(y), nFRMat, pair_ind, 'unif', 0);

% nCorr_all = cell2mat(cellfun(@zscore, nCorrPair, 'unif', 0)');
% nFR_all = cell2mat(cellfun(@zscore, nFRPair, 'unif', 0)');
nCorr_all = cell2mat(nCorrPair');
nFR_all = cell2mat(nFRPair');

nbins = 10;
[n,~,bin_ind] = histcounts(nFR_all,fr_bins);
binnedFR = splitapply(@mean, nFR_all, findgroups(bin_ind));
binnedCorr = splitapply(@mean, nCorr_all, findgroups(bin_ind));

figure;
h = scatter(binnedFR,binnedCorr);
h.MarkerFaceColor = 'k';
h.MarkerEdgeColor = [1,1,1];
figXAxis([],'neuron-neuron mean FR',[])
figYAxis([],'neuron-neuron corr',corr_lim)
figAxForm([])

%%
fr_bins = linspace(0,max(nFR_all),nbins);
figure;
for i = 1:nexp
    c = nCorrMat{i};
    cp_ind = logical(tril(ones(nc_expt{i}),-1));
    cp = c(cp_ind);
    frp = nFRMat{i}(cp_ind);

    nbins = 5;
%     fr_bins = linspace(0,max(frp),nbins);
    [n,~,bin_ind] = histcounts(frp,fr_bins);
    binnedFR = splitapply(@mean, frp, findgroups(bin_ind));
    binnedCorr = splitapply(@mean, cp, findgroups(bin_ind));
    
    hold on
    scatter(binnedFR,binnedCorr,100,'.')
end
figXAxis([],'neuron-neuron mean FR',resp_lim)
figYAxis([],'neuron-neuron corr',[0 .4])
figAxForm([])