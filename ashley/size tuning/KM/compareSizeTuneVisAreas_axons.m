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

rc = behavConstsAV;
fout = fullfile(rc.ashleyAnalysis, 'SizeTuning');
if ~exist(fout)
    mkdir(fout)
end

expfile = fullfile(fout, 'axonSz_experiment_list.txt');
fID = fopen(expfile);
head = textscan(fID,'%s%s%s%s%s%s',1,'delimiter',',');
head = vertcat(head{:});
temp = textscan(fID,'%s%s%s%s%s%s','delimiter',',','HeaderLines',1);
temp = horzcat(temp{:});
expdata = cell2table(temp,'VariableNames',head);
nExp = size(expdata,1);
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);

fprintf(['Size-tuning axon-area comparison analysis - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

%% load each experiment and concatenate data

fprintf('\nBegin loading and concatentating experiment data...\n')
% sizeTuneData
fprintf('Loading retData (raw ret data)\n')

nCellsRetExp = zeros(1,nExp);
retFits_all = struct([]); % no cells, 4 cons
lbub_fits_ret_all = [];
goodfit_ind_ret_all = [];
expInd_ret = [];
expArea_ret = [];

for i=1:nExp
    fprintf(['Ret Exp: ' num2str(i) '/' num2str(nExp) '...'])
    expDate = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = expdata.ret_run_str{i};
    
    fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    
    filename = fullfile(fn,run_str, [mouse '_' expDate '_input.mat']);
    if ~exist(filename, 'file')
        fprintf([[mouse '_' expDate '_input.mat'] ' not found! Please remove from list\n'])
    end
    load(filename);
    filename = fullfile(fn,run_str, [mouse '_' expDate '_Tuning.mat']);
    load(filename);
    filename = fullfile(fn,run_str, [mouse '_' expDate '_Fit_struct.mat']);
    if exist(filename, 'file')
        load(filename)
        Fit_struct = rmfield(Fit_struct, 'Shuf');
    end
    filename = fullfile(fn,run_str, [mouse '_' expDate '_Fit_struct_sub.mat']);
    if exist(filename, 'file')
        load(filename)
        Fit_struct = Fit_struct_sub;
    end
    filename = fullfile(fn,run_str, [mouse '_' expDate '_lbub_fits.mat']);
    load(filename, 'lbub_fits', 'goodfit_ind')

    nCellsRetExp(i) = size(lbub_fits,1);
    retFits_all = [retFits_all Fit_struct];
    lbub_fits_ret_all = cat(1,lbub_fits_ret_all,lbub_fits);
    tempinds = sum(nCellsRetExp(1:i-1)) + goodfit_ind; % offset by # cells in previous exps
    goodfit_ind_ret_all = [goodfit_ind_ret_all tempinds];
    
    expInd_ret = [expInd_ret repmat(i,1,nCellsRetExp(i))];
    expArea_ret = [expArea_ret repmat(expdata(i,5).area, 1,nCellsRetExp(i))];
    
    fprintf('done\n')
end

fprintf('Loading sizeTuneData (raw size-tuning data)\n')
sizeMean_all = [];
sizeSEM_all = [];
cellDists_all = [];
nCellsSizeExp = zeros(1,nExp);
sizeFits_all = struct([]); % no cells, 4 cons
lbub_fits_size_all = [];
goodfit_ind_size_all = [];
expInd_size = [];
expArea_size = [];

for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    expDate = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = expdata.size_run_str{i};
    
    fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    
    filename = fullfile(fn,run_str, [mouse '_' expDate '_Tuning.mat']);
    load(filename, 'sizeTune', 'tuning_mat', 'cellDists', 'tc_dfof')
    filename = fullfile(fn,run_str, [mouse '_' expDate '_input.mat']);
    load(filename);
    filename = fullfile(fn,run_str, [mouse '_' expDate '_Fit_struct.mat']);
    if exist(filename, 'file')
        load(filename)
        if isfield(Fit_struct,'Shuf')
            Fit_struct = rmfield(Fit_struct, 'Shuf');
        end
    end
    filename = fullfile(fn,run_str, [mouse '_' expDate '_Fit_struct_sub.mat']);
    if exist(filename, 'file')
        load(filename)
        Fit_struct = Fit_struct_sub;
    end
    filename = fullfile(fn,run_str, [mouse '_' expDate '_lbub_fits.mat']);
    load(filename, 'lbub_fits', 'goodfit_ind_size')
    
    nCellsSizeExp(i) = size(Fit_struct,2);
    sizeMean_all = [sizeMean_all squeeze(tuning_mat(:,1,1:nCellsSizeExp(i)))];
    sizeSEM_all = [sizeSEM_all squeeze(tuning_mat(:,2,1:nCellsSizeExp(i)))];
    cellDists_all = [cellDists_all cellDists(1:nCellsSizeExp(i),:)'];
    
    sizeFits_all = [sizeFits_all Fit_struct];
    lbub_fits_size_all = cat(1,lbub_fits_size_all,lbub_fits(1:nCellsSizeExp(i),:,:));
    tempinds = sum(nCellsSizeExp(1:i-1)) + goodfit_ind_size; % offset by # cells in previous exps
    goodfit_ind_size_all = [goodfit_ind_size_all tempinds];
    
    expInd_size = [expInd_size repmat(i,1,nCellsSizeExp(i))];
    expArea_size = [expArea_size repmat(expdata(i,5).area, 1,nCellsSizeExp(i))];
    
    fprintf('done\n')
end

nCellsRetTot = sum(nCellsRetExp);
nCellsSizeTot = sum(nCellsSizeExp);
fprintf('axons: %d boutons loaded, %d goodfit_ret, %d goodfit_size\n',nCellsRetTot,length(goodfit_ind_ret_all), length(goodfit_ind_size_all))

szs = unique(celleqel2mat_padded(input.tGratingDiameterDeg)); 
nSize = length(szs);
szRng = linspace(0,max(szs));

%% compare areas 
% making figs 1-7
fprintf('Examine cells from each area:\n')
areas = ["LM","AL","PM"];
nArea = length(areas);
prefSize_all = [];
suppInd_all = [];
RFsize = cell(1,nArea);
RFsize_all = [];

nExp_area = zeros(size(areas));
nCells_area = nExp_area;

close all
choosefig = [1:4];
% choose figs: 1= RF size; 2=modelcounts; 3=averagecurves; 4=prefSize and suppInd
retLegStrs = strings(1,length(areas));
sizeLegStrs = strings(1,length(areas)); 
for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind_ret = intersect(find(strcmp(expArea_ret, areas(i))),goodfit_ind_ret_all);
    ind_size = intersect(find(strcmp(expArea_size, areas(i))),goodfit_ind_size_all);
    
    % cutoff by cellDist
    % try looking with different cutoffs
    switch areas(i)
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
            excells = [631 2128];
            % case {'LM','AL'}
            %     cutoff = 15; %alm cutoff at 15
        case 'LM'
            cutoff = 10; %lm cutoff at 15
            excells = [1861 1863];
        case 'AL'
            cutoff = 10; %al cutoff at 15
            excells = [2395 1777];
        case 'PM'
            cutoff = 10; %pm cutoff at 20
            excells = [1952 2292];
    end
    ind_size = intersect(ind_size,find(cellDists_all<cutoff));
    
    nExpi = length(expIndi);
    nMousei = length(unique(expdata.mouse(expIndi)));
    nCells_Reti = length(ind_ret);
    nCells_Sizei = length(ind_size);
    nExp_area(i) = nExpi;
    nMouse_area(i) = nMousei;
    nCells_area_ret(i) = nCells_Reti;
    nCells_area_size(i) = nCells_Sizei;
    
    sigmax = lbub_fits_ret_all(ind_ret,2,4);
    sigmay = lbub_fits_ret_all(ind_ret,3,4);
    RFsize{i} = sqrt(2*log(2))*geo_mean([sigmax sigmay],2)';
    
    sizeMean = sizeMean_all(:,ind_size);
    sizeSEM = sizeSEM_all(:,ind_size);    
    sizeFits = struct;
    t = [];
    for ii = 1:length(ind_size)
        ix= ind_size(ii);
        if ~isempty(sizeFits_all(ix).True)
            sizeFits(ii).Ftest = sizeFits_all(ix).True.s_.Ftest; %cell
            sizeFits(ii).prefSize = sizeFits_all(ix).True.s_.prefSize; %cell
            sizeFits(ii).suppInd = sizeFits_all(ix).True.s_.suppInd; %cell
            sizeFits(ii).fitout1 = sizeFits_all(ix).True.s_.fitout1; %cell
            sizeFits(ii).fitout2 = sizeFits_all(ix).True.s_.fitout2; %cell
        else
            t = [t ix];
            sizeFits(ii).Ftest = NaN; %cell
            sizeFits(ii).prefSize = NaN; %cell
            sizeFits(ii).suppInd = NaN; %cell
            sizeFits(ii).fitout1 = NaN; %cell
            sizeFits(ii).fitout2 = NaN; %cell
        end
    end
    ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
    ism2 = reshape([sizeFits.Ftest],size(sizeFits));
    
    szs_c = cellstr(num2str(chop(szs,2)'))';
    
    retLegStrs(i)=sprintf('%s (n=%d, n_{exp}=%d, n_{mice}=%d)',areas(i),nCells_Reti,nExpi,nMousei);
    
    if sum(choosefig==1)
        figure(1);if i==1;clf;end %figure 1 = proportions of model2 by con
        subplot(2,1,1)
        errorbar(i, mean(RFsize{i},2), std(RFsize{i},[],2),'o');%./length(ind_ret),'o');
        RFsize_all = [RFsize_all [RFsize{i}; i.*ones(size(RFsize{i}))]];
        title('Average RF size')
        ylabel('RF size (deg)')
        xlim([0 nArea+1])
        hold on
        if i==nArea;legend(retLegStrs,'location','southeast');end %'location','southoutside','Orientation','horizontal' for bottom
    end
    
    if sum(choosefig==2)
        figure(2);if i==1;clf;end %figure 1 = proportions of model2 by con
        subplot(1,nArea,i)
        modelcounts = sum(ism2)/nCells_Sizei;
        bar(modelcounts)
        title({sprintf('Area:%s',areas(i));['(n=' num2str(nCells_Sizei) ', n_{exp}=' num2str(nExpi) ', n_{mice}=' num2str(nMousei) ')']})
        ylabel('Frac. cells M2')
        ylim([0 1])
    end
    
    if sum(choosefig==3)
        figure(3);if i==1;clf;end %figure 2 = average size tuning curves (normalized)
        % change now to normalize all cells, then plot 3 subplots of m1/m2/all
        sizeMean_norm = sizeMean*0; sizeSEM_norm = sizeSEM*0;
        for iCell = 1:nCells_Sizei
            dum = sizeMean(:,iCell); % take all sizeMean values for cell
            %dum = sizeMean(:,nCon,iCell); % only at highest con
            norm = max(dum(:)); % take max of all dF/F's including all cons
            sizeMean_norm(:,iCell) = sizeMean(:,iCell)/norm; % normalize by this max for the individual cell
            sizeSEM_norm(:,iCell) = sizeSEM(:,iCell)/norm;
        end
        sizeMean_normall = mean(sizeMean_norm,2);
        norm = max(sizeMean_normall(:));
        sizeMean_normall = sizeMean_normall/norm;
        sizeSEM_normall = geo_mean(sizeSEM_norm,2)/norm;
        subplot(2,nArea,i)
        errorbar(szs,sizeMean_normall,sizeSEM_normall)
        title({sprintf('Area:%s',areas(i));['(n=' num2str(nCells_Sizei) ', n_{exp}=' num2str(nExpi) ', n_{mice}=' num2str(nMousei) ')']})
        xlabel('Size (deg)')
        ylabel('dF/F (norm)')
        ylim([0 1.25])
        subplot(2,nArea,i+3)
        sizeFit1 = reshape([sizeFits.fitout1], [size(szRng,2) size(sizeFits,2)]);
        sizeFit2 = reshape([sizeFits.fitout2], [size(szRng,2) size(sizeFits,2)]);
        sizeFit_all = [sizeFit1(:,find(ism1)) sizeFit2(:,find(ism2))];
        sizeFit_all_norm = sizeFit_all./(max(sizeFit_all,[],1));
        shadedErrorBar(szRng,mean(sizeFit_all_norm,2),std(sizeFit_all_norm,[],2)./sqrt(size(sizeFit_all,2)));
        xlabel('Size (deg)')
        ylabel('dF/F (norm)')
        ylim([0 1.25])
    end
    
    if sum(choosefig==4)
        figure(4);if i==1;clf;end %figure 3 = prefSize vs con
        subplot(2,1,1)
        prefSize = reshape([sizeFits.prefSize],size(sizeFits));
        prefMean = mean(prefSize);
        prefSEM = std(prefSize);%./sqrt(nCells_Sizei);
        errorbar(i, prefMean,prefSEM,'o');
        hold on
        prefSize_all = [prefSize_all [prefSize; i.*ones(size(prefSize))]];
        title('Mean Preferred Size by Area')
        ylabel('PrefSize')
        xlim([0 nArea+1])
        ylim([0 60])
        axis square
        subplot(2,1,2)
        suppInd = reshape([sizeFits.suppInd],size(sizeFits));
        suppInd(suppInd<0)=0;suppInd(suppInd>1)=1;
        suppMean = mean(suppInd);
        suppSEM = std(suppInd);%./sqrt(nCells_Sizei);
        errorbar(i,suppMean,suppSEM,'o');
        hold on
        suppInd_all = [suppInd_all [suppInd; i.*ones(size(suppInd))]];
        title('Mean Suppression Index by area')
        ylabel('SI')
        xlim([0 nArea+1])
        ylim([0 1])
        axis square
        sizeLegStrs(i)=sprintf('%s (n=%d, n_{exp}=%d, n_{mice}=%d)',areas(i),nCells_Sizei,nExpi,nMousei);
        if i==nArea;legend(sizeLegStrs,'location','southeast');end %'location','southoutside','Orientation','horizontal' for bottom
    end
    eval(strcat('suppInd_axons_', areas(i), ' = suppInd;'))
    eval(strcat('prefSize_axons_', areas(i), ' = prefSize;'))
    eval(strcat('RFsize_axons_', areas(i), ' = RFsize{i};'))
end
save(fullfile(fout, 'axonSizeSummary.mat'),'suppInd_axons_AL','suppInd_axons_LM','suppInd_axons_PM','prefSize_axons_AL','prefSize_axons_LM','prefSize_axons_PM','RFsize_axons_AL','RFsize_axons_LM','RFsize_axons_PM')

figure(1)
subplot(2,1,2)
boxplot(RFsize_all(1,:),RFsize_all(2,:))
hold on
plot(1:3,cellfun(@mean, RFsize),'x')
[p_RFsize t_RFsize stats_RFsize] = anovan(RFsize_all(1,:),{RFsize_all(2,:)},'Display','off');
ylim([0 20])
ylabel('RF size(deg)')
text(0.5,15,['p = ' num2str(chop(p_RFsize,2))])
if p_RFsize<0.05
    RFsize_out = multcompare(stats_RFsize,'Display','off');
end
print(fullfile(fout, 'axonRet_sizeSummary.pdf'),'-dpdf','-bestfit')
figure(2)
print(fullfile(fout, 'axonSize_modelSummary.pdf'),'-dpdf','-bestfit')
figure(3)
print(fullfile(fout, 'axonSize_avgTuningSummary.pdf'),'-dpdf','-bestfit')
figure(4)
[p_size t_size stats_size] = anovan(prefSize_all(1,:),{prefSize_all(2,:)},'Display','off');
subplot(2,1,1)
text(0.5,50,['p = ' num2str(chop(p_size,2))])
if p_size<0.05
    size_out = multcompare(stats_size,'Display','off');
end
[p_supp t_supp stats_supp] = anovan(suppInd_all(1,:),{suppInd_all(2,:)},'Display','off');
subplot(2,1,2)
text(0.5,0.9,['p = ' num2str(chop(p_supp,2))])
if p_supp<0.05
    supp_out = multcompare(stats_supp,'Display','off');
end
print(fullfile(fout, 'axonSize_SizeSuppSummary.pdf'),'-dpdf','-bestfit')



