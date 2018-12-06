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
expfile = '\\CRASH.dhe.duke.edu\data\home\kevin\Code\axonSz_experiment_list.txt';
fID = fopen(expfile);
head = textscan(fID,'%s%s%s%s%s',1,'delimiter',',');
head = vertcat(head{:});
temp = textscan(fID,'%s%s%s%s%s','delimiter',',','HeaderLines',1);
temp = horzcat(temp{:});
expdata = cell2table(temp,'VariableNames',head);
nExp = size(expdata,1);
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);

rc = behavConstsAV;
fout = fullfile(rc.ashleyAnalysis, 'SizeTuning');
if ~exist(fout)
    mkdir(fout)
end
fprintf(['Size-tuning axon-area comparison analysis - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

%% load each experiment and concatenate data

fprintf('\nBegin loading and concatentating experiment data...\n')
% sizeTuneData
fprintf('Loading sizeTuneData (raw size-tuning data)\n')
sizeTune_all = cell(0);
sizeMean_all = [];
sizeSEM_all = [];
cellDists_all = [];
nCellsExp = zeros(1,nExp);
sizeFits_all = struct([]); % no cells, 4 cons
lbub_fits_all = [];
goodfit_ind_size_all = [];

t = 0;
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    expDate = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = expdata.run_str{i};
    
    fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    
    filename = fullfile(fn,run_str, [mouse '_' expDate '_Tuning.mat']);
    if ~exist(filename, 'file')
        fprintf([[mouse '_' expDate '_Tuning.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'sizeTune', 'tuning_mat', 'cellDists')
    filename = fullfile(fn,run_str, [mouse '_' expDate '_input.mat']);
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
    if ~exist(filename, 'file')
        fprintf([[mouse '_' expDate '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind_size')
    
    sizeMean_all = [sizeMean_all squeeze(tuning_mat(:,1,:))];
    sizeSEM_all = [sizeSEM_all squeeze(tuning_mat(:,2,:))];
    cellDists_all = [cellDists_all cellDists'];
    nCellsExp(i) = length(cellDists);
    sizeFits_all = [sizeFits_all Fit_struct];
    lbub_fits_all = cat(1,lbub_fits_all,lbub_fits);
    tempinds = sum(nCellsExp(1:i-1)) + goodfit_ind_size; % offset by # cells in previous exps
    goodfit_ind_size_all = [goodfit_ind_size_all tempinds];
    
    fprintf('done\n')
end

nCellsTot = sum(nCellsExp);
fprintf('axons: %d cells loaded (%d goodfit_size)\n',nCellsTot,length(goodfit_ind_size_all))
expInd = [];
expArea = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExp(i))];
    expArea = [expArea repmat(expdata(i,4).area, 1,nCellsExp(i))];
end

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

nExp_area = zeros(size(areas));
nCells_area = nExp_area;

close all
choosefig = [1:3];
% choose figs: 1=modelcounts; 2=averagecurves; 3=prefSize and suppInd
legStrs = strings(1,length(areas)); legStrs2=legStrs;
for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_all);
    
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
    ind = intersect(ind,find(cellDists_all<cutoff));
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    nExp_area(i) = nExpi;
    nCells_area(i) = nCellsi;
    sizeMean = sizeMean_all(:,ind);
    sizeSEM = sizeSEM_all(:,ind);
    sizeFits = [];
    for ii = 1:length(ind)
        ix= ind(ii);
        sizeFits(ii).Ftest = sizeFits_all(ix).True.s_.Ftest; %cell
        sizeFits(ii).prefSize = sizeFits_all(ix).True.s_.prefSize; %cell
        sizeFits(ii).suppInd = sizeFits_all(ix).True.s_.suppInd; %cell
    end
    lbub_fits = lbub_fits_all(ind,:,:); %cell,par,val (low up mean true stdev)
    ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
    ism2 = reshape([sizeFits.Ftest],size(sizeFits));
    
    szs_c = cellstr(num2str(chop(szs,2)'))';
    
    legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nCellsi,nExpi);
       
    if sum(choosefig==1)
        figure(1);if i==1;clf;end %figure 1 = proportions of model2 by con
        subplot(1,nArea,i)
        modelcounts = sum(ism2)/nCellsi;
        bar(modelcounts)
        title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        ylabel('Frac. cells M2')
        ylim([0 1])
    end
    
    if sum(choosefig==2)
        figure(2);if i==1;clf;end %figure 2 = average size tuning curves (normalized)
        % change now to normalize all cells, then plot 3 subplots of m1/m2/all
        sizeMean_norm = sizeMean*0; sizeSEM_norm = sizeSEM*0;
        for iCell = 1:nCellsi
            dum = sizeMean(:,iCell); % take all sizeMean values for cell
            %dum = sizeMean(:,nCon,iCell); % only at highest con
            norm = max(dum(:)); % take max of all dF/F's including all cons
            sizeMean_norm(:,iCell) = sizeMean(:,iCell)/norm; % normalize by this max for the individual cell
            sizeSEM_norm(:,iCell) = sizeSEM(:,iCell)/norm;
        end
        sizeMean_normall = mean(sizeMean_norm,2);
        norm = max(sizeMean_normall(:));
        sizeMean_normall = sizeMean_normall/norm;
        sizeSEM_normall = geomean(sizeSEM_norm,2)/norm;
        subplot(1,nArea,i)
        errorbar(szs,sizeMean_normall,sizeSEM_normall)
        title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Size (deg)')
        ylabel('dF/F (norm)')
        ylim([0 1.25])
    end
    
    if sum(choosefig==3)
        figure(3);if i==1;clf;end %figure 3 = prefSize vs con
        subplot(2,1,1)
        prefSize = reshape([sizeFits.prefSize],size(sizeFits));
        prefMean = mean(prefSize);
        prefSEM = std(prefSize)./sqrt(nCellsi);
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
        suppSEM = std(suppInd)./sqrt(nCellsi);
        errorbar(i,suppMean,suppSEM,'o');
        hold on
        suppInd_all = [suppInd_all [suppInd; i.*ones(size(suppInd))]];
        title('Mean Suppression Index by area')
        ylabel('SI')
        xlim([0 nArea+1])
        ylim([0 1])
        axis square
        legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nCellsi,nExpi);
        if i==nArea;legend(legStrs,'location','southeast');end %'location','southoutside','Orientation','horizontal' for bottom
    end
end
figure(1)
print(fullfile(fout, 'axonSize_modelSummary.pdf'),'-dpdf','-bestfit')
figure(2)
print(fullfile(fout, 'axonSize_avgTuningSummary.pdf'),'-dpdf','-bestfit')
figure(3)
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
