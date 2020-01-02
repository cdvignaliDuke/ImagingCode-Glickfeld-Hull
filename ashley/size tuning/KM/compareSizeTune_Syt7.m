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
fn_in = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';
fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Syt7';
expfile = fullfile(fn_out, 'syt7Sz_experiment_list.txt');
fID = fopen(expfile);
head = textscan(fID,'%s%s%s%s%s%s',1,'delimiter',',');
head = vertcat(head{:});
temp = textscan(fID,'%s%s%s%s%s%s','delimiter',',','HeaderLines',1);
temp = horzcat(temp{:});
expdata = cell2table(temp,'VariableNames',head);
nExp = size(expdata,1);
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);
fprintf(['Size-tuning Syt7 comparison analysis - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

%% load each experiment and concatenate data

fprintf('\nBegin loading and concatentating experiment data...\n')
% sizeTuneData
fprintf('Loading sizeTuneData (raw size-tuning data)\n')
sizeTune_all = cell(0);
sizeMean_all = [];
sizeSEM_all = [];
cellDists_all = [];
nCellsExp = zeros(1,nExp);
sizeFits_all = struct([]); 
lbub_fits_all = [];
goodfit_ind_size_all = [];

t = 0;
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    size_run_str = expdata.size_run_str{i};
    filename = fullfile(fn_in, [date '_i' mouse], [size_run_str], [mouse '_' date '_Tuning.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' size_run_str '_sizeTuneData.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'tuning_mat', 'sizeTune', 'cellDists')
    sizeMean = squeeze(tuning_mat(:,:,1,:));
    sizeSEM = squeeze(tuning_mat(:,:,2,:));
    filename = fullfile(fn_in, [date '_i' mouse], [size_run_str], [mouse '_' date '_input.mat']);
    load(filename);
    filename = fullfile(fn_in, [date '_i' mouse], [size_run_str], [mouse '_' date '_Fit_struct.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' size_run_str '_Fit_struct.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'Fit_struct')
    filename = fullfile(fn_in, [date '_i' mouse], [size_run_str], [mouse '_' date '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' size_run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind_size')
    
    
    sizeTune_all = cat(3,sizeTune_all,sizeTune);
    sizeMean_all = cat(3,sizeMean_all,sizeMean);
    sizeSEM_all = cat(3,sizeSEM_all,sizeSEM);
    cellDists_all = [cellDists_all;cellDists];
    nCellsExp(i) = length(cellDists);
    sizeFits_all = cat(1,sizeFits_all,Fit_struct);
    lbub_fits_all = cat(1,lbub_fits_all,lbub_fits);
    tempinds = sum(nCellsExp(1:i-1)) + goodfit_ind_size; % offset by # cells in previous exps
    goodfit_ind_size_all = [goodfit_ind_size_all tempinds];
    
    fprintf('done\n')
end

nCellsTot = sum(nCellsExp);
fprintf('%d cells loaded (%d goodfit_size)\n',nCellsTot,length(goodfit_ind_size_all))
expInd = [];
expArea = [];
expGenotype = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExp(i))];
    expArea = [expArea repmat(expdata(i,4).area, 1,nCellsExp(i))];
    expGenotype = [expGenotype repmat(expdata(i,6).genotype, 1,nCellsExp(i))];
end

szs = unique(celleqel2mat_padded(input.tGratingDiameterDeg)); 
nSize = length(szs);
szRng = linspace(0,max(szs));

clear Fit_struct
%% find low con responsive cells
p = zeros(1,nCellsTot);
for iCell = 1:nCellsTot
    x = [];
    y = [];
    for iS = 1:nSize
        x = [x; iS.*ones(size(cell2mat(sizeTune_all(iS,1,iCell))))];
        y = [y; cell2mat(sizeTune_all(iS,1,iCell))];
    end
    [p(iCell)] = anova1(y,x,'off');
end
ind_lowcon_resp = find(p<0.05);
        
%% con resp
if ~exist(fullfile(fn_out, 'conStruct23.mat'))
fprintf('Extracting contrast response of each cell at prefSize\n')
cons = unique(celleqel2mat_padded(input.tGratingContrast));
nCon = length(cons);
conStruct_all = struct('resp',zeros(1,nCon),'fit',zeros(1,nCon),'C50r',zeros(1),'Rsq',zeros(1),'x0',zeros(1,nCon));
conStruct_all(nCellsTot) = conStruct_all;
conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));
conRng = 0:0.001:1;
opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
usePrefSize = 1;
useRFSize = 0;

for iCell=1:nCellsTot
    if ~sum(iCell==goodfit_ind_size_all)
        if sum(iCell==[1 nCellsTot]) % check first and last to reset zeros to blank
            conStruct_all(iCell).resp = [];conStruct_all(iCell).fit = [];conStruct_all(iCell).C50r = [];conStruct_all(iCell).Rsq = [];conStruct_all(iCell).x0 = [];
        end
        continue % do not fit unless goodfit_size
    end
    
    if usePrefSize
        pS = sizeFits_all(iCell,nCon).True.s_.prefSize;
        pSind = find(szRng==pS);
    elseif useRFSize
        switch cell2mat(expArea(iCell))
            case 'V1'
                RFsize = 10; %v1 cutoff at 10
            case 'LM'
                RFsize = 15; %lm cutoff at 15
            case 'AL'
                RFsize = 15; %al cutoff at 15
            case 'PM'
                RFsize = 20; %pm cutoff at 20
        end
        pS = RFsize;
        [val pSind] = min(abs(szRng-pS));
    end
    
    for iCon = 1:nCon
        if sizeFits_all(iCell,iCon).True.s_.Ftest
            conStruct_all(iCell).resp(iCon) = sizeFits_all(iCell,iCon).True.s_.fitout2(pSind);
        else
            conStruct_all(iCell).resp(iCon) = sizeFits_all(iCell,iCon).True.s_.fitout1(pSind);
        end
    end
    
    cRi = conStruct_all(iCell).resp;
    lb = [0 0 0.1 1];
    ub = [Inf Inf 0.8 Inf];
    SStot = sum((cRi-mean(cRi)).^2);
    R2best = -Inf;
    for i=1%1:4
        x0 = [cRi(1) max(cRi) 0.1+0.1*i 3]; %BL Rmax C50 n
        [cF, res] = lsqcurvefit(conModelH,x0,cons,cRi,lb,ub,opts);
        R2 = 1-res/SStot;
        if R2>R2best
            R2best = R2;
            cFbest = cF;
            x0best = x0;
        end
    end
    cF = cFbest;
    R2 = R2best;
    
    conStruct_all(iCell).fit = cF;
    conStruct_all(iCell).Rsq = R2;
    conStruct_all(iCell).x0 = x0best;
    
    fitout = conModelH(cF,conRng);
    R50 = fitout(1)+(fitout(end)-fitout(1))/2;
    fitout50rect = abs(fitout - R50);
    i50 = find(fitout50rect == min(fitout50rect),1);
    C50 = conRng(i50);
    conStruct_all(iCell).C50r = C50;
    
    fprintf('Cell %d/%d fit: BL=%.3f Rmax=%.3f C50=%.3f n=%.2f : Rsq=%.3f C50r=%.3f\n',iCell,nCellsTot,cF(1),cF(2),cF(3),cF(4),R2,C50)
   
end
fprintf('Done, saving...\n')
filename = fullfile(fn_out, 'conStruct23.mat');
save(filename,'conStruct_all');
else
%% load contrast response instead of compute
filename = fullfile(fn_out, 'conStruct23.mat');
load(filename);
conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));
cons = unique(celleqel2mat_padded(input.tGratingContrast));
nCon = length(cons);
% 
% %% present each exp to examine cells+fits and choose example cells
% % show each cell rawdata w/ model1+2 overlaid

end
%% compare genotypes
% making figs 1-7
genotypes = unique(expGenotype);
nGenotype = length(genotypes);

nExp_genotype = zeros(1,nGenotype);
nCells_genotype = zeros(1,nGenotype);

area = 'V1';

close all
choosefig = [1:6];
% choose figs: 1=modelcounts; 2=averagecurves; 3=prefSize; 4=suppInd; 5=conresp; 6=ex.cells; 7=medianfits;
legStrs = strings(1,nGenotype); legStrs2=legStrs;
for i = 1:nGenotype
    fprintf(['Genotype #' num2str(i) ' : ' char(genotypes(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,genotypes(i)), expdata.genotype, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_all);
    
    
    % cutoff by cellDist
    % try looking with different cutoffs
    switch area
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
        case 'LM'
            cutoff = 15; %lm cutoff at 15
        case 'AL'
            cutoff = 15; %al cutoff at 15
        case 'PM'
            cutoff = 20; %pm cutoff at 20
    end
    ind = intersect(ind,find(cellDists_all<cutoff));
    ind_lowcon = intersect(ind,ind_lowcon_resp);
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    nExp_genotype(i) = nExpi;
    nCells_genotype(i) = nCellsi;
    sizeTune = sizeTune_all{:,:,ind}; % (size,con,cell)
    sizeMean = sizeMean_all(:,:,ind);
    sizeSEM = sizeSEM_all(:,:,ind);
    sizeFits = sizeFits_all(ind,:); %cell,con
    lbub_fits = lbub_fits_all(ind,:,:); %cell,par,val (low up mean true stdev)
    
    Ftest = zeros(nCellsi,nCon);
    prefSize = zeros(nCellsi,nCon);
    suppInd = zeros(nCellsi,nCon);
    for icell = 1:nCellsi
        for icon = 1:nCon
            Ftest(icell,icon) = sizeFits(icell,icon).True.s_.Ftest;
            prefSize(icell,icon) = sizeFits(icell,icon).True.s_.prefSize;
            suppInd(icell,icon) = sizeFits(icell,icon).True.s_.suppInd;
        end
    end
    
    ism1 = ~Ftest;
    ism2 = Ftest;
    
    nCellsi_lowcon = length(ind_lowcon);
    Ftest_lowcon = zeros(nCellsi_lowcon,nCon);
    prefSize_lowcon = zeros(nCellsi_lowcon,nCon);
    suppInd_lowcon = zeros(nCellsi_lowcon,nCon);
    for icell = 1:nCellsi_lowcon
        for icon = 1:nCon
            Ftest_lowcon(icell,icon) = sizeFits(icell,icon).True.s_.Ftest;
            prefSize_lowcon(icell,icon) = sizeFits(icell,icon).True.s_.prefSize;
            suppInd_lowcon(icell,icon) = sizeFits(icell,icon).True.s_.suppInd;
        end
    end
    
    ism1_lowcon = ~Ftest_lowcon;
    ism2_lowcon = Ftest_lowcon;
    
    cons_c = categorical(cellstr(num2str(cons'))');
    szs_c = cellstr(num2str(chop(szs,2)'))';
    conStruct = conStruct_all(ind);
    
    legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',cell2mat(genotypes(i)),nCellsi,nExpi);
    
    if sum(choosefig==1)
        figure(1);if i==1;clf;end %figure 1 = proportions of model2 by con
        subplot(1,nGenotype,i)
        modelcounts = [sum(ism1); sum(ism2)]'/nCellsi;
        bar(cons_c,modelcounts,'stacked')
        title({sprintf('Genotype:%s',cell2mat(genotypes(i)));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Contrast')
        ylabel('Frac. cells')
        if i==nGenotype;legend('m1','m2','location','best');end
    end
    
    if sum(choosefig==2)
        figure(2);if i==1;clf;end %figure 2 = average size tuning curves (normalized)
        % change now to normalize all cells, then plot 3 subplots of m1/m2/all
        sizeMean_norm = sizeMean*0; sizeSEM_norm = sizeSEM*0;
        for iCell = 1:nCellsi
            dum = sizeMean(:,:,iCell); % take all sizeMean values for cell
            %dum = sizeMean(:,nCon,iCell); % only at highest con
            norm = max(dum(:)); % take max of all dF/F's including all cons
            sizeMean_norm(:,:,iCell) = sizeMean(:,:,iCell)/norm; % normalize by this max for the individual cell
            sizeSEM_norm(:,:,iCell) = sizeSEM(:,:,iCell)/norm;
        end
        sizeMean_normall = mean(sizeMean_norm,3);
        norm = max(sizeMean_normall(:,nCon));
        sizeMean_normall = sizeMean_normall/norm;
        sizeSEM_normall = geomean(sizeSEM_norm,3)/norm;
        % split by model
        %subplot(4,3,3*(i-1)+1)
        %subplot(2,4,2*(i-1)+1)
        % for iCon = 1:nCon
        %     errorbar(szs,mean(sizeMean_norm(:,iCon,find(ism1(:,iCon))),3),geomean(sizeSEM_norm(:,iCon,find(ism1(:,iCon))),3))
        %     hold on
        % end
        % title({sprintf('Model1: Area:%s',areas(i));['(n=' num2str(mean(sum(ism1))) ', n_{exp}=' num2str(nExpi) ')']})
        % xlabel('Size (deg)')
        % ylabel('dF/F (norm)')
        % ylim([0 1.2])
        % subplot(2,4,2*(i-1)+2)
        % for iCon = 1:nCon
        %     errorbar(szs,mean(sizeMean_norm(:,iCon,find(ism2(:,iCon))),3),geomean(sizeSEM_norm(:,iCon,find(ism2(:,iCon))),3))
        %     hold on
        % end
        % title({sprintf('Model2: Area:%s',areas(i));['(n=' num2str(mean(sum(ism2))) ', n_{exp}=' num2str(nExpi) ')']})
        % xlabel('Size (deg)')
        % ylabel('dF/F (norm)')
        % ylim([0 1.2])
        % collapse models
        subplot(1,nGenotype,i)
        for iCon = 1:nCon
            errorbar(szs,sizeMean_normall(:,iCon),sizeSEM_normall(:,iCon))
            hold on
        end
        title({sprintf('Area:%s',cell2mat(genotypes(i)));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Size (deg)')
        ylabel('dF/F (norm)')
        ylim([0 1.25])
        if i==nGenotype;legend(num2str(cons'),'location','best');end
    end
    
    if sum(choosefig==3)
        figure(3);if i==1;clf;end %figure 3 = prefSize vs con
        % subplot(2,2,i)
        % prefMean1=zeros(1,nCon);prefSEM1=prefMean1;prefMean2=prefMean1;prefSEM2=prefMean1;prefMeanAll=prefMean1;prefSEMAll=prefMean1;
        % for iCon=1:nCon
        %     prefMean1(iCon) = mean(prefSize(find(ism1(:,iCon)),iCon));
        %     prefSEM1(iCon) = std(prefSize(find(ism1(:,iCon)),iCon))./sqrt(sum(ism1(:,iCon)));
        %     prefMean2(iCon) = mean(prefSize(find(ism2(:,iCon)),iCon));
        %     prefSEM2(iCon) = std(prefSize(find(ism2(:,iCon)),iCon))./sqrt(sum(ism2(:,iCon)));
        %     prefMeanAll(iCon) = mean(prefSize(:,iCon));
        %     prefSEMAll(iCon) = std(prefSize(:,iCon))./sqrt(nCellsi);
        % end
        % errorbar(cons,prefMean1,prefSEM1,'s-');
        % hold on
        % errorbar(cons,prefMean2,prefSEM2,'^-');
        % errorbar(cons,prefMeanAll,prefSEMAll,'kx-');
        % hold off
        % title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        % xlabel('Contrast')
        % ylabel('PrefSize')
        % xlim([0 1])
        % ylim([0 60])
        % if i==4;legend('m1','m2','all','location','best');end
        prefMean=zeros(1,nCon);prefSEM=prefMean;
        for iCon=1:nCon
            prefMean(iCon) = mean(prefSize(:,iCon));
            prefSEM(iCon) = std(prefSize(:,iCon))./sqrt(nCellsi);
        end
        errorbar(cons,prefMean,prefSEM);
        hold on
        title('Mean Preferred Size by Area')
        xlabel('Contrast')
        ylabel('PrefSize')
        xlim([0 1])
        ylim([0 60])
        if i==nGenotype;legend(legStrs,'location','eastoutside');end %'location','southoutside','Orientation','horizontal' for bottom
    end
    
    if sum(choosefig==4)
        figure(4);if i==1;clf;end %figure 4 = suppInd vs con
        suppInd(suppInd<0)=0;suppInd(suppInd>1)=1;
        % subplot(2,2,i)
        % suppMean2=zeros(1,nCon);suppSEM2=suppMean2;suppMeanAll=suppMean2;suppSEMAll=suppMean2;
        % for iCon=1:nCon
        %     suppMean2(iCon) = mean(suppInd(find(ism2(:,iCon)),iCon));
        %     suppSEM2(iCon) = std(suppInd(find(ism2(:,iCon)),iCon))./sqrt(sum(ism2(:,iCon)));
        %     suppMeanAll(iCon) = mean(suppInd(:,iCon));
        %     suppSEMAll(iCon) = std(suppInd(:,iCon))./sqrt(nCellsi);
        % end
        % errorbar(cons,suppMean2,suppSEM2,'^-');
        % hold on
        % errorbar(cons,suppMeanAll,suppSEMAll,'kx-');
        % hold off
        % title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        % xlabel('Contrast')
        % ylabel('Supp Ind')
        % xlim([0 1])
        % ylim([0 1])
        % if i==4;legend('m2 only','all','location','best');end
        
        suppMean=zeros(1,nCon);suppSEM=suppMean;
        for iCon=1:nCon
            suppMean(iCon) = mean(suppInd(:,iCon));
            suppSEM(iCon) = std(suppInd(:,iCon))./sqrt(nCellsi);
        end
        errorbar(cons,suppMean,suppSEM);
        hold on
        title('Mean Suppression Index by area')
        xlabel('Contrast')
        ylabel('SI')
        xlim([0 1])
        ylim([0 1])
        legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',cell2mat(genotypes(i)),nCellsi,nExpi);
        if i==nGenotype;legend(legStrs,'location','eastoutside');end %'location','southoutside','Orientation','horizontal' for bottom
    end
    
    if sum(choosefig==5) %figure 5: average contrast response in each area
        conRng = 0.001:0.001:1;
        opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
        cut = find([conStruct.Rsq]>0.9);
        legStrs2(i)=sprintf('%s (n=%d)',cell2mat(genotypes(i)),length(cut));
        conResp = reshape([conStruct(cut).resp],nCon,length(cut))';
        conResp_norm = conResp./conResp(:,nCon);
        conMean = mean(conResp_norm,1);
        conSEM = std(conResp_norm,[],1)./sqrt(length(cut));
        figure(5);if i==1;clf;end
        ax = gca;
        ax.ColorOrderIndex = i;
        %subplot(2,2,i)
        %for iCell = 1:nCellsi
        %    p1 = plot(cons,conResp_norm(iCell,:),'r-');
        %    p1.Color(4) = 0.1;
        %    hold on
        %end
        hold on
        errorbar(cons,conMean,conSEM)
        %title({sprintf('Contrast response - Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        title('Mean contrast response by area')
        xlabel('Contrast')
        ylabel('norm. dF/F @ pref size')
        xlim([0 1])
        ylim([0 1.2])
        if i==nGenotype;legend(legStrs2,'location','southoutside','Orientation','horizontal');end %'location','southoutside','Orientation','horizontal' for bottom
    
        % fit
        conResp_norm = conResp_norm';
        cRi = conResp_norm(:);
        cons_exp = repmat(cons,1,length(cut));
        lb = [0 0 0.1 1];
        ub = [Inf Inf 0.8 Inf];
        SStot = sum((cRi-mean(cRi)).^2);
        R2best = -Inf;
        x0 = [cRi(1) mean(cRi) 0.2 3]; %BL Rmax C50 n
        [cF, res] = lsqcurvefit(conModelH,x0,cons_exp',cRi,lb,ub,opts);
        R2 = 1-res/SStot;
        
        fitout = conModelH(cF,conRng);
        R50 = fitout(1)+(fitout(end)-fitout(1))/2;
        fitout50rect = abs(fitout - R50);
        i50 = find(fitout50rect == min(fitout50rect),1);
        C50 = conRng(i50);
        
        ax = gca;
        ax.ColorOrderIndex = i;
        plot(conRng,fitout,':','HandleVisibility','off')
        ax = gca;
        ax.ColorOrderIndex = i;
        plot(C50,R50,'x','HandleVisibility','off')
        ax = gca;
        ax.ColorOrderIndex = i;
        plot([C50 C50],[0 R50],'--','HandleVisibility','off')
    end
    
    if sum(choosefig==6)
        figure(6);if i==1;clf;end %figure 1 = proportions of model2 by con
        subplot(1,nGenotype,i)
        modelcounts = [sum(ism1_lowcon); sum(ism2_lowcon)]'/nCellsi_lowcon;
        bar(cons_c,modelcounts,'stacked')
        title({sprintf('Genotype:%s',cell2mat(genotypes(i)));['(n=' num2str(nCellsi_lowcon) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Contrast')
        ylabel('Frac. cells')
        if i==nGenotype;legend('m1','m2','location','best');end
    end
    
end

%% print
figure(1)
print(fullfile(fn_out, 'syt7_FractSuppressed.pdf'),'-dpdf','-bestfit')
figure(2)
print(fullfile(fn_out, 'syt7_TuningCurves.pdf'),'-dpdf','-bestfit')
figure(3)
print(fullfile(fn_out, 'syt7_PrefSize.pdf'),'-dpdf','-bestfit')
figure(4)
print(fullfile(fn_out, 'syt7_SuppIndex.pdf'),'-dpdf','-bestfit')
figure(5)
print(fullfile(fn_out, 'syt7_ContrastResp.pdf'),'-dpdf','-bestfit')


