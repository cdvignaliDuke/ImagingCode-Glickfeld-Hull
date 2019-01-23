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
expfile = '\\CRASH.dhe.duke.edu\data\home\kevin\Code\Ai9x_experiment_list_6conOnly.txt';
fID = fopen(expfile);
head = textscan(fID,'%s%s%s%s%s',1,'delimiter',',');
head = vertcat(head{:});
temp = textscan(fID,'%s%s%s%s%s','delimiter',',','HeaderLines',1);
temp = horzcat(temp{:});
expdata = cell2table(temp,'VariableNames',head);
nExp = size(expdata,1);
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);

fprintf(['Size-tuning visual-area comparison analysis (6 con) - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

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
sizeFits_all = struct([]); % no cells, 6 cons
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

nCellsTot = sum(nCellsExp);
fprintf('%d cells loaded (%d goodfit_size)\n',nCellsTot,length(goodfit_ind_size_all))

expInd = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExp(i))];
end

cons = [0.05 0.1 0.2 0.4 0.8 1]; nCon = length(cons);
szs = 5*1.5.^[0:7]; nSize = length(szs);
szRng = linspace(0,max(szs));

%% con resp for all cells
override = 0;
conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));
conRng = 0.001:0.001:1;
opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
if exist(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', 'conStruct_6con.mat'),'file') && ~override
    % load contrast response instead of compute
    fprintf('Found previous fits, loading...\n')
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', 'conStruct_6con.mat');
    load(filename);
else
    fprintf('Extracting contrast response of each cell at prefSize\n')
    s4 = zeros(1,4);
    s = zeros(1);
    conStruct_all = struct('resp',s4,'resp20',s4,'fit',s4,'C50r',s,'Rsq',s,'x0',s4);
    conStruct_all(nCellsTot) = conStruct_all;
    
    parfor iCell=1:nCellsTot
        if ~sum(iCell==goodfit_ind_size_all)
            if sum(iCell==[1 nCellsTot]) % check first and last to reset zeros to blank
                conStruct_all(iCell).resp = [];conStruct_all(iCell).resp20 = [];conStruct_all(iCell).fit = [];conStruct_all(iCell).C50r = [];conStruct_all(iCell).Rsq = [];conStruct_all(iCell).x0 = [];
            end
            continue % do not fit unless goodfit_size
        end

        pS = sizeFits_all(iCell,nCon).prefSize;
        pSind = find(szRng==pS);
        ind20 = find(min(abs(szRng-20))==abs(szRng-20),1);
        for iCon = 1:nCon
            if sizeFits_all(iCell,iCon).Ftest
                conStruct_all(iCell).resp(iCon) = sizeFits_all(iCell,iCon).fitout2(pSind);
                conStruct_all(iCell).resp20(iCon) = sizeFits_all(iCell,iCon).fitout2(ind20);
            else
                conStruct_all(iCell).resp(iCon) = sizeFits_all(iCell,iCon).fitout1(pSind);
                conStruct_all(iCell).resp20(iCon) = sizeFits_all(iCell,iCon).fitout1(ind20);
            end
        end

        % fit at prefSize
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
        i50 = find(abs(fitout - R50) == min(abs(fitout - R50)),1);
        C50 = conRng(i50);
        conStruct_all(iCell).C50r = C50;
        %fprintf('Cell %d/%d fit: BL=%.3f Rmax=%.3f C50=%.3f n=%.2f : Rsq=%.3f C50r=%.3f\n',iCell,nCellsTot,cF(1),cF(2),cF(3),cF(4),R2,C50)
        
        % repeat at 20 deg
        cRi = conStruct_all(iCell).resp20;
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
        conStruct_all(iCell).fit20 = cF;
        conStruct_all(iCell).Rsq20 = R2;
        conStruct_all(iCell).x020 = x0best;
        
        fitout = conModelH(cF,conRng);
        R50 = fitout(1)+(fitout(end)-fitout(1))/2;
        i50 = find(abs(fitout - R50) == min(abs(fitout - R50)),1);
        C50 = conRng(i50);
        conStruct_all(iCell).C50r20 = C50;
        %fprintf('Cell %d/%d fit: BL=%.3f Rmax=%.3f C50=%.3f n=%.2f : Rsq=%.3f C50r=%.3f\n',iCell,nCellsTot,cF(1),cF(2),cF(3),cF(4),R2,C50)
    end
    fprintf('Done, saving...\n')
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', 'conStruct_6con.mat');
    save(filename,'conStruct_all');
end

%% present each exp to examine cells+fits and choose example cells
% show each cell rawdata w/ model1+2 overlaid
if 0 % skip plotting
    close all
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
                suptitle(['Exp: ' num2str(iExp) ', area: ' char(expdata.area{iExp}) ', mouse: ' char(expdata.mouse{iExp}) ', date: ' char(expdata.date{iExp})])
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
        suptitle(['Exp: ' num2str(iExp) ', area: ' char(expdata.area{iExp}) ', mouse: ' char(expdata.mouse{iExp}) ', date: ' char(expdata.date{iExp})])
        set(gcf, 'Position', [0 0 800 1000]);
        fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P\HVAcomparison_6con', ['exp_' num2str(iExp) '_SizeTuneFits' num2str(ifig) '.pdf']);
        print(fn_out,'-dpdf')
    end
end
%% now look at each area
close all
areas_c = categorical({'V1','PM'},{'V1','PM'});

areas = ["V1","PM"];
nCells_area = zeros(size(areas));
nExp_area = nCells_area;
for i = 1:length(areas)
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_all);
    switch areas(i)
        case 'V1'
            sizeFits_V1 = sizeFits_all(ind,:);
        case 'PM'
            sizeFits_PM = sizeFits_all(ind,:);
    end
    nCells_area(i) = length(ind);
    nExp_area(i) = length(expIndi);
end
% first examine cell counts, proportion of each model
figure(15);clf;
modelcounts = [sum(~[sizeFits_V1.Ftest]) sum([sizeFits_V1.Ftest]); ...
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
ylabel('Frac. [cells x cons]')
legend('model1','model2','Location','ne')


%% extract unique cells from full set to compare areas
% making figs 1-11
fprintf('Examine cells from each area:\n')
areas = ["V1","PM"];

nExp_area = zeros(size(areas));
nCells_area = nExp_area;

%close all
choosefig = [3 5 8 9 10 11];
% choose figs: 1=modelcounts; 2=averagecurves; 3=prefSize; 4=suppInd;
% 5=conresp; 6=ex.cells; 7=C50f vs C50r; 8=conresp matched size @20deg,
% 9=prefSize but PS within 10-30, 10=conresp but PS within 10-30
% 11= conResp @20 and PS within 10-30deg
legStrs = strings(1,length(areas)); legStrs2=legStrs; legStrsPC=legStrs;
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
            excells = [515 529]; %m1:515, 1153; m2:374
        case 'PM'
            cutoff = 20; %pm cutoff at 20
            excells = [1029 652]; %m1: 590; m2:577, 591
    end
    ind = intersect(ind,find(cellDists_all<cutoff));
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    nExp_area(i) = nExpi;
    nCells_area(i) = nCellsi;
    sizeMean = sizeMean_all(:,:,ind); % (size,con,cell)
    sizeSEM = sizeSEM_all(:,:,ind);
    sizeFits = sizeFits_all(ind,:); %cell,con
    
    prefSize = reshape([sizeFits.prefSize],size(sizeFits));
    suppInd = reshape([sizeFits.suppInd],size(sizeFits));
    suppInd(suppInd<0)=0;suppInd(suppInd>1)=1;
    
    cons_c = categorical({'0.05' '0.1' '0.2' '0.4' '0.8' '1'});
    conStruct = conStruct_all(ind);
    
    legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nCellsi,nExpi);
    
    prefCut = find((prefSize(:,nCon)>10).*(prefSize(:,nCon)<30));
    nPC = length(prefCut);
    legStrsPC(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nPC,nExpi);
    
    if sum(choosefig==1)
        ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
        ism2 = reshape([sizeFits.Ftest],size(sizeFits));
        
        figure(1);if i==1;clf;end %figure 1 = proportions of model2 by con
        subplot(1,2,i)
        modelcounts = [sum(ism1); sum(ism2)]'/nCellsi;
        bar(cons_c,modelcounts,'stacked')
        title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Contrast')
        ylabel('Frac. cells')
        if i==2;legend('m1','m2','location','best');end
    end
    
    if sum(choosefig==2)
        figure(2);if i==1;clf;end %figure 2 = average size tuning curves (normalized)
        % change now to normalize all cells, then plot 3 subplots of m1/m2/all
        szs = 5*1.5.^(0:7); nSz=length(szs);
        cons = 0.05*2.^(0:5); cons(6)=1; nCon=length(cons);
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
        subplot(1,2,i)
        for iCon = 1:nCon
            errorbar(szs,sizeMean_normall(:,iCon),sizeSEM_normall(:,iCon))
            hold on
        end
        title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Size (deg)')
        ylabel('dF/F (norm)')
        ylim([0 1.25])
        if i==2;legend(num2str(cons'),'location','best');end
    end
    
    if sum(choosefig==3)
        figure(3);if i==1;clf;end %figure 3 = prefSize vs con
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
        if i==2;legend(legStrs,'location','southoutside','Orientation','horizontal');end %'location','eastoutside'
    end
    
    if sum(choosefig==4)
        figure(4);if i==1;clf;end %figure 4 = suppInd vs con
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
        legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nCellsi,nExpi);
        if i==2;legend(legStrs,'location','southoutside','Orientation','horizontal');end %'location','eastoutside'
    end
    
    if sum(choosefig==5) %figure 5: average contrast response in each area
        conRng = 0.001:0.001:1;
        opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
        cut = find([conStruct.Rsq]>0.9);
        legStrs2(i)=sprintf('%s (n=%d)',areas(i),length(cut));
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
        title('Mean contrast response @prefSize')
        xlabel('Contrast')
        ylabel('norm. dF/F @ pref size')
        xlim([0 1])
        ylim([0 1.2])
        if i==2;legend(legStrs2,'location','southoutside','Orientation','horizontal');end %'location','southoutside','Orientation','horizontal' for bottom
    
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
        i50 = find(abs(fitout - R50) == min(abs(fitout - R50)),1);
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
    
    if sum(choosefig==6) %figure 6: example cells from each area, with fits
        figure(6);if i==1;clf;end
        subplot(2,2,2*(i-1)+1)
        dum = sizeMean_all(:,:,excells(1)); % take all sizeMean values for cell
        %dum = sizeMean(:,nCon,iCell); % only at highest con
        norm = max(dum(:)); % take max of all dF/F's including all cons
        sizeMean_norm = sizeMean_all(:,:,excells(1))/norm; % normalize by this max for the individual cell
        sizeSEM_norm = sizeSEM_all(:,:,excells(1))/norm;
        for iCon = 1:nCon
            errorbar(szs,sizeMean_norm(:,iCon),sizeSEM_norm(:,iCon))
            hold on
        end
        title(sprintf('Non-suppressed cell in %s',areas(i)))
        if sum(i==[3 4]);xlabel('Size (deg)');end
        if sum(i==[1 3]);ylabel('dF/F (norm)');end
        ylim([0 1.2])
        subplot(2,2,2*(i-1)+2)
        dum = sizeMean_all(:,:,excells(2)); % take all sizeMean values for cell
        %dum = sizeMean(:,nCon,iCell); % only at highest con
        norm = max(dum(:)); % take max of all dF/F's including all cons
        sizeMean_norm = sizeMean_all(:,:,excells(2))/norm; % normalize by this max for the individual cell
        sizeSEM_norm = sizeSEM_all(:,:,excells(2))/norm;
        for iCon = 1:nCon
            errorbar(szs,sizeMean_norm(:,iCon),sizeSEM_norm(:,iCon))
            hold on
        end
        title(sprintf('Suppressed cell in %s',areas(i)))
        if sum(i==[3 4]);xlabel('Size (deg)');end
        %ylabel('dF/F (norm)')
        ylim([0 1.2])
        if i==2;legend(num2str(cons'));end
    end
    
    if sum(choosefig==7) % contrast-resp fits
        % look at different data:
        % histogram of iC50 (best guess) vs C50fit and C50 rec
        % cross plots of these, C50rec-C50fit
        % Rsq histogram, choose cutoff, look at ratio by area
        % C50rec across areas (boxplot + mean)
        C50f = 0*ind;
        for iCell=1:length(ind)
            C50f(iCell) = conStruct(iCell).fit(3);
        end
        C50r = [conStruct.C50r];
        Rsq = [conStruct.Rsq];
        cut = find(Rsq>0.9);
        figure(7);if i==1;clf;end
        subplot(1,2,i)
        plot(C50f,C50r,'.')
        xlabel('C50f')
        ylabel('C50r')
        title({sprintf('Area:%s',areas(i));['(n=' num2str(length(cut)) ', n_{exp}=' num2str(nExpi) ')']})
    end
    
    if sum(choosefig==8) %figure 5: average contrast response at 20 deg in each area
        conRng = 0.001:0.001:1;
        opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
        cut = find([conStruct.Rsq20]>0.9);
        legStrs2(i)=sprintf('%s (n=%d)',areas(i),length(cut));
        conResp = reshape([conStruct(cut).resp20],nCon,length(cut))';
        conResp_norm = conResp./conResp(:,nCon);
        conMean = mean(conResp_norm,1);
        conSEM = std(conResp_norm,[],1)./sqrt(length(cut));
        figure(8);if i==1;clf;end
        ax = gca;
        %subplot(1,2,i)
        %for iCell = 1:nCellsi
        %   p1 = plot(cons,conResp_norm(iCell,:),'r-');
        %   p1.Color(4) = 0.1;
        %   hold on
        %end
        hold on
        ax.ColorOrderIndex = i;
        errorbar(cons,conMean,conSEM)
        %title({sprintf('Contrast response - Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        title('Mean contrast response @20deg')
        xlabel('Contrast')
        ylabel('norm. dF/F @20deg')
        xlim([0 1])
        ylim([0 1.2])
        if i==2;legend(legStrs2,'location','southoutside','Orientation','horizontal');end %'location','southoutside','Orientation','horizontal' for bottom
    
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
        i50 = find(abs(fitout - R50) == min(abs(fitout - R50)),1);
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
    
    if sum(choosefig==9)
        figure(9);if i==1;clf;end %figure 9 = prefSize vs con, cut 10-30
        
        prefMean=zeros(1,nCon);prefSEM=prefMean;
        for iCon=1:nCon
            prefMean(iCon) = mean(prefSize(prefCut,iCon));
            prefSEM(iCon) = std(prefSize(prefCut,iCon))./sqrt(nPC);
        end
        errorbar(cons,prefMean,prefSEM);
        hold on
        title('Mean Preferred Size by Area, only 10-30 deg prefSize')
        xlabel('Contrast')
        ylabel('PrefSize')
        xlim([0 1])
        ylim([0 60])
        if i==2;legend(legStrsPC,'location','southoutside','Orientation','horizontal');end %'location','eastoutside' %'location','southoutside','Orientation','horizontal' for bottom
    end
    
    if sum(choosefig==10) %figure 5: average contrast response in each area, only 10-30 deg prefSize
        conRng = 0.001:0.001:1;
        opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
        cut = find([conStruct.Rsq]>0.9);
        cut = intersect(cut,prefCut);
        legStrs2(i)=sprintf('%s (n=%d)',areas(i),length(cut));
        conResp = reshape([conStruct(cut).resp],nCon,length(cut))';
        conResp_norm = conResp./conResp(:,nCon);
        conMean = mean(conResp_norm,1);
        conSEM = std(conResp_norm,[],1)./sqrt(length(cut));
        figure(10);if i==1;clf;end
        ax = gca;
        %subplot(1,2,i)
        %for iCell = 1:nCellsi
        %   p1 = plot(cons,conResp_norm(iCell,:),'r-');
        %   p1.Color(4) = 0.1;
        %   hold on
        %end
        hold on
        ax.ColorOrderIndex = i;
        errorbar(cons,conMean,conSEM)
        %title({sprintf('Contrast response - Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        title('Mean contrast response @prefSize, only 10-30 deg prefSize')
        xlabel('Contrast')
        ylabel('norm. dF/F @ pref size')
        xlim([0 1])
        ylim([0 1.2])
        if i==2;legend(legStrs2,'location','southoutside','Orientation','horizontal');end %'location','southoutside','Orientation','horizontal' for bottom
    
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
        i50 = find(abs(fitout - R50) == min(abs(fitout - R50)),1);
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
    
    if sum(choosefig==11) %figure 5: average contrast response at 20 deg in each area
        conRng = 0.001:0.001:1;
        opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
        cut = find([conStruct.Rsq20]>0.9);
        cut = intersect(cut,prefCut);
        legStrs2(i)=sprintf('%s (n=%d)',areas(i),length(cut));
        conResp = reshape([conStruct(cut).resp20],nCon,length(cut))';
        conResp_norm = conResp./conResp(:,nCon);
        conMean = mean(conResp_norm,1);
        conSEM = std(conResp_norm,[],1)./sqrt(length(cut));
        figure(11);if i==1;clf;end
        ax = gca;
        %subplot(1,2,i)
        %for iCell = 1:nCellsi
        %   p1 = plot(cons,conResp_norm(iCell,:),'r-');
        %   p1.Color(4) = 0.1;
        %   hold on
        %end
        hold on
        ax.ColorOrderIndex = i;
        errorbar(cons,conMean,conSEM)
        %title({sprintf('Contrast response - Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        title('Mean contrast response @20deg, only 10-30 deg prefSize')
        xlabel('Contrast')
        ylabel('norm. dF/F @20deg')
        xlim([0 1])
        ylim([0 1.2])
        if i==2;legend(legStrs2,'location','southoutside','Orientation','horizontal');end %'location','southoutside','Orientation','horizontal' for bottom
    
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
        i50 = find(abs(fitout - R50) == min(abs(fitout - R50)),1);
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
end

%% fig x - contrast response on 8 random example cells
while 0 %collapse
    
    fprintf('Contrast-response in 8 random examples\n')
    % rn looks at random 18 cells in each area, one area at a time
    
    cons_c = categorical({'0.1' '0.2' '0.4' '0.8'});
    for i = 1:length(areas)
        fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
        % select exps matching area
        expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
        % find cells with correct exp inds, take only good fit cells
        ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_all);
        
        % cutoff by cellDist
        switch areas(i)
            case 'V1'
                cutoff = 10; %v1 cutoff at 10
            case {'LM','AL'}
                cutoff = 15; %alm cutoff at 15
            case 'PM'
                cutoff = 20; %pm cutoff at 20
        end
        ind = intersect(ind,find(cellDists_all<cutoff));
        
        nExpi = length(expIndi);
        nCellsi = length(ind);
        sizeTune = sizeTune_all{:,:,ind}; % (size,con,cell)
        sizeMean = sizeMean_all(:,:,ind);
        sizeSEM = sizeSEM_all(:,:,ind);
        sizeFits = sizeFits_all(ind,:); %cell,con
        lbub_fits = lbub_fits_all(ind,:,:); %cell,par,val (low up mean true stdev)
        ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
        ism2 = reshape([sizeFits.Ftest],size(sizeFits));
        
        
        sizeMean_norm = sizeMean*0; sizeSEM_norm = sizeSEM*0;
        for iCell = 1:nCellsi
            dum = sizeMean(:,:,iCell); % take all sizeMean values for cell
            %dum = sizeMean(:,nCon,iCell); % only at highest con
            norm = max(dum(:)); % take max of all dF/F's including all cons
            sizeMean_norm(:,:,iCell) = sizeMean(:,:,iCell)/norm; % normalize by this max for the individual cell
            sizeSEM_norm(:,:,iCell) = sizeSEM(:,:,iCell)/norm;
        end
        % contrast dependence at prefSize(0.8)
        % first plot example cells (random 16 from each area)
        % and plot size tuning curves at all cons with vertical line at prefSize(0.8)
        % next to plot of dF/F resp vs contrast at peak size
        figure(5);clf;
        count = 1;
        szRng = linspace(0,max(szs));
        nRand = 1;
        [n1, n2] = subplotn(2);
        cellshufs = randi(nCellsi,nRand,1);
        for iRand=1:nRand
            iCell = cellshufs(iRand);
            pS = sizeFits(iCell,nCon).prefSize;
            pSind = find(szRng==pS);
            subplot(n2,n1,count)%6,6,count)
            for iCon = 1:nCon
                errorbar(szs,sizeMean(:,iCon,iCell),sizeSEM(:,iCon,iCell))
                hold on
            end
            ax = gca;
            ax.ColorOrderIndex = 1;
            for iCon = 1:nCon
                if sizeFits(iCell,iCon).Ftest
                    plot(szRng,sizeFits(iCell,iCon).fitout2)
                else
                    plot(szRng,sizeFits(iCell,iCon).fitout1)
                end
                hold on
            end
            line([pS pS],[-1,2],'Color','red','LineStyle','--')
            if sizeFits(iCell,iCon).Ftest
                title([sprintf('Area:%s, cell# %d',areas(i),iCell) ' m2'])
            else
                title([sprintf('Area:%s, cell# %d',areas(i),iCell) ' m1'])
            end
            xlabel('Size (deg)')
            ylabel('dF/F')
            maxAmp = max([sizeFits(iCell,:).maxResp1 sizeFits(iCell,:).maxResp2]);
            ylim([-0.2*maxAmp 1.2*maxAmp])
            if iRand==nRand;legend(num2str(cons'));end
            count=count+1;
            subplot(n2,n1,count)%6,6,count)
            conResp = zeros(nCon,1);
            for iCon = 1:nCon
                if sizeFits(iCell,iCon).Ftest
                    conResp(iCon) = sizeFits(iCell,iCon).fitout2(pSind);
                else
                    conResp(iCon) = sizeFits(iCell,iCon).fitout1(pSind);
                end
            end
            plot(cons,conResp)
            title(['Cell ' num2str(iCell) ' resp vs con @ ' num2str(pS,3)])
            xlabel('Contrast')
            ylabel('dF/F')
            ylim([0 1.2*max(conResp)])
            count=count+1;
        end
        pause
        
        if sum(choosefig==6)
            conResp = conResp_all(ind,:);
            conResp_norm = conResp./conResp(:,nCon);
            ism1_n = find(ism1(:,nCon));
            ism2_n = find(ism2(:,nCon));
            conMean1 = mean(conResp_norm(ism1_n,:),1);
            conSEM1 = std(conResp_norm(ism1_n,:),[],1)./sqrt(length(ism2_n));
            conMean2 = mean(conResp_norm(ism2_n,:),1);
            conSEM2 = std(conResp_norm(ism2_n,:),[],1)./sqrt(length(ism2_n));
            figure(6);if i==1;clf;end %figure 6: example cells from each area, with fits
            for iCell = 1:nCellsi
                if sizeFits(iCell,nCon).Ftest
                    subplot(2,4,2*(i-1)+2)
                else
                    subplot(2,4,2*(i-1)+1)
                end
                p1 = plot(cons,conResp_norm(iCell,:),'r-');
                p1.Color(4) = 0.1;
                hold on
            end
            subplot(2,4,2*(i-1)+1)
            errorbar(cons,conMean1,conSEM1,'k')
            title({sprintf('ConResp M1 - Area:%s',areas(i));['(n=' num2str(length(ism1_n)) ', n_{exp}=' num2str(nExpi) ')']})
            xlabel('Contrast')
            ylabel('norm. dF/F @ pref size')
            xlim([0 1])
            ylim([0 1.5])
            subplot(2,4,2*(i-1)+2)
            errorbar(cons,conMean2,conSEM2,'k')
            title({sprintf('ConResp M2 - Area:%s',areas(i));['(n=' num2str(length(ism2_n)) ', n_{exp}=' num2str(nExpi) ')']})
            xlabel('Contrast')
            ylabel('norm. dF/F @ pref size')
            xlim([0 1])
            ylim([0 1.5])
        end
    end
end

%% plot example cell at highest con with models overlaid (fig. 5c)
iCell = 2054; %randi(sum(nCellsExp),1);
cellFit = sizeFits_all(iCell,nCon);
figure(8);clf;
plot(cellFit.szs0,cellFit.data,'k.')
hold on
plot(szRng,cellFit.fitout1,'r');
plot(szRng,cellFit.fitout2,'b');
errorbar(szs,sizeMean_all(:,nCon,iCell),sizeSEM_all(:,nCon,iCell),'k')
title(['Cell #' num2str(iCell)])
xlabel('Size (deg)')
ylabel('dF/F')
legend('data','model1','model2')
fprintf(['\nSSE1: ' num2str(cellFit.SSE1) ', SSE2: ' num2str(cellFit.SSE2) '\nF-score: ' num2str(cellFit.Fscore) ', F-crit(' num2str(length(cellFit.data)) ' pts): ' num2str(cellFit.Fcrit)])


%% plot 2 chosen cells in V1 and PM (fig. 8a)
cons_c = categorical({'0.1' '0.2' '0.4' '0.8'});
conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));
conRng = 0.001:0.001:1;
szRng = linspace(0,max(szs));

areas2 = ["V1" "PM"];
%chosen = [randi(339,1) randi(73,1)];
chosen = [1 3];
for i = 1:length(areas2)
    fprintf(['Area #' num2str(i) ' : ' char(areas2(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas2(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_all);
    
    % cutoff by cellDist
    switch areas(i)
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
        case {'LM','AL'}
            cutoff = 15; %alm cutoff at 15
        case 'PM'
            cutoff = 20; %pm cutoff at 20
    end
    ind = intersect(ind,find(cellDists_all<cutoff));
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    sizeTune = sizeTune_all{:,:,ind}; % (size,con,cell)
    sizeMean = sizeMean_all(:,:,ind);
    sizeSEM = sizeSEM_all(:,:,ind);
    sizeFits = sizeFits_all(ind,:); %cell,con
    lbub_fits = lbub_fits_all(ind,:,:); %cell,par,val (low up mean true stdev)
    ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
    ism2 = reshape([sizeFits.Ftest],size(sizeFits));
    
    conStruct = conStruct_all(ind);
    
    % contrast dependence at prefSize(0.8)
    % first plot example cells (random 16 from each area)
    % and plot size tuning curves at all cons with vertical line at prefSize(0.8)
    % next to plot of dF/F resp vs contrast at peak size
    iCell = chosen(i)
    figure(5);if i==1;clf;end
    pS = sizeFits(iCell,nCon).prefSize;
    pSind = find(szRng==pS);
    subplot(2,2,2*(i-1)+1)
    for iCon = 1:nCon
        errorbar(szs,sizeMean(:,iCon,iCell),sizeSEM(:,iCon,iCell))
        hold on
    end
    ax = gca;
    ax.ColorOrderIndex = 1;
    for iCon = 1:nCon
        if sizeFits(iCell,iCon).Ftest
            plot(szRng,sizeFits(iCell,iCon).fitout2)
        else
            plot(szRng,sizeFits(iCell,iCon).fitout1)
        end
        hold on
    end
    line([pS pS],[-1,2],'Color','red','LineStyle','--')
    if sizeFits(iCell,iCon).Ftest
        title([sprintf('Area:%s, cell# %d',areas2(i),iCell) ' m2'])
    else
        title([sprintf('Area:%s, cell# %d',areas2(i),iCell) ' m1'])
    end
    xlabel('Size (deg)')
    ylabel('dF/F')
    maxAmp = max([sizeFits(iCell,:).maxResp1 sizeFits(iCell,:).maxResp2]);
    ylim([-0.2*maxAmp 1.2*maxAmp])
    if i==2;legend(num2str(cons'),'location','best');end
    subplot(2,2,2*(i-1)+2)
    plot(cons,conStruct(iCell).resp,'ko')
    hold on
    plot(conRng,conModelH(conStruct(iCell).fit,conRng),'-r')
    C50 = conStruct(iCell).C50r;
    R50 = conModelH(conStruct(iCell).fit,C50);
    plot(C50,R50,'bx')
    line([C50 C50],[0 R50],'color','b')
    title(['Cell ' num2str(iCell) ' resp vs con @ ' num2str(pS,3)])
    xlabel('Contrast')
    ylabel('dF/F')
    ylim([0 1.2*max(conStruct(iCell).resp)])
    if i==2;legend('data','fit','C50','Location','best');end
end

%% examine contrast fits - C50 and BLr
Rsqcutoff = 0.9;
% extract C50i and
C50i = 0*goodfit_ind_size_all;
C50f = 0*goodfit_ind_size_all;
for i=1:length(goodfit_ind_size_all)
    C50i(i) = conStruct_all(goodfit_ind_size_all(i)).x0(3);
    C50f(i) = conStruct_all(goodfit_ind_size_all(i)).fit(3);
end
C50r = [conStruct_all.C50r];

Rsq = [conStruct_all.Rsq];
figure(9);clf;
subplot(1,2,1)
histogram(Rsq,-0.05:0.1:1.05)
title('R^2 all cells')
cut = find(Rsq>Rsqcutoff);
subplot(1,2,2)
histogram(Rsq(cut),(Rsqcutoff-0.01):0.01:1.01)
title(['R^2 cutoff at >' num2str(Rsqcutoff,2)])

nCut = 0*nCells_area;
rCut=nCut;
x =[]; yC50=x; yBLr=x;
for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_all);
    
    % cutoff by cellDist
    switch areas(i)
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
        case {'LM','AL'}
            cutoff = 15; %alm cutoff at 15
        case 'PM'
            cutoff = 20; %pm cutoff at 20
    end
    ind = intersect(ind,find(cellDists_all<cutoff));
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    sizeFits = sizeFits_all(ind,:); %cell,con
    lbub_fits = lbub_fits_all(ind,:,:); %cell,par,val (low up mean true stdev)
    ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
    ism2 = reshape([sizeFits.Ftest],size(sizeFits));
    
    conStruct = conStruct_all(ind);
    C50f = 0*ind;
    baseline = C50f;
    Rmax = C50f;
    for iCell = 1:nCellsi
        C50f(iCell) = conStruct(iCell).fit(3);
        baseline(iCell) = conStruct(iCell).fit(1);
        Rmax(iCell) = conStruct(iCell).fit(2);
    end
    C50r = [conStruct.C50r]';
    
    Rsq = [conStruct.Rsq];
    cut = find(Rsq>Rsqcutoff);
    nCut(i) = length(cut);
    rCut(i) = nCut(i)/nCellsi;
    
    x = [x; i*ones(size(C50r(cut)))];
    yC50 = [yC50; C50r(cut)];
    C50mean(i) = mean(C50r(cut));
    C50SEM(i) = std(C50r(cut))./sqrt(length(cut));
    BLr = baseline(cut)./(baseline(cut)+Rmax(cut));
    yBLr = [yBLr; BLr];
    BLrmean(i) = mean(BLr);
    BLrSEM(i) = std(BLr)./sqrt(length(cut));
    
    figure(11);if i==1;clf;end
    subplot(2,2,i)
    histogram(C50r,0.15:0.05:0.85)
    hold on
    histogram(C50r(cut),0.15:0.05:0.85)
    title(['C50r, area:' char(areas(i))])
    xlabel('C50_r')
    if i==4;legend('all',['R^2>' num2str(Rsqcutoff,2)]);end
    
    %figure(12);if i==1;clf;end
    %subplot(2,2,i)
    %histogram(C50r(cut)./C50f(cut),0.15:0.1:1.05)
    %title(['Ratio rC_{50}/fC_{50}, area:' char(areas(i))])
    %xlabel('rC50/fC50')
    %ylabel('#cells')
    
    %figure(13);if i==1;clf;end
    %subplot(2,2,i)
    %plot(Rsq(cut),C50r(cut)./C50f(cut),'.')
    %title(['Ratio rC_{50}/fC_{50} vs Rsq, area:' char(areas(i))])
    %xlabel('R^2')
    %ylabel('rC50/fC50')
    
    figure(14);if i==1;clf;end
    subplot(2,2,i)
    histogram(baseline(cut)./(baseline(cut)+Rmax(cut)),-0.05:0.1:1.05)
    title(['BL/(BL+R_{max}), area:' char(areas(i))])
    xlabel('BL/(BL+R_{max})')
    ylabel('#cells')
    
end

figure(15);clf
c_areas = categorical({'V1' 'PM' 'all'},{'V1' 'PM' 'all'});
nTot = sum(nCut./rCut);
nCut2 = [nCut nTot];
rCut2 = [rCut sum(nCut)/nTot];
bar(c_areas,rCut2)
ylim([0 1])
xlabel('Area')
ylabel('frac. cutoff')
title('Con-fit cutoffs in each area')

figure(16);clf;
boxplot(yC50,x)
hold on
errorbar(1:2,C50mean,C50SEM,'x')
plot([1 2],[0.82 0.82],'-k', 'LineWidth',2)
plot([1.45 1.55],[0.85 0.85],'*k')
hold off
set(gca,'XTick',1:2,'XTickLabel',{['V1 (n=' num2str(nCut(1)) ')'],['PM (n=' num2str(nCut(2)) ')']})
title('C_{50} by area')
ylabel('C_{50}')
legend('mean+SEM')
ylim([0 1])

figure(17);clf;
boxplot(yBLr,x)
hold on
errorbar(1:2,BLrmean,BLrSEM,'x')
plot([1 2],[0.83 0.83],'-k', 'LineWidth',2)
plot([1.4 1.5 1.6],[0.85 0.85 0.85],'*k')
hold off
set(gca,'XTick',1:2,'XTickLabel',{['V1 (n=' num2str(nCut(1)) ')'],['PM (n=' num2str(nCut(2)) ')']})
title('BL/(BL+R_{max}) by area')
ylabel('BL/(BL+R_{max})')
legend('mean+SEM')
ylim([-0.1 1])

%% stats on C50, BLr
[p,~,statC50] = anova1(yC50,x)
[results, means] = multcompare(statC50,'CType','hsd')

[p,~,statBLr] = anova1(yBLr,x)
[results, means] = multcompare(statBLr,'CType','hsd')