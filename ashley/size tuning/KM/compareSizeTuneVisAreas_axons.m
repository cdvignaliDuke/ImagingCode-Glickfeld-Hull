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
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = expdata.run_str{i};
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_sizeTuneData.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'sizeTune', 'sizeMean', 'sizeSEM', 'cellDists')
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']);
    load(filename);
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_SP.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_sizeFitResults_SP.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'sizeFits')
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind_size')
    
    sizeTune_all = cat(3,sizeTune_all,sizeTune(:,:));
    sizeMean_all = cat(3,sizeMean_all,sizeMean(:,:));
    sizeSEM_all = cat(3,sizeSEM_all,sizeSEM(:,:));
    cellDists_all = [cellDists_all;cellDists];
    nCellsExp(i) = length(cellDists);
    sizeFits_all = cat(1,sizeFits_all,sizeFits);
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

nCellsTot_6con = sum(nCellsExp_6con);
fprintf('6 con: %d cells loaded (%d goodfit_size)\n',nCellsTot_6con,length(goodfit_ind_size_6con))
expInd_6con = [];
expArea_6con = [];
for i = 1:nExp
    expInd_6con = [expInd_6con repmat(i,1,nCellsExp_6con(i))];
    expArea_6con = [expArea_6con repmat(expdata(i,4).area, 1,nCellsExp_6con(i))];
end

szs = unique(celleqel2mat_padded(input.tGratingDiameterDeg)); 
nSize = length(szs);
szRng = linspace(0,max(szs));

%% con resp for 4 con data
fprintf('Extracting contrast response of each cell at prefSize\n')
cons = con4;
nCon = length(cons);
s4 = zeros(1,nCon);
s = zeros(1);
conStruct_all = struct('resp',s4,'fit',s4,'C50r',s,'Rsq',s,'x0',s4);
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
        pS = sizeFits_all(iCell,nCon).prefSize;
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
        if sizeFits_all(iCell,iCon).Ftest
            conStruct_all(iCell).resp(iCon) = sizeFits_all(iCell,iCon).fitout2(pSind);
        else
            conStruct_all(iCell).resp(iCon) = sizeFits_all(iCell,iCon).fitout1(pSind);
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
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', 'conStruct23.mat');
save(filename,'conStruct_all');

%% con resp for 6 con data
fprintf('Extracting contrast response of each cell at prefSize\n')
cons = con6;
nCon = length(cons);
s4 = zeros(1,nCon);
s = zeros(1);
conStruct_6con = struct('resp',s4,'fit',s4,'C50r',s,'Rsq',s,'x0',s4);
conStruct_allsize_6con = struct('resp',s4,'fit',s4,'C50r',s,'Rsq',s,'x0',s4);
conStruct_6con(nCellsTot_6con) = conStruct_6con;
conStruct_allsize_6con(nCellsTot_6con,nSize) = conStruct_allsize_6con;
conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));
conRng = 0:0.001:1;
opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
goodSize_6con = zeros(nCellsTot_6con,nSize);
usePrefSize = 0;
useRFSize = 1;
doAllSize = 0;

for iCell=1:nCellsTot_6con
    if ~sum(iCell==goodfit_ind_size_6con)
        if sum(iCell==[1 nCellsTot_6con]) % check first and last to reset zeros to blank
            conStruct_6con(iCell).resp = [];conStruct_6con(iCell).fit = [];conStruct_6con(iCell).C50r = [];conStruct_6con(iCell).Rsq = [];conStruct_6con(iCell).x0 = [];
        end
        continue % do not fit unless goodfit_size
    end
    
    if usePrefSize
        pS = sizeFits_6con(iCell,nCon).prefSize;
        pSind = find(szRng==pS);
    elseif useRFSize
        switch cell2mat(expArea_6con(iCell))
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
        if sizeFits_6con(iCell,iCon).Ftest
            conStruct_6con(iCell).resp(iCon) = sizeFits_6con(iCell,iCon).fitout2(pSind);
        else
            conStruct_6con(iCell).resp(iCon) = sizeFits_6con(iCell,iCon).fitout1(pSind);
        end 
    end
    
    cRi = conStruct_6con(iCell).resp;
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
    
    conStruct_6con(iCell).fit = cF;
    conStruct_6con(iCell).Rsq = R2;
    conStruct_6con(iCell).x0 = x0best;
    
    fitout = conModelH(cF,conRng);
    R50 = fitout(1)+(fitout(end)-fitout(1))/2;
    fitout50rect = abs(fitout - R50);
    i50 = find(fitout50rect == min(fitout50rect),1);
    C50 = conRng(i50);
    conStruct_6con(iCell).C50r = C50;
    
    if doAllSize   
        for isz = 1:length(szs)
            szData(iCell,isz).data = [];
            szData(iCell,isz).cons0 = [];
            for iCon = 1:nCon
                ind = find(sizeFits_6con(iCell,iCon).szs0 == szs(isz));
                szData(iCell,isz).data = [szData(iCell,isz).data sizeFits_6con(iCell,iCon).data(ind)];
                szData(iCell,isz).cons0 = [szData(iCell,isz).cons0 cons(iCon).*ones(1,length(ind))];
            end
            p(isz) = anova1(szData(iCell,isz).data, szData(iCell,isz).cons0,'off');
        end
        h = find(p<(0.05));
        goodSize_6con(iCell,h) = 1;
        for isz = 1:nSize
            if find(h == isz)
                for iCon = 1:nCon
                    ind = find(szData(iCell,isz).cons0 == cons(iCon));
                    conStruct_allsize_6con(iCell,isz).resp(iCon) = mean(szData(iCell,isz).data(ind),2);
                end
                cRi = conStruct_allsize_6con(iCell,isz).resp;
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

                conStruct_allsize_6con(iCell,isz).fit = cF;
                conStruct_allsize_6con(iCell,isz).Rsq = R2;
                conStruct_allsize_6con(iCell,isz).x0 = x0best;

                fitout = conModelH(cF,conRng);
                R50 = fitout(1)+(fitout(end)-fitout(1))/2;
                fitout50rect = abs(fitout - R50);
                i50 = find(fitout50rect == min(fitout50rect),1);
                C50 = conRng(i50);
                conStruct_allsize_6con(iCell,isz).C50r = C50;
            else
                conStruct_allsize_6con(iCell,isz).resp = NaN;
                conStruct_allsize_6con(iCell,isz).fit = NaN;
                conStruct_allsize_6con(iCell,isz).Rsq = NaN;
                conStruct_allsize_6con(iCell,isz).x0 = NaN;
                conStruct_allsize_6con(iCell,isz).C50r = NaN;
            end
        end
    end
    
    fprintf('Cell %d/%d fit: BL=%.3f Rmax=%.3f C50=%.3f n=%.2f : Rsq=%.3f C50r=%.3f\n',iCell,nCellsTot_6con,cF(1),cF(2),cF(3),cF(4),R2,C50)
    
end
fprintf('Done, saving...\n')
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', 'conStruct23_6con.mat');
save(filename,'conStruct_6con', 'conStruct_allsize_6con');

%% load contrast response instead of compute
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', 'conStruct23.mat');
load(filename);
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', 'conStruct23_6con.mat');
load(filename);
conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));

%% present each exp to examine cells+fits and choose example cells
% show each cell rawdata w/ model1+2 overlaid

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
    fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P\HVAcomparison', ['exp_' num2str(iExp) '_SizeTuneFits' num2str(ifig) '.pdf']);
    print(fn_out,'-dpdf')
end

%% compare areas for 4 con data
% making figs 1-7
fprintf('Examine cells from each area:\n')
areas = ["V1","LM","AL","PM"];
nArea = length(areas);
cons = con4;
nCon = length(cons);

nExp_area = zeros(size(areas));
nCells_area = nExp_area;

close all
choosefig = [1:5];
% choose figs: 1=modelcounts; 2=averagecurves; 3=prefSize; 4=suppInd; 5=conresp; 6=ex.cells; 7=medianfits;
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
            cutoff = 15; %lm cutoff at 15
            excells = [1861 1863];
        case 'AL'
            cutoff = 15; %al cutoff at 15
            excells = [2395 1777];
        case 'PM'
            cutoff = 20; %pm cutoff at 20
            excells = [1952 2292];
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
    
    cons_c = categorical(cellstr(num2str(cons'))');
    szs_c = cellstr(num2str(chop(szs,2)'))';
    conStruct = conStruct_all(ind);
    
    legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nCellsi,nExpi);
    
    if sum(choosefig==1)
        figure(1);if i==1;clf;end %figure 1 = proportions of model2 by con
        subplot(1,nArea,i)
        modelcounts = [sum(ism1); sum(ism2)]'/nCellsi;
        bar(cons_c,modelcounts,'stacked')
        title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Contrast')
        ylabel('Frac. cells')
        if i==nArea;legend('m1','m2','location','best');end
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
        subplot(1,nArea,i)
        for iCon = 1:nCon
            errorbar(szs,sizeMean_normall(:,iCon),sizeSEM_normall(:,iCon))
            hold on
        end
        title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Size (deg)')
        ylabel('dF/F (norm)')
        ylim([0 1.25])
        if i==nArea;legend(num2str(cons'),'location','best');end
    end
    
    if sum(choosefig==3)
        figure(3);if i==1;clf;end %figure 3 = prefSize vs con
        prefSize = reshape([sizeFits.prefSize],size(sizeFits));
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
        if i==nArea;legend(legStrs,'location','eastoutside');end %'location','southoutside','Orientation','horizontal' for bottom
    end
    
    if sum(choosefig==4)
        figure(4);if i==1;clf;end %figure 4 = suppInd vs con
        suppInd = reshape([sizeFits.suppInd],size(sizeFits));
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
        legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nCellsi,nExpi);
        if i==nArea;legend(legStrs,'location','eastoutside');end %'location','southoutside','Orientation','horizontal' for bottom
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
        title('Mean contrast response by area')
        xlabel('Contrast')
        ylabel('norm. dF/F @ pref size')
        xlim([0 1])
        ylim([0 1.2])
        if i==nArea;legend(legStrs2,'location','southoutside','Orientation','horizontal');end %'location','southoutside','Orientation','horizontal' for bottom
    
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
    
    
    if sum(choosefig==6) %figure 6: example cells from each area, with fits
        figure(6);if i==1;clf;end
        subplot(2,4,2*(i-1)+1)
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
        subplot(2,4,2*(i-1)+2)
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
        if i==4;legend(num2str(cons'));end
    end
    
    if sum(choosefig==7) % median model parameters curve
        par1 = zeros(nCellsi,nCon,3);
        par2 = zeros(nCellsi,nCon,6);
        medpar1 = zeros(nCon,3);
        medpar2 = zeros(nCon,6);
        for iCell = 1:nCellsi
            for iCon = 1:nCon
                par1(iCell,iCon,:) = sizeFits(iCell,iCon).fit1.c1;
                par2(iCell,iCon,:) = sizeFits(iCell,iCon).fit2.c2;
            end
        end
        for iCon = 1:nCon
            ism1_i = find(~[sizeFits(:,iCon).Ftest]);
            ism2_i = find([sizeFits(:,iCon).Ftest]);
            medpar1(iCon,:) = squeeze(median(par1(ism1_i,iCon,:),1));
            medpar2(iCon,:) = squeeze(median(par2(ism2_i,iCon,:),1));
            meanpar1(iCon,:) = squeeze(mean(par1(ism1_i,iCon,:),1));
            meanpar2(iCon,:) = squeeze(mean(par2(ism2_i,iCon,:),1));
        end
        medpar1(:,1) = medpar1(:,1)./max(medpar1(nCon,1));
        medpar2(:,1) = medpar2(:,1)./max(medpar2(nCon,1));
        meanpar1(:,1) = meanpar1(:,1)./max(meanpar1(nCon,1));
        meanpar2(:,1) = meanpar2(:,1)./max(meanpar2(nCon,1));
        medpar1 = (medpar1+meanpar1)/2;
        medpar2 = (medpar2+meanpar2)/2;
        m1 = str2func(sizeFits(1,1).m1);
        m2 = str2func(sizeFits(1,1).m2);
        
        figure(7);if i==1;clf;end
        subplot(2,2,i)
        for iCon=1:nCon
            plot(szRng,m1(medpar1(iCon,:),szRng))
            hold on
        end
        ax = gca;
        ax.ColorOrderIndex = 1;
        for iCon=1:nCon
            plot(szRng,m2(medpar2(iCon,:),szRng))
            hold on
        end
        legend(num2str(cons'))
        title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Size (deg)')
        ylabel('dF/F')
        xlim([0 max(szs)+1])
    end
    
    if sum(choosefig==8) % contrast-resp fits
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
        figure(8);if i==1;clf;end
        subplot(2,2,i)
        plot(C50f,C50r,'.')
        xlabel('C50f')
        ylabel('C50r')
        title({sprintf('Area:%s',areas(i));['(n=' num2str(length(cut)) ', n_{exp}=' num2str(nExpi) ')']})
        
    end
    
end

%% compare areas for 6 con data
% making figs 1-7
fprintf('Examine cells from each area:\n')
areas = ["V1","PM"];
nArea = length(areas);
cons = con6;
nCon = length(cons);

nExp_area = zeros(size(areas));
nCells_area = nExp_area;

close all
choosefig = [5];
% choose figs: 1=modelcounts; 2=averagecurves; 3=prefSize; 4=suppInd; 5=conresp; 6=ex.cells; 7=medianfits;
legStrs = strings(1,length(areas)); legStrs2=legStrs;
for i = 1:nArea
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area with 6 con data
    expIndi = intersect(find(nCellsExp_6con), find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1)));
    if length(expIndi > 0)
        % find cells with correct exp inds, take only good fit cells
        ind = intersect(find(ismember(expInd_6con,expIndi)),goodfit_ind_size_6con);

        % cutoff by cellDist
        % try looking with different cutoffs
        switch areas(i)
            case 'V1'
                cutoff = 10; %v1 cutoff at 10
                excells = [631 2128];
                % case {'LM','AL'}
                %     cutoff = 15; %alm cutoff at 15
            case 'LM'
                cutoff = 15; %alm cutoff at 15
                excells = [1861 1863];
            case 'AL'
                cutoff = 15; %alm cutoff at 15
                excells = [2395 1777];
            case 'PM'
                cutoff = 20; %pm cutoff at 20
                excells = [1952 2292];
        end
        ind = intersect(ind,find(cellDists_6con<cutoff));

        nExpi = length(expIndi);
        nCellsi = length(ind);
        nExp_area(i) = nExpi;
        nCells_area(i) = nCellsi;
        sizeTune = sizeTune_6con{:,:,ind}; % (size,con,cell)
        sizeMean = sizeMean_6con(:,:,ind);
        sizeSEM = sizeSEM_6con(:,:,ind);
        sizeFits = sizeFits_6con(ind,:); %cell,con
        lbub_fits = lbub_fits_6con(ind,:,:); %cell,par,val (low up mean true stdev)
        ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
        ism2 = reshape([sizeFits.Ftest],size(sizeFits));

        cons_c = categorical(cellstr(num2str(cons'))');
        szs_c = cellstr(num2str(chop(szs,2)'))';
        conStruct = conStruct_6con(ind);
        conStruct_allsize = conStruct_allsize_6con(ind,:);

        legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nCellsi,nExpi);

        if sum(choosefig==1)
            figure(1);if i==1;clf;end %figure 1 = proportions of model2 by con
            subplot(1,nArea,i)
            modelcounts = [sum(ism1); sum(ism2)]'/nCellsi;
            bar(cons_c,modelcounts,'stacked')
            title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
            xlabel('Contrast')
            ylabel('Frac. cells')
            if i==nArea;legend('m1','m2','location','best');end
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
            subplot(1,nArea,i)
            for iCon = 1:nCon
                errorbar(szs,sizeMean_normall(:,iCon),sizeSEM_normall(:,iCon))
                hold on
            end
            title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
            xlabel('Size (deg)')
            ylabel('dF/F (norm)')
            ylim([0 1.25])
            if i==nArea;legend(num2str(cons'),'location','best');end
        end

        if sum(choosefig==3)
            figure(3);if i==1;clf;end %figure 3 = prefSize vs con
            prefSize = reshape([sizeFits.prefSize],size(sizeFits));
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
            if i==nArea;legend(legStrs,'location','eastoutside');end %'location','southoutside','Orientation','horizontal' for bottom
        end

        if sum(choosefig==4)
            figure(4);if i==1;clf;end %figure 4 = suppInd vs con
            suppInd = reshape([sizeFits.suppInd],size(sizeFits));
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
            legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nCellsi,nExpi);
            if i==nArea;legend(legStrs,'location','eastoutside');end %'location','southoutside','Orientation','horizontal' for bottom
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
            title('Mean contrast response by area')
            xlabel('Contrast')
            if usePrefSize
                ylabel('norm. dF/F @ pref size')
            elseif useRFSize
                ylabel('norm. dF/F @ RF size')
            end
            xlim([0 1])
            ylim([0 1.2])
            if i==nArea;legend(legStrs2,'location','southoutside','Orientation','horizontal');end

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
        if sum(choosefig==6) %figure 6: fraction of contrast modulated cells by size by area
            figure(6); if i==1;clf;end
            subplot(1,nArea,i)
            n = length(ind);
            fractResp = sum(goodSize_6con(ind,:),1)./n;
            bar(fractResp);
            xticklabels(szs_c); xticks(1:1:nSize);
            title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
            xlabel('Size (Deg)')
            ylabel('Frac. cells')
        end
        if sum(choosefig==7) %figure 7: average c50 by size by area
            figure(7); if i==1;clf;end
            c50_allsize = nan(n,nSize);
            for iC = 1:n
                for isz = 1:nSize
                    if find(isz == find(goodSize_6con(ind(iC),:)))
                        if conStruct_allsize(iC,isz).Rsq>0.9
                        	c50_allsize(iC,isz) = conStruct_allsize(iC,isz).C50r;
                        end
                    end
                end
            end
            errorbar(szs, nanmean(c50_allsize,1), nanstd(c50_allsize,[],1)./sqrt(sum(~isnan(c50_allsize),1)),'-o');
            hold on
            title('Mean C50 by area')
            xlabel('Size (Deg)')
            ylabel('c50')
            xlim([0 100])
            ylim([0 1])
            legStrs(i)=sprintf('%s (n=%s, n_{exp}=%d)',areas(i),num2str(sum(~isnan(c50_allsize),1)),nExpi);
            if i==nArea;legend(legStrs,'location','northeast');end %'location','southoutside','Orientation','horizontal' for bottom
        end
    end
end
            %%
           
        if sum(choosefig==6) %figure 6: example cells from each area, with fits
            figure(6);if i==1;clf;end
            subplot(1,nArea,2*(i-1)+1)
            dum = sizeMean_6con(:,:,excells(1)); % take all sizeMean values for cell
            %dum = sizeMean(:,nCon,iCell); % only at highest con
            norm = max(dum(:)); % take max of all dF/F's including all cons
            sizeMean_norm = sizeMean_6con(:,:,excells(1))/norm; % normalize by this max for the individual cell
            sizeSEM_norm = sizeSEM_6con(:,:,excells(1))/norm;
            for iCon = 1:nCon
                errorbar(szs,sizeMean_norm(:,iCon),sizeSEM_norm(:,iCon))
                hold on
            end
            title(sprintf('Non-suppressed cell in %s',areas(i)))
            if sum(i==[3 4]);xlabel('Size (deg)');end
            if sum(i==[1 3]);ylabel('dF/F (norm)');end
            ylim([0 1.2])
            subplot(2,nArea,2*(i-1)+2)
            dum = sizeMean_6con(:,:,excells(2)); % take all sizeMean values for cell
            %dum = sizeMean(:,nCon,iCell); % only at highest con
            norm = max(dum(:)); % take max of all dF/F's including all cons
            sizeMean_norm = sizeMean_6con(:,:,excells(2))/norm; % normalize by this max for the individual cell
            sizeSEM_norm = sizeSEM_6con(:,:,excells(2))/norm;
            for iCon = 1:nCon
                errorbar(szs,sizeMean_norm(:,iCon),sizeSEM_norm(:,iCon))
                hold on
            end
            title(sprintf('Suppressed cell in %s',areas(i)))
            if sum(i==[3 4]);xlabel('Size (deg)');end
            %ylabel('dF/F (norm)')
            ylim([0 1.2])
            if i==nArea;legend(num2str(cons'));end
        end

        if sum(choosefig==7) % median model parameters curve
            par1 = zeros(nCellsi,nCon,3);
            par2 = zeros(nCellsi,nCon,6);
            medpar1 = zeros(nCon,3);
            medpar2 = zeros(nCon,6);
            for iCell = 1:nCellsi
                for iCon = 1:nCon
                    par1(iCell,iCon,:) = sizeFits(iCell,iCon).fit1.c1;
                    par2(iCell,iCon,:) = sizeFits(iCell,iCon).fit2.c2;
                end
            end
            for iCon = 1:nCon
                ism1_i = find(~[sizeFits(:,iCon).Ftest]);
                ism2_i = find([sizeFits(:,iCon).Ftest]);
                medpar1(iCon,:) = squeeze(median(par1(ism1_i,iCon,:),1));
                medpar2(iCon,:) = squeeze(median(par2(ism2_i,iCon,:),1));
                meanpar1(iCon,:) = squeeze(mean(par1(ism1_i,iCon,:),1));
                meanpar2(iCon,:) = squeeze(mean(par2(ism2_i,iCon,:),1));
            end
            medpar1(:,1) = medpar1(:,1)./max(medpar1(nCon,1));
            medpar2(:,1) = medpar2(:,1)./max(medpar2(nCon,1));
            meanpar1(:,1) = meanpar1(:,1)./max(meanpar1(nCon,1));
            meanpar2(:,1) = meanpar2(:,1)./max(meanpar2(nCon,1));
            medpar1 = (medpar1+meanpar1)/2;
            medpar2 = (medpar2+meanpar2)/2;
            m1 = str2func(sizeFits(1,1).m1);
            m2 = str2func(sizeFits(1,1).m2);

            figure(7);if i==1;clf;end
            subplot(2,2,i)
            for iCon=1:nCon
                plot(szRng,m1(medpar1(iCon,:),szRng))
                hold on
            end
            ax = gca;
            ax.ColorOrderIndex = 1;
            for iCon=1:nCon
                plot(szRng,m2(medpar2(iCon,:),szRng))
                hold on
            end
            legend(num2str(cons'))
            title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
            xlabel('Size (deg)')
            ylabel('dF/F')
            xlim([0 max(szs)+1])
        end

        if sum(choosefig==8) % contrast-resp fits
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
            figure(8);if i==1;clf;end
            subplot(2,2,i)
            plot(C50f,C50r,'.')
            xlabel('C50f')
            ylabel('C50r')
            title({sprintf('Area:%s',areas(i));['(n=' num2str(length(cut)) ', n_{exp}=' num2str(nExpi) ')']})

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
        ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_6con);
        
        % cutoff by cellDist
        switch areas(i)
            case 'V1'
                cutoff = 10; %v1 cutoff at 10
            case {'LM','AL'}
                cutoff = 15; %alm cutoff at 15
            case 'PM'
                cutoff = 20; %pm cutoff at 20
        end
        ind = intersect(ind,find(cellDists_6con<cutoff));
        
        nExpi = length(expIndi);
        nCellsi = length(ind);
        sizeTune = sizeTune_6con{:,:,ind}; % (size,con,cell)
        sizeMean = sizeMean_6con(:,:,ind);
        sizeSEM = sizeSEM_6con(:,:,ind);
        sizeFits = sizeFits_6con(ind,:); %cell,con
        lbub_fits = lbub_fits_6con(ind,:,:); %cell,par,val (low up mean true stdev)
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
            conResp = conResp_6con(ind,:);
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
iCell = 12; %randi(sum(nCellsExp),1);
cellFit = sizeFits_all(iCell,nCon);
figure(8);clf;
plot(cellFit.szs0,cellFit.data,'k.')
hold on
plot(szRng,cellFit.fitout1,'r');
plot(szRng,cellFit.fitout2,'b');
errorbar(szs,sizeMean_all(:,nCon,iCell),sizeSEM_all(:,iCon,iCell),'k')
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
    BLr = baseline(cut)./(baseline(cut)+Rmax(cut))
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
c_areas = categorical({'V1' 'LM' 'AL' 'PM' 'all'},{'V1' 'LM' 'AL' 'PM' 'all'});
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
errorbar(1:4,C50mean,C50SEM,'x')
plot([1 4],[0.82 0.82],'-k', 'LineWidth',2)
plot(2.5,0.85,'*k')
hold off
set(gca,'XTick',1:4,'XTickLabel',{['V1 (n=' num2str(nCut(1)) ')'],['LM (n=' num2str(nCut(2)) ')'],['AL (n=' num2str(nCut(3)) ')'],['PM (n=' num2str(nCut(4)) ')']})
title('C_{50} by area')
ylabel('C_{50}')
legend('mean+SEM')
ylim([0 1])

figure(17);clf;
boxplot(yBLr,x)
hold on
errorbar(1:4,BLrmean,BLrSEM,'x')
plot([1 4],[0.93 0.93],'-k', 'LineWidth',2)
plot(2.5,0.95,'*k')
plot([2 4],[0.88 0.88],'-k', 'LineWidth',2)
plot(3,0.9,'*k')
plot([3 4],[0.83 0.83],'-k', 'LineWidth',2)
plot([3.4 3.5 3.6],[0.85 0.85 0.85],'*k')
hold off
set(gca,'XTick',1:4,'XTickLabel',{['V1 (n=' num2str(nCut(1)) ')'],['LM (n=' num2str(nCut(2)) ')'],['AL (n=' num2str(nCut(3)) ')'],['PM (n=' num2str(nCut(4)) ')']})
title('BL/(BL+R_{max}) by area')
ylabel('BL/(BL+R_{max})')
legend('mean+SEM')
ylim([-0.1 1])

%% stats on C50, BLr
[p,~,statC50] = anova1(yC50,x)
[results, means] = multcompare(statC50,'CType','hsd')

[p,~,statBLr] = anova1(yBLr,x)
[results, means] = multcompare(statBLr,'CType','hsd')
%%

for i=1:4
    for j=i+1:4
        comb = [i j];
        df = nCut(i)+nCut(j)-2;
        tC50 = (C50mean(i)-C50mean(j))./sqrt(C50SEM(i).^2+C50SEM(j).^2);
        tBLr = (BLrmean(i)-BLrmean(j))./sqrt(BLrSEM(i).^2+BLrSEM(j).^2);
        
        fprintf('\n%s v %s\n Deg. f = %d\ntC50 = %.4f\ntBLr = %.4f\ntCrit = %.4f (0.05), %.4f (0.01), %.4f (0.001)\n',areas(i),areas(j),df,tC50,tBLr,tinv(0.05,df),tinv(0.01,df),tinv(0.001,df))
    end
end