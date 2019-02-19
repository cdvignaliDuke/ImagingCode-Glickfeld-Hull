% new contrast response analysis
% only take curve at 20 degrees
% using all cells regardless of size-tuning curve fits

%% load csv to define exp

clear all; clc;
%expfile = '\\CRASH.dhe.duke.edu\data\home\kevin\Code\Ai9x_experiment_list.txt';
%expfile = '\\CRASH.dhe.duke.edu\data\home\kevin\Code\Ai9x_experiment_list_6conOnly_LMAL.txt';
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
    if size(sizeTune,2) == 6 % if 6 con, only take middle 4 cons (2-5) % CHANGE TO 7 INSTEAD OF 6 FOR 6 CON DATA ONLY
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

nCellsTot = sum(nCellsExp);
fprintf('%d cells loaded\n',nCellsTot)

expInd = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExp(i))];
end

cons = [0.1 0.2 0.4 0.8]; nCon = length(cons);
%cons = [0.05 0.1 0.2 0.4 0.8 1]; nCon = length(cons);
szs = 5*1.5.^[0:7]; nSize = length(szs);
szRng = linspace(0,max(szs));

%% con resp for all cells
override = 1;
conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));
conRng = 0.001:0.001:1;
opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
% if exist(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', 'conStruct_4con.mat'),'file') && ~override
%     % load contrast response instead of compute
%     fprintf('Found previous fits, loading...\n')
%     filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', 'conStruct_4con.mat');
%     load(filename);
% else
fprintf('Extracting contrast response of each cell at 25 deg\n')
s4 = zeros(1,4);
s = zeros(1);
conStruct_all = struct('resp25',s4,'fit',s4,'C50r',s,'Rsq',s,'x0',s4);
conStruct_all(nCellsTot) = conStruct_all;

shuf_ind = cell(nCon,1);
%nTr = zeros(1,nCon);
iSz = 5; %corresponds to 25.3 deg
tic
for iCell=1:nCellsTot
    iCell
    if sum(iCell==[1 nCellsTot]) % check first and last to reset zeros to blank
        conStruct_all(iCell).resp25 = [];conStruct_all(iCell).fit = [];conStruct_all(iCell).C50r = [];conStruct_all(iCell).Rsq = [];conStruct_all(iCell).x0 = [];conStruct_all(iCell).sig = [];
    end
    
    nTr = cellfun(@(c) length(c), sizeTune_all(iSz,:,iCell));
    sig_resp = zeros(1,nCon);
    
    % recreate con data with resampled indices
    nPts = floor(mean(nTr));
    dumdum = zeros(1,nPts);%[0]; % define zero point for dF/F
    cons0 = zeros(1,nPts);%[0.1]; % and size
    for iCon = 1:nCon
        nPts = nPts + nTr(iCon);
        dum = sizeTune_all{iSz,iCon,iCell}';
        dumdum = [dumdum dum];
        cons0 = [cons0 cons(iCon)*ones(1,nTr(iCon))];
        %sig_resp(iCon) = ttest(zeros(1,nTr(iCon)),dum,'Tail','left','Alpha',0.05/(nCon-1));
    end
    sig_resp(nCon) = ttest(zeros(1,nTr(nCon)),dum,'Tail','left','Alpha',0.05);
    if sum(sig_resp)
        conStruct_all(iCell).sig = 1;
    else
        conStruct_all(iCell).sig = 0;
    end
    
    conStruct_all(iCell).resp25 = cellfun(@(c) mean(c), sizeTune_all(iSz,:,iCell));
    cRi = conStruct_all(iCell).resp25;
    lb = [0 0 0 1]; %[0 0 0.1 1]
    ub = [Inf Inf 1 Inf]; %[Inf Inf 0.8 Inf]
    SStot = sum((cRi-mean(cRi)).^2);
    R2best = -Inf;
    for i=1%1:4
        x0 = [cRi(1) max(cRi) 0.1+0.1*i 3]; %BL Rmax C50 n
        %[cF, res] = lsqcurvefit(conModelH,x0,cons,cRi,lb,ub,opts);
        [cF, res] = lsqcurvefit(conModelH,x0,cons0,dumdum,lb,ub,opts);
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
    conStruct_all(iCell).R50 = R50;
    
    %store raw data
    conStruct_all(iCell).cons0 = cons0;
    conStruct_all(iCell).dumdum = dumdum;
    
    %plot con-resp with fit and C50
    %         plot(conRng,fitout,'r-')
    %         hold on
    %         %plot(cons,cRi,'ok')
    %         for iCon=1:nCon
    %             if sig_resp(iCon)
    %                 plot(cons(iCon),cRi(iCon),'or')
    %             else
    %                 plot(cons(iCon),cRi(iCon),'ok')
    %             end
    %         end
    %         plot(cons0,dumdum,'k.')
    %         line([C50 C50],[0 R50],'Color','g')
    %         plot(C50,R50,'gx')
    %         hold off
    %         title(num2str(iCell))
    %         xlabel('Con')
    %         ylabel('resp')
    %         xlim([0 1])
    %         pause
end
fprintf('Done, not saving...\n')
toc

% currently have it run
% next steps, do bootstrapping to determine those with consistent fits
% i.e. con within one octave?

% for now just look at average C50 by area
%% average c50 by area

fprintf('Examine cells from each area:\n')
areas = ["V1","LM","AL","PM"];

expInd = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExp(i))];
end

x=[];y_C50=[];

C50_all = [conStruct_all.C50r];

sig_ind = find([conStruct_all.sig]);
for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = find(ismember(expInd,expIndi));
    ind=intersect(ind,sig_ind);
    
    % no cutoff for now
    % cutoff by cellDist
    % try looking with different cutoffs
    switch areas(i)
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
        case {'LM','AL'}
            cutoff = 15; %alm cutoff at 15
        case 'PM'
            cutoff = 20; %pm cutoff at 20
    end
    %ind = intersect(ind,find(cellDists_all<cutoff));
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    cellDists = cellDists_all(ind);
    C50area = C50_all(ind)';
    
    cons_c = categorical({'0.1' '0.2' '0.4' '0.8'});
    
    C50_area{i} = C50area;
    x = [x; i*ones(size(C50area))];
    y_C50 = [y_C50; C50area];
    
    % C50 histogram
    figure(2);if i==1;clf;end
    subplot(2,2,i)
    histogram(C50area,[0:0.05:1])
    hold on
    title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    xlabel('C50')
    ylabel('num cells')
    ylim([0 nCellsi/6])
    
end

figure(1);clf;
boxplot(y_C50,x)
hold on
y_mean = [mean(C50_area{1}) mean(C50_area{2}) mean(C50_area{3}) mean(C50_area{4})];
y_std = [std(C50_area{1}) std(C50_area{2}) std(C50_area{3}) std(C50_area{4})];
plot(1:4,y_mean,'x')
hold off
set(gca,'XTick',1:4,'XTickLabel',{['V1 (n=' num2str(length(C50_area{1})) ')'],['LM (n=' num2str(length(C50_area{2})) ')'],['AL (n=' num2str(length(C50_area{3})) ')'],['PM (n=' num2str(length(C50_area{4})) ')']})
set(gcf,'Color','w')
%set(gca,'XAxisLocation','bottom','YAxisLocation','left','TickDir','out')
set(gca,'box','off','TickDir','out')
%title('RF radius')
ylabel('C50')
%%
[p,~,stat_C50] = anova1(y_C50,x)
[results, means] = multcompare(stat_C50,'CType','hsd')
areas
%% bootstrap
iSz = 5; %corresponds to 25.3 deg

Fit_struct = [];
stemp = struct('fit',s4,'C50r',s,'Rsq',s,'x0',s4);
Nshuf = 500;
% for iCell=[1:nCellsTot] % check first and last to reset zeros to blank
%     Fit_struct(iCell).Shuf(1).s = stemp;
% end

% store # trials at each size and highest con (same for all cells)


cd('K:\Code')
opts = optimset('Display','off');

fprintf('\nBegin shuffling...\n')
%figure;
for count_shuf = 1:Nshuf
    fprintf(['count_shuf: ' num2str(count_shuf) '/' num2str(Nshuf) '\n'])
    %stemp = struct('fit',s4,'C50r',s,'Rsq',s,'x0',s4);
    for iCell = 1:nCellsTot
        if sum(iCell==[1 nCellsTot])
            stemp(iCell).fit = zeros(1,4); stemp(iCell).C50r = []; stemp(iCell).Rsq = [];
        end
        if ~conStruct_all(iCell).sig
            continue
        end
        nTr = cellfun(@(c) length(c), sizeTune_all(iSz,:,iCell));
        shuf_ind = cell(nCon,1);
        for iCon = 1:nCon
            shuf_ind{iCon} = randsample(nTr(iCon),nTr(iCon),1);
            % fprintf([num2str(iCon) ': ' num2str(shuf_ind{iCon}') ' num2str(shuf ind)\n'])
        end
        
        % recreate con data with resampled indices
        nPts = floor(mean(nTr(:)));
        dumdum = zeros(1,nPts);%[0]; % define zero point for dF/F
        cons0 = zeros(1,nPts);%[0.1]; % and size
        for iCon = 1:nCon
            nPts = nPts + nTr(iCon);
            dum = sizeTune_all{iSz,iCon,iCell}';
            dumdum = [dumdum dum(shuf_ind{iCon})];
            cons0 = [cons0 cons(iCon)*ones(1,nTr(iCon))];
        end
        
        cRi = conStruct_all(iCell).resp25;
        lb = [0 0 0 1];
        ub = [Inf Inf 1 Inf];
        %SStot = sum((cRi-mean(cRi)).^2);
        %R2best = -Inf;
        %for i=1%1:4
        x0 = [cRi(1) max(cRi) 0.1+0.1*1 3]; %BL Rmax C50 n
        [cF, res] = lsqcurvefit(conModelH,x0,cons0,dumdum,lb,ub,opts);
        R2 = 1-res/SStot;
        %if R2>R2best
        %    R2best = R2;
        %    cFbest = cF;
        %    x0best = x0;
        %end
        %end
        %cF = cFbest;
        %R2 = R2best;
        %conStruct_shuf(iCell).fit = cF;
        %conStruct_shuf(iCell).Rsq = R2;
        %conStruct_shuf(iCell).x0 = x0best;
        stemp(iCell).fit = cF;
        stemp(iCell).Rsq = R2;
        stemp(iCell).x0 = x0best;
        
        fitout = conModelH(cF,conRng);
        R50 = fitout(1)+(fitout(end)-fitout(1))/2;
        i50 = find(abs(fitout - R50) == min(abs(fitout - R50)),1);
        C50 = conRng(i50);
        %conStruct_shuf(iCell).C50r = C50;
        stemp(iCell).C50r=C50;
    end
    Fit_struct.Shuf(count_shuf).s = stemp;
    
end
fprintf('\nShuffling done, now saving fit results\n')

Fit_struct.True.s = conStruct_all;
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', 'Fit_struct_4con.mat');
save(filename, 'Fit_struct')

%% load data
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', 'Fit_struct_4con.mat');
load(filename, 'Fit_struct')
%conStruct_all = Fit_struct.True.s;
nCellsTot = length(conStruct_all);
%cons = [0.05 0.1 0.2 0.4 0.8 1];
cons = [0.1 0.2 0.4 0.8];
nCon = length(cons);
fprintf('Loaded %d-con data, %d cells\n',nCon,nCellsTot);
Nshuf=length(Fit_struct.Shuf);

%% assess shuffles
fprintf('Assessing goodness of fit\n')
if Nshuf>1
    fprintf('Reading in variables of interest\n')
    fit_true_vec = zeros(nCellsTot, 6);
    for iCell = 1:nCellsTot
        %if ~conStruct_all(iCell).sig
        %    continue
        %end
        if ~isempty(Fit_struct.True.s(iCell))
            eval('tmp = Fit_struct.True.s(iCell).fit(:);');
            eval('tmp = [tmp; Fit_struct.True.s(iCell).C50r];');
            eval('tmp = [tmp; Fit_struct.True.s(iCell).Rsq];');
            
            % BL Rmax C50f n C50r Rsq
            fit_true_vec(iCell,:) = tmp';
        end
    end
    
    fit_shuf_vec = NaN(nCellsTot,6,Nshuf);
    for count_shuf = 1:Nshuf
        for iCell = 1:nCellsTot
            if ~isempty(Fit_struct.Shuf(Nshuf).s(iCell).C50r)
                eval('tmp = Fit_struct.Shuf(count_shuf).s(iCell).fit(:);');
                eval('tmp = [tmp; Fit_struct.Shuf(count_shuf).s(iCell).C50r];');
                eval('tmp = [tmp; Fit_struct.Shuf(count_shuf).s(iCell).Rsq];');
                
                % prefSize PS1 PS2 suppInd SI1 SI2 Fscore Ftest
                fit_shuf_vec(iCell,:,count_shuf) = tmp';
            end
        end
    end
    
    %% plot cells with size tuning curve and shuffle results
    chosen = find([conStruct_all.sig]);
    Npars = size(fit_shuf_vec,2);
    lbub_fits = NaN(nCellsTot,Npars,5);
    alpha_bound = .025;
    ind_shuf_lb = ceil(Nshuf*alpha_bound); % 0.025 percentile
    ind_shuf_ub = ceil(Nshuf*(1-alpha_bound)); % 0.975 percentile
    for iCell = 1:nCellsTot
        for count2 = 1:Npars
            tmp = squeeze(fit_shuf_vec(iCell,count2,:));
            [i,j] = sort(tmp); % sort in order
            lbub_fits(iCell,count2,1) = i(ind_shuf_lb); %lower 0.025
            lbub_fits(iCell,count2,2) = i(ind_shuf_ub); %upper 0.975
            lbub_fits(iCell,count2,3) = mean(i); %mean
            lbub_fits(iCell,count2,5) = std(i); %stdev
        end
        lbub_fits(iCell,:,4) = fit_true_vec(iCell,:); % true (no shuffle)
        
        if 0%sum(iCell==chosen)
            iCell
            conStruct_all(iCell).sig
            stemp = Fit_struct.True.s(iCell);
            figure(1);clf;
            subplot(2,4,1)
            % insert plot of con-resp and fit
            fitout = conModelH(stemp.fit,conRng);
            plot(conRng,fitout,'r-')
            hold on
            plot(cons,stemp.resp25,'ok')
            plot(stemp.cons0,stemp.dumdum,'k.')
            line([stemp.C50r stemp.C50r],[0 stemp.R50],'Color','g')
            plot(stemp.C50r,stemp.R50,'gx')
            hold off
            title(num2str(iCell))
            xlabel('Con')
            ylabel('resp')
            xlim([0 1])
            subplot(2,4,3)
            %insert histogram of BL values, mark confidence interval
            histogram(squeeze(fit_shuf_vec(iCell,1,:)),50);
            hold on
            line([lbub_fits(iCell,1,1) lbub_fits(iCell,1,1)], [0 Nshuf],'Color','red')
            line([lbub_fits(iCell,1,2) lbub_fits(iCell,1,2)], [0 Nshuf],'Color','red')
            line([lbub_fits(iCell,1,4) lbub_fits(iCell,1,4)], [0 Nshuf],'Color','green')
            hold off
            title('BL')
            subplot(2,4,4)
            % histogram of Rmax
            histogram(squeeze(fit_shuf_vec(iCell,2,:)),50);
            hold on
            line([lbub_fits(iCell,2,1) lbub_fits(iCell,2,1)], [0 Nshuf],'Color','red')
            line([lbub_fits(iCell,2,2) lbub_fits(iCell,2,2)], [0 Nshuf],'Color','red')
            line([lbub_fits(iCell,2,4) lbub_fits(iCell,2,4)], [0 Nshuf],'Color','green')
            hold off
            title('Rmax')
            subplot(2,4,5)
            % histogram of C50f
            histogram(squeeze(fit_shuf_vec(iCell,3,:)),50);
            hold on
            line([0.5*lbub_fits(iCell,3,4) 0.5*lbub_fits(iCell,3,4)], [0 Nshuf],'Color','red')
            line([2*lbub_fits(iCell,3,4) 2*lbub_fits(iCell,3,4)], [0 Nshuf],'Color','red')
            line([lbub_fits(iCell,3,1) lbub_fits(iCell,3,1)], [0 Nshuf],'Color','green')
            line([lbub_fits(iCell,3,2) lbub_fits(iCell,3,2)], [0 Nshuf],'Color','green')
            line([lbub_fits(iCell,3,4) lbub_fits(iCell,3,4)], [0 Nshuf],'Color','green')
            hold off
            title('C50f')
            xlim([0 1])
            subplot(2,4,6)
            % histogram of n
            histogram(squeeze(fit_shuf_vec(iCell,4,:)),50);
            hold on
            line([lbub_fits(iCell,4,1) lbub_fits(iCell,4,1)], [0 Nshuf],'Color','red')
            line([lbub_fits(iCell,4,2) lbub_fits(iCell,4,2)], [0 Nshuf],'Color','red')
            line([lbub_fits(iCell,4,4) lbub_fits(iCell,4,4)], [0 Nshuf],'Color','green')
            hold off
            title('n')
            subplot(2,4,7)
            % histogram of C50r
            histogram(squeeze(fit_shuf_vec(iCell,5,:)),50);
            hold on
            line([0.5*lbub_fits(iCell,5,4) 0.5*lbub_fits(iCell,5,4)], [0 Nshuf],'Color','red')
            line([2*lbub_fits(iCell,5,4) 2*lbub_fits(iCell,5,4)], [0 Nshuf],'Color','red')
            line([lbub_fits(iCell,5,1) lbub_fits(iCell,5,1)], [0 Nshuf],'Color','green')
            line([lbub_fits(iCell,5,2) lbub_fits(iCell,5,2)], [0 Nshuf],'Color','green')
            line([lbub_fits(iCell,5,4) lbub_fits(iCell,5,4)], [0 Nshuf],'Color','green')
            hold off
            title('C50r')
            xlim([0 1])
            subplot(2,4,8)
            % histogram of Rsq
            histogram(squeeze(fit_shuf_vec(iCell,2,:)),50);
            hold on
            line([lbub_fits(iCell,6,1) lbub_fits(iCell,6,1)], [0 Nshuf],'Color','red')
            line([lbub_fits(iCell,6,2) lbub_fits(iCell,6,2)], [0 Nshuf],'Color','red')
            line([lbub_fits(iCell,6,4) lbub_fits(iCell,6,4)], [0 Nshuf],'Color','green')
            hold off
            title('Rsq')
            xlim([0 1])
            
            pause
        end
    end
end

%% determine good fits
% first base on significance, then check C50
% check bounds of confidence interval are within 1 octave of "True" fit
% e.g. lower(2.5th)>0.5*C50true and upper(97.5th)<2*C50true
goodfit_ind_con = [];
for iCell = 1:nCellsTot
    if ~conStruct_all(iCell).sig
        continue
    end
    lOct = 0.5*lbub_fits(iCell,5,4); %5=C50r, 4=True fit, lower octave
    hOct = 2*lbub_fits(iCell,5,4); %5=C50r, 4=True fit, upper octave
    
    if (lbub_fits(iCell,5,1)>lOct) && (lbub_fits(iCell,5,2)<hOct) %5=C50r, 1=lb/2=ub
        if ~(conStruct_all(iCell).resp25(1)<conStruct_all(iCell).resp25(4))
            continue
        end
        goodfit_ind_con = [goodfit_ind_con iCell];
    end
end

fprintf(['#Good cells = ' num2str(length(goodfit_ind_con)) '\nNot Saving good fits\n'])

%% examine all cells

figure(4);clf;
histogram(log(lbub_fits(goodfit_ind_con,2,4)./lbub_fits(goodfit_ind_con,1,4)))
xlim([0 1])

fprintf('Examine cells from each area:\n')
areas = ["V1","LM","AL","PM"];

expInd = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExp(i))];
end

x=[];y_C50=[];

C50_all = [conStruct_all.C50r];

sig_ind = find([conStruct_all.sig]);
    notgood_ind = setdiff(sig_ind,goodfit_ind_con);
for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = find(ismember(expInd,expIndi));
    ind=intersect(ind,goodfit_ind_con);
    %ind = intersect(ind,notgood_ind);
    
    % no cutoff for now
    % cutoff by cellDist
    % try looking with different cutoffs
    switch areas(i)
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
        case {'LM','AL'}
            cutoff = 15; %alm cutoff at 15
        case 'PM'
            cutoff = 20; %pm cutoff at 20
    end
    %ind = intersect(ind,find(cellDists_all<cutoff));
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    cellDists = cellDists_all(ind);
    C50area = C50_all(ind)';
    
    cons_c = categorical({'0.1' '0.2' '0.4' '0.8'});
    
    C50_area{i} = C50area;
    x = [x; i*ones(size(C50area))];
    y_C50 = [y_C50; C50area];
    
    % C50 histogram
    figure(2);if i==1;clf;end
    subplot(2,2,i)
    histogram(C50area,[0:0.05:1])
    hold on
    title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    xlabel('C50')
    ylabel('num cells')
    ylim([0 nCellsi/6])
    
    resp = reshape([conStruct_all(ind).resp25],nCon, length(ind))';
    respnorm = resp./max(resp,[],2); % normalize by max
    %respnorm = resp./(resp(:,nCon)); % normalize by highest con
    respsem = std(respnorm)./sqrt(nCellsi);
    figure(3);if i==1;clf;end % plot of all good cell's raw contrast data normalized
    subplot(2,2,i)
    for iCell = 1:nCellsi
        if lbub_fits(ind(iCell),2,4)/lbub_fits(ind(iCell),1,4)<1
            plot1 = plot(cons,respnorm(iCell,:),'r--');
        else
            plot1 = plot(cons,respnorm(iCell,:),'b');
        end
        hold on
        plot1.Color(4) = 0.2;
    end
    hold on
    %plot(cons,mean(xnorm),'r') %overall mean in red
    errorbar(cons,mean(respnorm),respsem,'r')
    xlim([0 1])
    ylim([-0.1 1.1])
    title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    xlabel('Contrast')
    ylabel('Resp. norm')
    
    figure(4);if i==1;clf;end % plot of all good cell's raw contrast data normalized
    subplot(2,2,i)
    for iCell = 1:nCellsi
        fitout = conModelH(conStruct_all(ind(iCell)).fit,conRng);
        fitoutnorm = fitout./max(fitout);
        if lbub_fits(ind(iCell),2,4)/lbub_fits(ind(iCell),1,4)<1
            plot1 = plot(conRng,fitoutnorm,'r--');
        else
            plot1 = plot(conRng,fitoutnorm,'b');
        end
        hold on
        plot1.Color(4) = 0.2;
    end
    hold on
    %plot(cons,mean(xnorm),'r') %overall mean in red
    errorbar(cons,mean(respnorm),respsem,'r')
    xlim([0 1])
    ylim([-0.1 1.1])
    title({sprintf('Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    xlabel('Contrast')
    ylabel('Fit norm')
    
    %     figure(4);if i==1;clf;end
    %     subplot(4,2,2*i-1)
    %     histogram(lbub_fits(ind,1,4))
    %     xlim([0 1])
    %     title('BL')
    %     subplot(4,2,2*i)
    %     histogram(lbub_fits(ind,2,4))
    %     title('Rmax')
    %     xlim([0 1])
    %
    %     figure(5);if i==1;clf;end
    %     subplot(2,2,i)
    %     histogram(lbub_fits(ind,2,4)./lbub_fits(ind,1,4),[0.1:0.1:2])
    %     title('Rmax/BL')
    %     xlim([0 2])
    
end

figure(1);clf;
boxplot(y_C50,x)
hold on
y_mean = [mean(C50_area{1}) mean(C50_area{2}) mean(C50_area{3}) mean(C50_area{4})];
y_std = [std(C50_area{1}) std(C50_area{2}) std(C50_area{3}) std(C50_area{4})];
plot(1:4,y_mean,'x')
hold off
set(gca,'XTick',1:4,'XTickLabel',{['V1 (n=' num2str(length(C50_area{1})) ')'],['LM (n=' num2str(length(C50_area{2})) ')'],['AL (n=' num2str(length(C50_area{3})) ')'],['PM (n=' num2str(length(C50_area{4})) ')']})
set(gcf,'Color','w')
%set(gca,'XAxisLocation','bottom','YAxisLocation','left','TickDir','out')
set(gca,'box','off','TickDir','out')
%title('RF radius')
ylabel('C50')
%%
[p,~,stat_C50] = anova1(y_C50,x)
[results, means] = multcompare(stat_C50,'CType','hsd')
areas

%%
figure(6);clf;
CIdiff = lbub_fits(sig_ind,5,2) - lbub_fits(sig_ind,5,1);
respMax = max(reshape([conStruct_all(sig_ind).resp25],nCon, length(sig_ind))',[],2);
C50x = lbub_fits(sig_ind,5,4);
for iCell = 1:length(sig_ind)
    if intersect(sig_ind(iCell),goodfit_ind_con)
        plot(CIdiff(iCell), C50x(iCell), '.r')
        hold on
    else
        plot(CIdiff(iCell), C50x(iCell), '.b')
        hold on
    end
end
xlabel('C50 CIdiff')
ylabel('C50 true')
ylim([0 1])

figure(7);clf;
respMax = max(reshape([conStruct_all(sig_ind).resp25],nCon, length(sig_ind))',[],2);
ntrue = lbub_fits(sig_ind,4,4);
for iCell = 1:length(sig_ind)
    if intersect(sig_ind(iCell),goodfit_ind_con)
        plot(respMax(iCell), ntrue(iCell), '.r')
        hold on
    else
        plot(respMax(iCell), ntrue(iCell), '.b')
        hold on
    end
end
xlabel('resp Max')
ylabel('n true') 

%next step, implement smoothness penalty or n restriction
% think start with n restriction, just examine some made up graphs (based
% on fits) but change n to be higher than 100 or low etc, just find a good
% range
% third option is to just fit the mean data points, idk, need to look at
% retinotopy fitting methods and size tuning