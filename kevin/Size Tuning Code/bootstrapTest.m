% bootstrap method for fitting size-tuning curves
% flow:
% 1. load size tuning data (sizeTuneData.mat -> sizeTune, sizeMean, sizeSEM, cellDists
% 2. then fit first run through, plot cells with mean curves in different
% colors and fits in corresponding colors? or do four subplots for each con
% 3. then begin shuffling... storing prefSize?
%

fprintf('\nBegin fitting size-tuning data...\n')

Fit_struct = [];

Nshuf = 500;
fprintf(['Nshuf = ' num2str(Nshuf) '\n'])

% store # trials at each size and highest con (same for all cells)
nTr = zeros(nSize,1);
for iSz = 1:nSize
    nTr(iSz) = length(sizeTune{iSz,nCon,1});
end
fprintf(['Sizes: ' num2str(szs) '\n# trials: ' num2str(nTr')])
shuf_ind = cell(nSize,1);


cd('K:\Code')
opts = optimset('Display','off');

fprintf('\nBegin shuffling...\n')
for count_shuf = 0:Nshuf
    fprintf(['count_shuf: ' num2str(count_shuf) '/' num2str(Nshuf) '\n'])
    for iSz = 1:nSize
        if count_shuf > 0
            shuf_ind{iSz} = randsample(nTr(iSz),nTr(iSz),1); % resample with replacement
        else
            shuf_ind{iSz} = 1:nTr(iSz); % shuf_count==0 use all trials
        end
    end
    
    ifig = 1;
    start = 1;
    for iCell = 1:nCells
        %fprintf(num2str(iCell));
        % recreate sizeTune data with resampled indices
        nPts = floor(mean(nTr));
        dumdum = zeros(1,nPts);%[0]; % define zero point for dF/F
        szs0 = zeros(1,nPts);%[0.1]; % and size
        for iSz = 1:nSize
            nPts = nPts + nTr(iSz);
            dum = sizeTune{iSz,nCon,iCell}';
            dumdum = [dumdum dum(shuf_ind{iSz})];
            szs0 = [szs0 szs(iSz)*ones(1,nTr(iSz))];
        end
        
        % max of each size mean for use in initial guesses
        maxMean = max(sizeMean(:,nCon,iCell));
        
        if count_shuf == 0
            PLOTIT_FIT = 1;
            SAVEALLDATA = 1;
            Fit_SizeTuneSmooth_KM % call fit script, returns fit structure s, saves plots to kevin analysis folder
            eval(['Fit_struct(iCell).True.s_',' = s;']);
        else
            SAVEALLDATA = 0;
            PLOTIT_FIT = 0;
            Fit_SizeTuneSmooth_KM
            eval(['Fit_struct(iCell).Shuf(count_shuf).s_',' = s;']);
        end
    end
    if count_shuf == 0
        set(gcf, 'Position', [0 0 800 1000]);
        fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SizeTuneFits' num2str(ifig) '.pdf']);
        print(fn_out,'-dpdf')
    end
end
fprintf('\nShuffling done, saving fit results\n')

fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Fit_struct.mat']);
save(fn_out, 'Fit_struct')

%% assess
% extract values for prefSize (use for decision), Ftest
% also prefSize(1/2), suppInd(1,2), fit1.c1/OF1, fit2.c2/OF2, Rsq12

fprintf('Assessing goodness of fit\n')
if Nshuf>1
    fprintf('Reading in variables of interest\n')
    for iCell = 1:nCells
        if ~isempty(Fit_struct(iCell).True)
            eval('tmp = Fit_struct(iCell).True.s_.prefSize;');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.prefSize1];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.prefSize2];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.suppInd];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.suppInd1];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.suppInd2];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.Fscore];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.Ftest];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.maxResp1];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.maxResp2];');
            
            % prefSize PS1 PS2 suppInd SI1 SI2 Fscore Ftest
            fit_true_vec(iCell,:) = tmp;
        end
    end
    
    fit_shuf_vec = NaN(nCells,10,Nshuf);
    for count_shuf = 1:Nshuf
        for iCell = 1:nCells
            if ~isempty(Fit_struct(iCell).Shuf)
                eval('tmp = Fit_struct(iCell).Shuf(count_shuf).s_.prefSize;');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.prefSize1];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.prefSize2];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.suppInd];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.suppInd1];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.suppInd2];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Fscore];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Ftest];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.maxResp1];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.maxResp2];');
                
                % prefSize PS1 PS2 suppInd SI1 SI2 Fscore Ftest
                fit_shuf_vec(iCell,:,count_shuf) = tmp;
            end
        end
    end
    %%
    chosen=[44 54]; %[89 71 75 67 72 64 83 31 41 73 77 45 46 52 79];
    Npars = size(fit_shuf_vec,2);
    lbub_fits = NaN(nCells,Npars,5);
    alpha_bound = .025;
    ind_shuf_lb = ceil(Nshuf*alpha_bound); % 0.025 percentile
    ind_shuf_ub = ceil(Nshuf*(1-alpha_bound)); % 0.975 percentile
    for iCell = 1:nCells
        if sum(iCell==chosen)
            s = Fit_struct(iCell).True.s_;
            figure(1);clf;
            subplot(3,3,1)
            errorbar([0 szs],[0 sizeMean(:,nCon,iCell)'],[0 sizeSEM(:,nCon,iCell)'])
            hold on
            plot(s.szs0,s.data,'.b')
            plot(szRng,s.fitout1,'-r')
            plot(szRng,s.fitout2,'-g')
            hold off
            ylim([min([-0.5*s.maxResp1 min(s.data)]) 1.2*max([s.maxResp2 max(s.data)])])
            title(['Cell #' num2str(iCell) ' Size Tuning @Con 0.8 (Ftest=' num2str(fit_true_vec(iCell,8)) ')']);
            xlabel('Stimulus size (deg)')
            ylabel('dF/F')
        end
        
        for count2 = 1:Npars
            tmp = squeeze(fit_shuf_vec(iCell,count2,:));
            [i,j] = sort(tmp); % sort in order
            lbub_fits(iCell,count2,1) = i(ind_shuf_lb); %lower 0.025
            lbub_fits(iCell,count2,2) = i(ind_shuf_ub); %upper 0.975
            lbub_fits(iCell,count2,3) = mean(i); %mean
            lbub_fits(iCell,count2,5) = std(i); %stdev
            
            if sum(iCell==chosen)
                switch count2
                    case 1
                        subplot(3,3,4)
                        histogram(i,0:max(szs)+1)
                        xlim([0 100]);
                        line([0.5*mean(i) 0.5*mean(i)], ylim, 'color','red')
                        line([2*mean(i) 2*mean(i)], ylim, 'color','red')
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['PrefSize shuffles cell #' num2str(iCell)])
                        xlabel('Size (deg)')
                        ylabel('count')
                    case 2
                        subplot(3,3,5)
                        histogram(i,0:max(szs)+1)
                        xlim([0 100]);
                        line([0.5*mean(i) 0.5*mean(i)], ylim, 'color','red')
                        line([2*mean(i) 2*mean(i)], ylim, 'color','red')
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['PrefSize1 shuffles cell #' num2str(iCell)])
                        xlabel('Size (deg)')
                    case 3
                        subplot(3,3,6)
                        histogram(i,0:max(szs)+1)
                        xlim([0 100]);
                        line([0.5*mean(i) 0.5*mean(i)], ylim, 'color','red')
                        line([2*mean(i) 2*mean(i)], ylim, 'color','red')
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['PrefSize2 shuffles cell #' num2str(iCell)])
                        xlabel('Size (deg)')
                    case 4
                        subplot(3,3,7)
                        histogram(i,0:0.01:2)
                        xlim([0 1])
                        line([mean(i) mean(i)], ylim)
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['SuppInd shuffles cell #' num2str(iCell)])
                        xlabel('SI')
                        ylabel('count')
                    case 5
                        subplot(3,3,8)
                        histogram(i,0:0.01:2)
                        xlim([0 1])
                        line([mean(i) mean(i)], ylim)
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['SuppInd1 shuffles cell #' num2str(iCell)])
                        xlabel('SI1')
                    case 6
                        subplot(3,3,9)
                        histogram(i,0:0.01:2)
                        xlim([0 1])
                        line([mean(i) mean(i)], ylim)
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['SuppInd2 shuffles cell #' num2str(iCell)])
                        xlabel('SI2')
                    case 7
                        Fcrit = finv(0.95,nPts-6,nPts-3);
                        subplot(3,3,2)
                        histogram(i,50)
                        xlim([0 max(i)]);
                        line([mean(i) mean(i)], ylim)
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        line([Fcrit Fcrit], ylim, 'LineWidth', 3, 'color', 'red')
                        title(['Fscore shuffles cell #' num2str(iCell)])
                        xlabel('Fscore')
                        ylabel('count')
                    case 8
                        subplot(3,3,3)
                        histogram(i,[-0.5 0.5 1.5])
                        xlim([-0.5 1.5])
                        ylim([0 Nshuf])
                        title(['Ftest shuffles cell #' num2str(iCell)])
                        xlabel('Ftest')
                        pause
                end
            end
        end
        lbub_fits(iCell,:,4) = fit_true_vec(iCell,:); % true (no shuffle)
    end
end

%% determine good fits
% first sort cells based on Ftest same as True fit in >50% of shuffles
% then use that model's prefSize confidence interval (e.g. prefSize1 or 2)
% check bounds of confidence interval are within 1 octave of "True" fit
% e.g. lower(2.5th)>0.5*prefSizeTrue and upper(97.5th)<2*prefSizeTrue
goodfit_ind_size = [];
for iCell = 1:nCells
    Ftest = lbub_fits(iCell,8,4); %8=Ftest, 4=True fit
    Ftestshuf = lbub_fits(iCell,8,3); %8=Ftest, 3=mean of shuffles
    lOct = 0.5*lbub_fits(iCell,1,4); %1=prefSize, 4=True fit, lower octave
    hOct = 2*lbub_fits(iCell,1,4); %1=prefSize, 4=True fit, upper octave
    switch Ftest
        case 0 % model1
            if Ftestshuf<0.5 %model1 >50% of shuffles
                if (lbub_fits(iCell,2,1)>lOct) && (lbub_fits(iCell,2,2)<hOct) %2=prefSize1, 1=lb/2=ub
                    goodfit_ind_size = [goodfit_ind_size iCell];
                end
            end
        case 1 % model2
            if Ftestshuf>0.5 % model2 >50% of shuffles
                if (lbub_fits(iCell,3,1)>lOct) && (lbub_fits(iCell,3,2)<hOct) %3=prefSize2, 1=lb/2=ub
                    goodfit_ind_size = [goodfit_ind_size iCell];
                end
            end
    end
end
fprintf(['#Good cells = ' num2str(length(goodfit_ind_size)) '\nSaving good fits\n'])

fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
save(fn_out, 'lbub_fits', 'goodfit_ind_size')

%% load

mouse = 'i842';
date = '180330';
ImgFolder = char('003');
time = char('1647');
nrun = size(ImgFolder,1);
run_str = catRunName(ImgFolder, nrun);
RetImgFolder = char('002');
nret = size(RetImgFolder,1);
ret_str = catRunName(RetImgFolder, nret);

% ret goodfits
fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_lbub_fits.mat']);
load(fn_out,'goodfit_ind')

% sizetune data
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
load(filename)

% size tune bootstrap Fit_struct
fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Fit_struct.mat']);
load(fn_out)

% bootstrap goodfits
fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
load(fn_out)

szs=5*1.5.^(0:7);

%% plots, change for size tuning
% ideas:
% prefSize + SI histograms (separate model1/2?)
% prefSize vs RF-distance
% SI vs RF-distance
% maxResp vs RF-distance
% maxResp vs prefSize? 
% ?

% is model1 + is model2
ism1 = find(~lbub_fits(goodfit_ind_size,8,4));
ism2 = find(lbub_fits(goodfit_ind_size,8,4));

cellDists_m1 = cellDists(goodfit_ind(goodfit_ind_size(ism1)));
cellDists_m2 = cellDists(goodfit_ind(goodfit_ind_size(ism2)));
prefSize_m1 = lbub_fits(goodfit_ind_size(ism1),1,4);
prefSize_m2 = lbub_fits(goodfit_ind_size(ism2),1,4);
suppInd_m1 = lbub_fits(goodfit_ind_size(ism1),4,4);
suppInd_m2 = lbub_fits(goodfit_ind_size(ism2),4,4);

% examine Ftest proportions across all cells
Ftest_allcells = lbub_fits(:,8,3); % 8=Ftest, 3=mean
figure(2);clf;
subplot(2,1,1)
histogram(Ftest_allcells,0:0.02:1)
title('Ftest proportion for all cells')
xlabel('Fraction of shuffles passing Ftest')
ylabel('Count')
subplot(2,1,2)
h1 = histogram(lbub_fits(goodfit_ind_size(ism1),8,3),0:0.02:1);
hold on
h2 = histogram(lbub_fits(goodfit_ind_size(ism2),8,3),0:0.02:1);
hold off
title('Ftest proportion for goodfit cells')
xlabel('Fraction of shuffles passing Ftest')
ylabel('Count')
legend('Model1','Model2','location','best')

% prefSize + SI
figure(3);clf;
subplot(2,1,1)
histogram(prefSize_m1,0:5:90)
hold on
histogram(prefSize_m2,0:5:90)
hold off
title('Preferred Size histogram')
xlabel('Pref Size (deg)')
ylabel('Count')
legend('Model1','Model2','location','best')
subplot(2,1,2)
histogram(suppInd_m1,0:0.05:2)
hold on
histogram(suppInd_m2,0:0.05:2)
hold off
title('Suppression Index histogram')
xlabel('SI')
ylabel('Count')
legend('Model1','Model2','location','best')

% prefSize vs RF-distance
figure(4);clf;
plot(cellDists_m1,prefSize_m1,'o')
hold on
plot(cellDists_m2,prefSize_m2,'o')
plot([0 max(szs)],[0 max(szs)],'--')
coefs1 = polyfit(cellDists_m1,prefSize_m1,1);
coefs2 = polyfit(cellDists_m2,prefSize_m2,1);
x1 = min(cellDists_m1):max(cellDists_m1);
x2 = min(cellDists_m2):max(cellDists_m2);
plot(x1,polyval(coefs1,x1),'b')
plot(x2,polyval(coefs2,x2),'r')
r1 = corrcoef(cellDists_m1,prefSize_m1);
r2 = corrcoef(cellDists_m2,prefSize_m2);
text(max(x1)+10,polyval(coefs1,max(x1)+10)+20,['y_1=' num2str(coefs1(1)) 'x+' num2str(coefs1(2))])
text(max(x2)+10,polyval(coefs2,max(x2)+10),['y_2=' num2str(coefs2(1)) 'x+' num2str(coefs2(2))])
text(max(x1)+10,polyval(coefs1,max(x1)+10)+15,['R^2_1=' num2str(r1(2))])
text(max(x2)+10,polyval(coefs2,max(x2)+10)-5,['R^2_2=' num2str(r2(2))])
hold off
title('Preferred Size vs RF-Stim distance')
xlabel('RF-stim dist (deg)')
ylabel('PrefSize (deg)')
legend('Model1','Model2','location','best')


% SI vs RF-distance
figure(5);clf;
plot(cellDists_m1,suppInd_m1,'o')
hold on
plot(cellDists_m2,suppInd_m2,'o')
coefs2 = polyfit(cellDists_m2,suppInd_m2,1);
x2 = min(cellDists_m2):max(cellDists_m2);
plot(x2,polyval(coefs2,x2),'r')
r2 = corrcoef(cellDists_m2,suppInd_m2);
text(max(x2)-5,polyval(coefs2,max(x2)-5)+0.2,['y_2=' num2str(coefs2(1)) 'x+' num2str(coefs2(2))])
text(max(x2)-5,polyval(coefs2,max(x2)-5)+0.1,['R^2_2=' num2str(r2(2))])
hold off
title('Suppression Index vs RF-Stim distance')
xlabel('RF-stim dist (deg)')
ylabel('SI')
legend('Model1','Model2','location','best')

% SI vs prefSize
figure(6);clf;
plot(prefSize_m1,suppInd_m1,'o')
hold on
plot(prefSize_m2,suppInd_m2,'o')
coefs2 = polyfit(prefSize_m2,suppInd_m2,1);
x2 = min(prefSize_m2):max(prefSize_m2);
plot(x2,polyval(coefs2,x2),'r')
r2 = corrcoef(prefSize_m2,suppInd_m2);
text(max(x2)-10,polyval(coefs2,max(x2)-10)+0.2,['y_2=' num2str(coefs2(1)) 'x+' num2str(coefs2(2))])
text(max(x2)-10,polyval(coefs2,max(x2)-10)+0.1,['R^2_2=' num2str(r2(2))])
hold off
title('Suppression Index vs Preferred Size')
xlabel('PrefSize (deg)')
ylabel('SI')
legend('Model1','Model2','location','best')


% to print
% set(gcf, 'Position', [0 0 800 1000]);
% fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_RFs.pdf']);
% print(fn_out,'-dpdf')