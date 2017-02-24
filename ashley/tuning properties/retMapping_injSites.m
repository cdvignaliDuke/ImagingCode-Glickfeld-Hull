%% load data and extract timecourses from dir tuning mask

awFSAVdatasets % reference for datasets
iexp = 10; % choose dataset from reference

try % create an analysis folder
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
catch
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging');
    cd(filedir)
    mkdir(date,ImgFolder)
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
end
fnout = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);

% register data to dir tuning FOV and use mask to get cell TCs
down = 5;
data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down

%register data
load(fullfile('Z:\analysis\',mouse,'two-photon imaging',date,expt(iexp).dirtuning,'regImg.mat'))
[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

%get time-courses
load(fullfile('Z:\analysis\',mouse,'two-photon imaging',date,expt(iexp).dirtuning,'mask&TCDir.mat'))
dataTC = stackGetTimeCourses(data_reg, mask_cell);

% neurpil subtraction
nCells = size(dataTC,2);
buf = 4;
np = 6;
neuropil = imCellNeuropil(mask_cell,buf,np);

npTC = zeros(size(dataTC));
for i = 1:size(dataTC,2)
    tempNPmask = squeeze(neuropil(:,:,i));
    if sum(sum(tempNPmask)) > 0
    npTC(:,i) = stackGetTimeCourses(data_reg,tempNPmask);
    end
end

%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(dataTC-tcRemoveDC(npTC*ii(i)));
end
[max_skew ind] =  max(x,[],1);
% skew(buf,:) = max_skew;
np_w = 0.01*ind;
npSubTC = dataTC-bsxfun(@times,tcRemoveDC(npTC),np_w);

%% trial-by-trial dF/F
nOn = double(input.nScansOn)./down;
nOff = double(input.nScansOff)./down;
ntrials = size(input.tGratingDirectionDeg,2);

data_reg = data_reg(:,:,1:ntrials*(nOn+nOff));
npSubTC = npSubTC(1:ntrials*(nOn+nOff),:);

sz = size(data_reg);
data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
data_df = bsxfun(@minus, double(data_tr), data_f); 
data_dfof = bsxfun(@rdivide,data_df, data_f); 
clear data_f data_df data_tr

szTC = size(npSubTC);
dataTC_tr = reshape(npSubTC,[nOn+nOff, ntrials, szTC(2)]);
dataTC_f = mean(dataTC_tr(nOff/2:nOff,:,:),1);
dataTC_df = bsxfun(@minus, double(dataTC_tr), dataTC_f); 
dataTC_dfof = bsxfun(@rdivide,dataTC_df, dataTC_f); 
clear dataTC_f dataTC_df dataTC_tr
%% variables
nCells = size(dataTC_dfof,3);

Az = celleqel2mat_padded(input.tGratingAzimuthDeg);
El = celleqel2mat_padded(input.tGratingElevationDeg);
Azs = unique(Az);
Els = unique(El);
if min(Els,[],2)<0
    Els = fliplr(Els);
end
nStim = length(Azs).*length(Els);

%% ret tuning
Stims = [];
data_dfof_avg = zeros(sz(1), sz(2), nOn+nOff, length(Azs).*length(Els));
start = 1;
for iEl = 1:length(Els)
    ind1 = find(El == Els(iEl));
    for iAz = 1:length(Azs)
        Stims = [Stims; Els(iEl) Azs(iAz)];
        ind2 = find(Az == Azs(iAz));
        ind = intersect(ind1,ind2);
        data_dfof_avg(:,:,:,start) = mean(data_dfof(:,:,:,ind),4);
        start = start +1;
    end
end
clear data_dfof
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_avg_ret = squeeze(mean(imfilter(data_dfof_avg(:,:,nOff:nOff+nOn,:),myfilter),3));

% time-courses for each cell, for each stim
tuning_mat = zeros(nStim, 2, nCells);
Ind_struct = [];
if nCells<37
    [n, n2] = subplotn(nCells);
else
    [n, n2] = subplotn(36);
end
tt= (1-nOff:nOn)*(1000./(expt(iexp).frame_rate/down));
figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
%         print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_TCs' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    for iStim = 1:nStim
        ind1 = find(Az == Stims(iStim,2));
        ind2 = find(El == Stims(iStim,1));
        ind = intersect(ind1,ind2);
        plot(tt', squeeze(mean(dataTC_dfof(:,ind,iCell),2)))
        hold on
        tuning_mat(iStim,1,iCell) = mean(mean(dataTC_dfof(nOff+1:nOn+nOff,ind,iCell),2),1);
        tuning_mat(iStim,2,iCell) = std(mean(dataTC_dfof(nOff+1:nOn+nOff,ind,iCell),1),[],2)./sqrt(length(ind));
        Ind_struct(iStim).all_trials = ind;
    end
    ylim([-0.05 0.5])
    ylabel('dF/F')
    xlabel('time from stim on (ms)')
    title(['cell ' num2str(iCell)])
    vline(nOff)
    start = start + 1;
end
print(fullfile(fnout, [date '_' mouse '_TCs' num2str(f) '.pdf']), '-dpdf')

% heatmap of responses
figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
        print(fullfile(fnout,[date '_' mouse '_Tuning' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    ret_mat = reshape(tuning_mat(:,1,iCell), [length(Azs) length(Els)]);
    ret_mat = ret_mat';
    h = imagesc(ret_mat);
    h.Parent.XTickLabel = strread(num2str(Azs),'%s');
    h.Parent.YTickLabel = strread(num2str(Els),'%s');
    colormap gray
    clim([0 max(max(tuning_mat(:,1,:),[],1),[],3)])
    colorbar
    title(num2str(chop(max(tuning_mat(:,1,iCell),[],1),2)))  
    start = start +1;
end
print(fullfile(fnout,[date '_' mouse '_Tuning' num2str(f) '.pdf']), '-dpdf')

%
Fit_struct = [];
[AzAz, ElEl] = meshgrid(Azs,Els); 
grid2.AzAz = AzAz;
grid2.ElEl = ElEl;

dAz = median(diff(Azs));
dEl = median(diff(Els));
Az_vec00 = Azs(1):(dAz/10):Azs(end);
El_vec00 = Els(1):(dEl/10):Els(end);
[AzAz00,ElEl00]=meshgrid(Az_vec00,El_vec00);
grid2.AzAz00 = AzAz00;
grid2.ElEl00 = ElEl00;
Nshuf = 500;
resp_dFoverF = squeeze(mean(dataTC_dfof(nOff:nOn+nOff,:,:),1))';
base_dFoverF = squeeze(mean(dataTC_dfof(nOff/2:nOff,:,:),1))';
p_ttest = zeros(nCells,nStim);
h_ttest = zeros(nCells,nStim);
h_all = zeros(1,nCells);
for count_shuf = 0:Nshuf
    fprintf('.')
    Im_mat_USE = zeros(nCells, nStim);
    for iCond = 1:nStim        
        ind_all = Ind_struct(iCond).all_trials;
        if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
            ind_all_1 = ind_all(randsample(length(ind_all),length(ind_all),1));
        else
            ind_all_1 = ind_all;
            [h_ttest(:,iCond) p_ttest(:,iCond)] = ttest(resp_dFoverF(:,ind_all), base_dFoverF(:,ind_all), 'tail', 'right', 'dim', 2, 'alpha', 0.05./(nStim-1));
        end
        Im_mat_USE(:,iCond) = mean(resp_dFoverF(:,ind_all_1),2);
    end

    start = 1;
    for iCell = 1:nCells;
        if count_shuf == 0
            if sum(h_ttest(iCell,:),2) == 0 
                ind_p = find(p_ttest(iCell,:)< 0.05./((nStim-1)/2));
                if length(ind_p)<2
                    ind_p = find(p_ttest(iCell,:)< 0.05./((nStim-1)/3));
                    if length(ind_p)<3
                        ind_p = find(p_ttest(iCell,:)< 0.05./((nStim-1)/4));
                        if length(ind_p)<4
                            h_all(1,iCell) = 0;
                        else
                            h_all(1,iCell) = 1;
                        end
                    else
                        h_all(1,iCell) = 1;
                    end
                else
                    h_all(1,iCell) = 1;
                end
            else
                h_all(1,iCell) = 1;
            end
        end
        if count_shuf>0
            if h_all(1,iCell) == 0
                continue
            end
        end
        a = Im_mat_USE(iCell,:);
        if max(a,[],2) > 0
            b = reshape(a',length(Azs),length(Els));
            data = b';
            if count_shuf == 0
                PLOTIT_FIT = 1;
                SAVEALLDATA = 1;
                Fit_2Dellipse_LG_Ret
                eval(['Fit_struct(iCell).True.s_',' = s;']);
            else
                SAVEALLDATA = 0;
                PLOTIT_FIT = 0;
                Fit_2Dellipse_LG_Ret
                eval(['Fit_struct(iCell).Shuf(count_shuf).s_',' = s;']);
            end
        end               
    end
    if count_shuf == 0
        fn_out = fullfile(fnout, ['RFfits' num2str(ifig) '.pdf']);   
        print(fn_out,'-dpdf')
    end
end

fn_out = fullfile(fnout, ['Fit_struct.mat']);   
save(fn_out, 'Fit_struct')

resp_ind = find(h_all);

 if Nshuf>1;
    for iCell = 1:nCells
        if ~isempty(Fit_struct(iCell).True)                
            eval(['tmp = Fit_struct(iCell).True.s_.x;']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.Elhicut_50];']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.Azhicut_50];']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.Elhicut_10];']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.Azhicut_10];']);
            fit_true_vec(iCell,:) = tmp;
        end
    end
    
    fit_shuf_vec = NaN(nCells,10,Nshuf);
    for count_shuf = 1:Nshuf
        for iCell = 1:nCells
            if ~isempty(Fit_struct(iCell).Shuf)
                eval(['tmp = Fit_struct(iCell).Shuf(count_shuf).s_.x;']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_50];']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_50];']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_10];']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_10];']);
                %fit is: %A sigma_Az sigma_El Az0 El0 xi
                fit_shuf_vec(iCell,:,count_shuf) = tmp;
            end
        end
    end

    Npars = size(fit_shuf_vec,2);
    lbub_fits = NaN(nCells,Npars,5);
    alpha_bound = .025;
    for iCell = 1:nCells
        for count2 = 1:Npars
            tmp = squeeze(fit_shuf_vec(iCell,count2,:));
            [i,j] = sort(tmp);
            ind_shuf_lb = ceil(Nshuf*alpha_bound);
            ind_shuf_ub = ceil(Nshuf*(1-alpha_bound));
            lbub_fits(iCell,count2,1) = i(ind_shuf_lb);
            lbub_fits(iCell,count2,2) = i(ind_shuf_ub);
            lbub_fits(iCell,count2,3) = mean(i); 
            lbub_fits(iCell,count2,5) = std(i);
        end
        %now take means from truedata fit:
        lbub_fits(iCell,:,4) = fit_true_vec(iCell,:);
    end
end

lbub_diff = lbub_fits(:,:,2)-lbub_fits(:,:,1);

goodfit_ind = [];
for iCell = 1:nCells
    if lbub_diff(iCell,4)<input.gratingAzimuthStepDeg*2
        if lbub_diff(iCell,5)<input.gratingAzimuthStepDeg*2
            goodfit_ind = [goodfit_ind iCell];
        end
    end
end

fn_out = fullfile(fnout, ['lbub_fits.mat']);   
save(fn_out, 'lbub_fits', 'lbub_diff', 'goodfit_ind', 'resp_ind')


figure
subplot(2,2,1)
for i = goodfit_ind
    plot(lbub_fits(i,4,4), lbub_fits(i,5,4), 'o')
    hold on;
    xlim([min(Azs,[],2) max(Azs,[],2)])
    ylim([min(Els,[],2) max(Els,[],2)])
end
axis equal
title('RF center')
subplot(2,2,2)
for i = goodfit_ind
    ellipse(lbub_fits(i,2,4), lbub_fits(i,3,4), 0, lbub_fits(i,4,4), lbub_fits(i,5,4));
    hold on;
end
axis equal
title('RF- 1 sigma')
fn_out = fullfile(fnout, ['RFs.pdf']);   
print(fn_out,'-dpdf')

figure;
subplot(2,2,1)
hist(lbub_fits(goodfit_ind,2,4))
title('Sigma azimuth')
subplot(2,2,2)
hist(lbub_fits(goodfit_ind,3,4))
title('Sigma elevation')
subplot(2,2,3)
a = lbub_fits(goodfit_ind,3,4).*lbub_fits(goodfit_ind,2,4).*pi;
hist(a)
title('Area')
subplot(2,2,4)
scatter(a, lbub_fits(goodfit_ind,1,4),'o')
xlabel('Area')
ylabel('Peak dF/F')
fn_out = fullfile(fnout, ['RFdists.pdf']);   
print(fn_out,'-dpdf')
