%% Parameters
orig_rate = 30;
final_rate = 3;
down = orig_rate./final_rate;
nON = 150./down;
nOFF = 150./down;
pre_win = [floor(0.75*nOFF) nOFF];
post_win = [nOFF+1 nOFF+nON];

dirs = unique(cell2mat(input.tGratingDirectionDeg));
nCond = length(dirs);

alpha = .05./(nCond);
b= 5;
pix = 5;
bouton_diam = 5;

% Az = [0 15 30];
% El = [0];

date = '140209';
mouse = 'G008';

%% load data
data = readrawfile;

%if two-channel, choose just green channel
data = data(:,:,:,1);

%% reshape data
%average signals in time

data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down
%% register
data_avg = mean(data_sub(:,:,1:100),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

%% sort data
nRep_mat = zeros(1,nCond);
siz = size(data_reg);
data_sort = zeros(siz);
start = 1;
for iCond = 1:nCond
    ids = find(cell2mat(input.tGratingDirectionDeg) == dirs(iCond));
     nRep_mat(1,iCond) = length(ids);
    for iRep = 1:length(ids)
        trial = ids(iRep);
        data_sort(:,:,start:start+nON+nOFF-1) = data_reg(:,:,1+((trial-1)*(nON+nOFF)):trial*(nON+nOFF));
        start = start+nON+nOFF;
    end
end

%% create resp structure
resp = struct;
start = 1;
for iCond = 1:nCond; 
    resp(iCond).on = [];
    resp(iCond).off = [];
    roi_stim = zeros(siz(1),siz(2),nOFF+nON,nRep_mat(1,iCond));
    for iRep = 1:nRep_mat(1,iCond);
        roi_stim(:,:,:,iRep) = data_sort(:,:,start:start-1+nON+nOFF);
        start = start+nON+nOFF;
    end
    resp(iCond).off = squeeze(mean(roi_stim(:,:,pre_win(1):pre_win(2),:),3));
    resp(iCond).on = squeeze(mean(roi_stim(:,:,post_win(1):post_win(2),:),3));
end

%% pixel based ttest
Info_ttest_mat = zeros(siz(1),siz(2),nCond);

f1 = fspecial('average');
for iCond = 1:nCond
    resp(iCond).on_long = reshape(resp(iCond).on, [siz(1) siz(2)*nRep_mat(1,iCond)]);
    resp(iCond).off_long = reshape(resp(iCond).off, [siz(1) siz(2)*nRep_mat(1,iCond)]);
    resp(iCond).on_sm = reshape(filter2(f1,resp(iCond).on_long),[siz(1) siz(2) nRep_mat(1,iCond)]);
    resp(iCond).off_sm = reshape(filter2(f1,resp(iCond).off_long),[siz(1) siz(2) nRep_mat(1,iCond)]);
end

for iy = b+1:siz(1)-b
    fprintf([num2str(iy) ' '])
    for ix = b+1:siz(2)-b
        p_ttestB = zeros(1,1,nCond);
        for iCond = 1:nCond;
            [h_ttestB1,p_ttestB1] = ttest(resp(iCond).off_sm(iy,ix,:),resp(iCond).on_sm(iy,ix,:),alpha,'left');
            p_ttestB(1,1,iCond) = p_ttestB1;
        end
    Info_ttest_mat(iy,ix,:) = p_ttestB;
    end
end

Info_ttest_mat_long = interp2(reshape(Info_ttest_mat, [siz(1) siz(2)*nCond]));
Info_ttest_smooth = filter2(f1,Info_ttest_mat_long);
siz_interp = size(Info_ttest_smooth);
Xi = 1:2:siz_interp(1);
Yi = 1:2:siz_interp(2);
ttest_smooth_siz = interp2(Info_ttest_smooth, Yi', Xi);
ttest_smooth = min(reshape(ttest_smooth_siz, [siz(1) siz(2) nCond]),[],3);
ttest_mask = min(ttest_smooth,[],3) < alpha;

ttest_mask(1:b,1:end) = 0;
ttest_mask(1:end, 1:b) = 0;
ttest_mask(1:end, siz(2)-b:end) = 0;
ttest_mask(siz(1)-b:end,1:end) = 0;

figure; imagesq(ttest_mask);
        
%% find dF/F for iCond to find active boutons
resp_dFoverF = zeros(siz(1),siz(2),nCond);
for iCond = 1:nCond
    resp_dFoverF(:,:,iCond) = (mean(resp(iCond).on,3)-mean(resp(iCond).off,3))./(mean(resp(iCond).off,3));
end

max_dF = max(resp_dFoverF,[],3);
figure; imagesq(max_dF); colormap(gray)

%% use max dF/F to find ROIS- local maxima

max_interp = interp2(max_dF);
f1 = fspecial('average');   
max_interp_sm = filter2(f1, max_interp);
siz2 = size(max_interp);
Xi = 1:2:siz2(1);
Yi = 1:2:siz2(2);
stack_max_interp_sm = interp2(max_interp_sm, Yi', Xi);

local_max = zeros(siz(1), siz(2));
for iy = b+1:(siz(1)-b);
    for ix = b+1:(siz(2)-b);            
        sub = stack_max_interp_sm(iy-pix:iy+pix,ix-pix:ix+pix);
        sub_long = reshape(sub, [1 (pix*2+1)^2]);
        [sub_long_order ind_order] = sort(sub_long);
        if ind_order(end)==ceil(((pix*2+1)^2)/2)
            local_max(iy,ix) = 1;
        end
    end
end
ind_local_max = find(local_max == 1);
figure; imagesq(local_max);

%% combine ttest and local maxima
ttest_long = reshape(ttest_smooth, [siz(1)*siz(2) 1]);
ind_highP = find(ttest_long(ind_local_max,:)>=(alpha));
local_max_long(ind_local_max(ind_highP,:),:) = 0;
local_max_sig = reshape(local_max_long, [siz(1) siz(2)]);

n_pix = sum(sum(local_max_sig));
[i, j] = find(local_max_sig ==1);

%expand local maxima
pix_add = floor(bouton_diam/2);
FOV = zeros(size(local_max_sig));
for ipix = 1:n_pix
    FOV(i(ipix)-pix_add:i(ipix)+pix_add,j(ipix)-pix_add:j(ipix)+pix_add) = 1;
end

figure; imagesq(FOV)

%% Get timecourses

boutons = bwlabel(FOV);
figure; imagesq(boutons);
data_TC = stackGetTimeCourses(data_reg,boutons);

figure; tcOffsetPlot(data_TC)

%% plot responses
stim_mat = zeros(nCond,nRep,length(pre_win(1):post_win(2)));
for iRep = 1:nRep
    for iCond = 1:nCond 
        stim_mat(iCond,iRep,:) = pre_win(1)+((iCond-1)*(nON+nOFF))+((iRep-1)*((nON+nOFF)*nCond)): nOFF+nON +((iCond-1)*(nON+nOFF))+((iRep-1)*((nON+nOFF)*nCond));
    end
end

%for retinotopy
for iCell = 7:15;
    figure;
    for iCond = 1:nCond
        subplot(1,nCond,iCond)
        rep_mat = zeros(length(pre_win(1):post_win(2)),nRep);
        for iRep = 1:nRep
            plot(pre_win(1)-post_win(1):post_win(2)-post_win(1),data_TC(squeeze(stim_mat(iCond,iRep,:))',iCell), 'k');
            hold on
            rep_mat(:,iRep) = data_TC(squeeze(stim_mat(iCond,iRep,:))',iCell);
        end
        plot(pre_win(1)-post_win(1):post_win(2)-post_win(1), mean(rep_mat,2), 'r');
        hold on
        ylim([min(data_TC(:,iCell),[],1) max(data_TC(:,iCell),[],1)])
        stim_Az = rem(iCond,size(Az,2));
        if stim_Az == 0
            stim_Az = size(Az,2);
        end
        stim_El= ceil(iCond./size(Az,2));
        title(['Az = ' num2str(Az(1,stim_Az)) ' El = ' num2str(El(1,stim_El))]);
    end
end

%for speed
for iCell = 7:15;
    figure;
    for iCond = 1:nCond
        subplot(1,nCond,iCond)
        rep_mat = zeros(length(pre_win(1):post_win(2)),nRep);
        for iRep = 1:nRep
            plot(pre_win(1)-post_win(1):post_win(2)-post_win(1),data_TC(squeeze(stim_mat(iCond,iRep,:))',iCell), 'k');
            hold on
            rep_mat(:,iRep) = data_TC(squeeze(stim_mat(iCond,iRep,:))',iCell);
        end
        plot(pre_win(1)-post_win(1):post_win(2)-post_win(1), mean(rep_mat,2), 'r');
        hold on
        ylim([min(data_TC(:,iCell),[],1) max(data_TC(:,iCell),[],1)])
        stim_TF = TF(iCond);
        stim_SF= SF(iCond);
        stim_speed = stim_TF./stim_SF;
        title(['Speed = ' num2str(stim_speed)]);
    end
end
