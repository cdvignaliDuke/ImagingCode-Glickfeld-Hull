
% nON = #frames of stimulus
% nOFF = #frames of ISI
% nCond = #unique stimuli
% nRep_mat = [1 nCond] matrix giving #repeats per stimulus type (iCond)
% alpha = 0.05/nCond % for pixel-based ttest- may need to tune stringency 
% b= 5; % sets size of border to compensate for registration jitter
% pix = 5; % sets size of area to search for local maxima for boutons and
% minimum distance between boutons- depends on your resolution
% bouton_diam = 5; %sets size of square ROI to be averaged for timecourse


%start with a motion-registered, condition-sorted tiff stack "data"
siz = size(data);

% create resp structure where for each stimulus type (iCond) you have an
% average fov ([siz(1) siz(2)]) for every stimulus repeat for "on" and
% "off"
resp = struct;
% resp(iCond).on = [siz(1) siz(2) nRep_mat(1,iCond)];
% resp(iCond).off = [siz(1) siz(2) nRep_mat(1,iCond)];

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

Info_ttest_mat_long = reshape(Info_ttest_mat, [siz(1) siz(2)*nCond]);
ttest_smooth = reshape(filter2(f1,Info_ttest_mat_long), [siz(1) siz(2) nCond]);
ttest_mask = min(ttest_smooth,[],3) < alpha;

ttest_mask(1:b,1:end) = 0;
ttest_mask(1:end, 1:b) = 0;
ttest_mask(1:end, siz(2)-b:end) = 0;
ttest_mask(siz(1)-b:end,1:end) = 0;

figure; imagesq(ttest_mask);

%% Find dF/F for each condition to find active boutons

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
local_max_long = reshape(local_max, [siz(1)*siz(2) 1]);
ind_local_max = find(local_max_long==1);
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

        
