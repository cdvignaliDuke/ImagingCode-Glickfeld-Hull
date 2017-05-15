clear all
close all
rc = behavConstsAV;
awFSAVdatasets_V1axonsPM
for iexp = 1:size(expt,2)

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
dirFolder = expt(iexp).dirtuning;
dirTime = expt(iexp).dirtuning_time;
down = 10;

fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate, dirFolder);
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate);
%% load tuning data
fName = [dirFolder '_000_000'];
[input, data] = Load_SBXdataPlusMWorksData(SubNum,expDate,dirTime,mouse,dirFolder,fName);  

% down-sample
data_down = stackGroupProject(data,down);
clear data

% remove negative data by subtraction
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down


%% load outs and re-regester data
load(fullfile(fn,'regOuts&Img.mat'));

[out data_reg] = stackRegister_MA(data_sub,[],[],out_tun);

xpix = size(data_reg,2);
ypix = size(data_reg,1);
%% output goes here
fnout = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate);

%% sort data by trial type
[resp,data_tr,nstim] = sortDirectionData(data_reg,down,input);
nRep = cellfun(@(x) size(x,3),{resp.on});
%% pixel-based ttest
alpha = .05./(nstim);
b= 5;
pix = 5;
bouton_diam = 5; % 5 for 16x, 8 for 25x

f1 = fspecial('average');
for istim = 1:nstim
    resp(istim).on_long = reshape(resp(istim).on, [ypix xpix*nRep(istim)]);
    resp(istim).off_long = reshape(resp(istim).off, [ypix xpix*nRep(istim)]);
    resp(istim).on_sm = reshape(filter2(f1,resp(istim).on_long),[ypix xpix nRep(istim)]);
    resp(istim).off_sm = reshape(filter2(f1,resp(istim).off_long),[ypix xpix nRep(istim)]);
end


Info_ttest_mat = zeros(ypix,xpix,nstim);
for iy = b+1:ypix-b
    fprintf([num2str(iy) ' '])
    for ix = b+1:xpix-b
        p_ttestB = zeros(1,1,nstim);
        for istim = 1:nstim;
            [h_ttestB1,p_ttestB1] = ttest(resp(istim).off_sm(iy,ix,:),resp(istim).on_sm(iy,ix,:),alpha,'left');
            p_ttestB(1,1,istim) = p_ttestB1;
        end
    Info_ttest_mat(iy,ix,:) = p_ttestB;
    end
end

Info_ttest_mat_long = interp2(reshape(Info_ttest_mat, [ypix xpix*nstim]));
Info_ttest_smooth = filter2(f1,Info_ttest_mat_long);
siz_interp = size(Info_ttest_smooth);
Xi = 1:2:siz_interp(1);
Yi = 1:2:siz_interp(2);
ttest_smooth_siz = interp2(Info_ttest_smooth, Yi', Xi);
ttest_smooth = min(reshape(ttest_smooth_siz, [ypix xpix nstim]),[],3);
ttest_mask = min(ttest_smooth,[],3) < alpha;

ttest_mask(1:b,1:end) = 0;
ttest_mask(1:end, 1:b) = 0;
ttest_mask(1:end, xpix-b:end) = 0;
ttest_mask(ypix-b:end,1:end) = 0;

figure; imagesq(ttest_mask);

%% find local activity maxima to get ROIs

resp_dff = cellfun(@(x,y) mean(y,3)-mean(x,3)./mean(x,3), {resp.off},{resp.on},'unif',false);
resp_dff = reshape(cell2mat(resp_dff),ypix,xpix,nstim);
max_dff = max(resp_dff,[],3);

max_interp = interp2(max_dff);
f1 = fspecial('average');   
max_interp_sm = filter2(f1, max_interp);
siz2 = size(max_interp);
Xi = 1:2:siz2(1);
Yi = 1:2:siz2(2);
stack_max_interp_sm = interp2(max_interp_sm, Yi', Xi);

local_max = zeros(ypix, xpix);
for iy = b+1:(ypix-b);
    for ix = b+1:(xpix-b);            
        sub = stack_max_interp_sm(iy-pix:iy+pix,ix-pix:ix+pix);
        sub_long = reshape(sub, [1 (pix*2+1)^2]);
        [sub_long_order, ind_order] = sort(sub_long);
        if ind_order(end)==ceil(((pix*2+1)^2)/2)
            local_max(iy,ix) = 1;
        end
    end
end
ind_local_max = find(local_max == 1);
figure; imagesq(local_max);

% combine ttest and local maxima
ttest_long = reshape(ttest_smooth, [ypix*xpix 1]);
ind_highP = find(ttest_long(ind_local_max,:)>=(alpha));
local_max_long = reshape(local_max, [ypix*xpix 1]);
local_max_long(ind_local_max(ind_highP,:),:) = 0;
local_max_sig = reshape(local_max_long, [ypix xpix]);

n_pix = sum(sum(local_max_sig));
[i, j] = find(local_max_sig ==1);

%expand local maxima
pix_add = floor(bouton_diam/2);
FOV = zeros(size(local_max_sig));
for ipix = 1:n_pix
    FOV(i(ipix)-pix_add:i(ipix)+pix_add,j(ipix)-pix_add:j(ipix)+pix_add) = 1;
end

figure; imagesq(FOV);

mask_boutons = bwlabel(FOV);

save(fullfile(fnout,'final_mask.mat'),'mask_boutons');

end