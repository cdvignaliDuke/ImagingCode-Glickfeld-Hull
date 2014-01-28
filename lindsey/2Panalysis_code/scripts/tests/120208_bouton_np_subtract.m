P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120208';
mouse = 'AC45';
date = '110820';
userun = [1:4];
nCond = 50;
dirs = 2;

base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse,date);



fn_stim = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_allstim.mat']);
load(fn_stim);

fn_reps = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
load(fn_reps);

fn_resp = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
load(fn_resp);

bouton_map = zeros(size(local_max));
for ipix = 1:n_pix
    sub_y = [i(ipix)-1:i(ipix)+1]; 
    sub_x = [j(ipix)-1:j(ipix)+1];
    bouton_map(sub_y, sub_x) = ipix;
end

siz = size(bouton_map); 
npmasks = zeros(siz(1),siz(2),n_pix);
mask_dil = zeros(siz(1),siz(2),n_pix);
for ipix = 1:n_pix
    mask_dil(i(ipix)-2:i(ipix)+2, j(ipix)-2:j(ipix)+2,ipix) = 1;
end
all_mask_dil = sum(mask_dil,3);

for ipix = 1:n_pix
    np_plus = zeros(siz(1),siz(2));
    np_plus(i(ipix)-4:i(ipix)+4, j(ipix)-4:j(ipix)+4) = 1;
	npmask = np_plus-mask_dil(:,:,ipix);
    npmask(find(all_mask_dil>0)) = 0;
    npmask(1:5,:)=0;
    npmask(236:240,:)=0;
    npmask(:,1:5)=0;
    npmask(:,252:256)=0;
    npmasks(:,:,ipix) = npmask;
end

resp_dF_np = zeros(n_pix, size(stim_on,3));
resp_dFoverF_np = zeros(n_pix, size(stim_on,3));

siz = size(stim_on); 
stim_on_long = reshape(stim_on, siz(1)*siz(2), siz(3));
stim_off_long = reshape(stim_off, siz(1)*siz(2), siz(3));
npmask_sum = sum(npmasks,3);
npmask_all = zeros(size(npmask_sum));
npmask_all(find(npmask_sum>0))=1;
ind_all = find(npmask_all==1);
stim_on_npmask_all = mean(stim_on_long(ind_all,:),1);    
stim_off_npmask_all = mean(stim_off_long(ind_all,:),1);
nponly_dFoverF = (stim_on_npmask_all-stim_off_npmask_all)./stim_off_npmask_all;

np_all_avg = mean([stim_on_npmask_all stim_off_npmask_all],2);
roi_on_np = zeros(n_pix, sum(stim_reps));
roi_off_np = zeros(n_pix, sum(stim_reps));
for ipix = 1:n_pix
    sub_y = [i(ipix)-1:i(ipix)+1]; 
    sub_x = [j(ipix)-1:j(ipix)+1];
    stim_on_long = reshape(stim_on, siz(1)*siz(2),siz(3));
    stim_off_long = reshape(stim_off, siz(1)*siz(2),siz(3));
    ind = find(npmasks(:,:,ipix)==1);
    stim_on_npmask = mean(stim_on_long(ind,:),1);    
    stim_off_npmask = mean(stim_off_long(ind,:),1);
    np_avg = mean([stim_on_npmask stim_off_npmask],2);
    np_on = stim_on_npmask_all*(np_avg/np_all_avg);
    np_off = stim_off_npmask_all*(np_avg/np_all_avg);
    roi_on_np(ipix,:) = bsxfun(@minus, squeeze(mean(mean(stim_on(sub_y,sub_x,:),1),2)), np_on') +np_avg;
    roi_off_np(ipix,:) = bsxfun(@minus, squeeze(mean(mean(stim_off(sub_y,sub_x,:),1),2)), np_off') +np_avg;
    resp_dF_np(ipix,:) = (roi_on_np(ipix,:)-roi_off_np(ipix,:))';
    resp_dFoverF_np(ipix,:) = ((roi_on_np(ipix,:)-roi_off_np(ipix,:))./roi_off_np(ipix,:))';
end


nReps = sum(stim_reps(1,:));
if dirs == 2
    add = 1;
    stim_reps_dir = zeros(1,26);
    for iCond = 1:nCond/2
    stim_reps_dir(1,iCond) = sum(stim_reps(1,add:add+1));
    add = add+2;
    end
    stim_reps_dir(1,end) = stim_reps(1,end);
    stim_reps = stim_reps_dir;
end        

Ind_struct = [];
start = 1;
for iCond = 1:25
    nRep = stim_reps(iCond);
    Ind_struct(iCond).all_trials = [start:start-1+nRep];
    start = start+nRep;
end

Im_mat_nponly = zeros(1, 25);
Im_mat_np = zeros(n_pix, 25);
Im_mat_USE = zeros(n_pix, 25);
for iCond = 1:25        
    ind_all = Ind_struct(iCond).all_trials;
    Im_mat_nponly(:,iCond) = mean(nponly_dFoverF(:,ind_all),2);
    Im_mat_np(:,iCond) = mean(resp_dFoverF_np(:,ind_all),2);
    Im_mat_USE(:,iCond) = mean(resp_dFoverF(:,ind_all),2);
end

Fit_struct_np2 = [];
for count_shuf = 0:nShuf
    fprintf('.')
    Im_mat_USE = zeros(size(resp_dFoverF,1), 25);
    Im_mat_std = zeros(size(resp_dFoverF,1), 25);
    dF_mat = zeros(size(resp_dFoverF,1), 25);
    for iCond = 1:25        
        ind_all = Ind_struct(iCond).all_trials;
        if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
            ind_all_1 = ind_all(randsample(length(ind_all),length(ind_all),1));
        else
            ind_all_1 = ind_all;        
        end
        Im_mat_USE(:,iCond) = mean(resp_dFoverF_np2(:,ind_all_1),2);
        dF_mat(:,iCond) = mean(resp_dF_np2(:,ind_all_1),2);
    end

    start = 1;
    for iCell = 1:n_pix;
        a = Im_mat_USE(iCell,:);
        if max(a,[],2) > 0     
            b = reshape(a',length(SF_vec0),length(TF_vec0));
            %b2 = b( ind_SFuse(:,1),ind_TFuse(:,1));
            data = b';
            ind0 = find(data<0);
            data(ind0) = NaN;
            if count_shuf == 0
                PLOTIT_FIT = 0;
                SAVEALLDATA = 1;
                Fit_2Dellipse_LG
                eval(['Fit_struct_np2(iCell).True.s_',' = s;']);
            else
                SAVEALLDATA = 0;
                PLOTIT_FIT = 0;
                Fit_2Dellipse_LG
                eval(['Fit_struct_np2(iCell).Shuf(count_shuf).s_',' = s;']);
            end
        end               
    end
end

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct_np2.mat']);   
save(fn_out, 'Fit_struct_np2')
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct.mat']);   
        load(fn_out)
        if Nshuf>1;
        for iCell = 1:n_pix
            if ~isempty(Fit_struct_np(iCell).True)                
                eval(['tmp = Fit_struct_np2(iCell).True.s_.x;']);
                eval(['tmp = [tmp Fit_struct_np2(iCell).True.s_.SFhicut_50];']);
                eval(['tmp = [tmp Fit_struct_np2(iCell).True.s_.TFhicut_50];']);
                eval(['tmp = [tmp Fit_struct_np2(iCell).True.s_.SFhicut_10];']);
                eval(['tmp = [tmp Fit_struct_np2(iCell).True.s_.TFhicut_10];']);
                fit_true_vec(iCell,:) = tmp;
            end
        end

        for count_shuf = 1:Nshuf
            for iCell = 1:n_pix
                if ~isempty(Fit_struct_np2(iCell).Shuf)
                    eval(['tmp = Fit_struct_np2(iCell).Shuf(count_shuf).s_.x;']);
                    eval(['tmp = [tmp Fit_struct_np2(iCell).Shuf(count_shuf).s_.SFhicut_50];']);
                    eval(['tmp = [tmp Fit_struct_np2(iCell).Shuf(count_shuf).s_.TFhicut_50];']);
                    eval(['tmp = [tmp Fit_struct_np2(iCell).Shuf(count_shuf).s_.SFhicut_10];']);
                    eval(['tmp = [tmp Fit_struct_np2(iCell).Shuf(count_shuf).s_.TFhicut_10];']);
                    %fit is: %A sigma_SF sigma_TF sf0 tf0 xi
                    fit_shuf_vec(iCell,:,count_shuf) = tmp;
                end
            end
        end

        Npars = size(fit_shuf_vec,2);
        lbub_fits_np2 = zeros(n_pix,Npars,5);
        alpha_bound = .025;
        for iCell = 1:n_pix
            for count2 = 1:Npars
                tmp = squeeze(fit_shuf_vec(iCell,count2,:));
                [i,j] = sort(tmp);
                ind_shuf_lb = ceil(Nshuf*alpha_bound);
                ind_shuf_ub = ceil(Nshuf*(1-alpha_bound));
                lbub_fits_np2(iCell,count2,1) = i(ind_shuf_lb);
                lbub_fits_np2(iCell,count2,2) = i(ind_shuf_ub);
                lbub_fits_np2(iCell,count2,3) = mean(i); 
                lbub_fits_np2(iCell,count2,5) = std(i);
            end
            %now take means from truedata fit:
            lbub_fits_np2(iCell,:,4) = fit_true_vec(iCell,:);
        end
    end

    lbub_diff_np2 = lbub_fits_np2(:,:,2)-lbub_fits_np2(:,:,1);

    goodfit_ind_np2 = [];
    for iCell = 1:n_pix
        if lbub_diff_np2(iCell,4)<2 
            if lbub_diff_np2(iCell,5)<2
                goodfit_ind_np2 = [goodfit_ind_np2 iCell];
            end
        end
    end

    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits_np.mat']);   
    save(fn_out, 'lbub_fits_np', 'lbub_diff_np', 'goodfit_ind_np')
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']); 
    load(fn_out);
figure;
start = 1;
for ipix = 1:100
   if start >64
       figure;
       start = 1;
   end
   if find(goodfit_ind_np==ipix)>0 
       data = reshape(Im_mat_USE(ipix,:), 5,5)';
       data_np2 = reshape(Im_mat_np2(ipix,:), 5,5)';
       subplot(8,8,start)
       imagesq(data);
       caxis([0 1])
       subplot(8,8,start+1)
       imagesq(data_np2);
       caxis([0 1])
       start = start+2;
       colormap(gray)
   end
end

figure;
start = 1;
for ipix = 1:n_pix
   if start >64
       figure;
       start = 1;
   end
   subplot(8,8,start)
   imagesq(Fit_struct(ipix).True.s_.k2b_plot);
   if find(goodfit_ind==ipix)>0 
       title('**')
   end
   subplot(8,8,start+1)
   imagesq(Fit_struct_np2(ipix).True.s_.k2b_plot);
   if find(goodfit_ind_np2==ipix)>0 
       title('**')
   end
   start = start+2;
   colormap(gray)
end

plotfit = [];
plotfit_np2 = [];
for ipix = 1:n_pix
    if find(goodfit_ind==ipix)>0 
        plotfit = cat(3, plotfit, Fit_struct(ipix).True.s_.k2b_plot);
    end
    if find(goodfit_ind_np2==ipix)>0 
        plotfit_np2 = cat(3, plotfit_np2, Fit_struct_np2(ipix).True.s_.k2b_plot);
    end
end
plotfit_av = mean(plotfit,3);
plotfit_np2_av = mean(plotfit_np2,3);
figure;
subplot(1,2,1)
imagesq(plotfit_av)
title('No correction')
subplot(1,2,2)
imagesq(plotfit_np2_av)
title('NP subtracted')
colormap(gray)

ind = intersect(goodfit_ind, goodfit_ind_np2);
TF_comp = zeros(length(ind),2);
SF_comp = zeros(length(ind),2);
dF_comp = zeros(length(ind),2);
for ipix = 1:length(ind);
    TF_comp(ipix,1) = 2^Fit_struct(ind(1,ipix)).True.s_.x(5);
    TF_comp(ipix,2) = 2^Fit_struct_np2(ind(1,ipix)).True.s_.x(5);
    SF_comp(ipix,1) = 2^Fit_struct(ind(1,ipix)).True.s_.x(4);
    SF_comp(ipix,2) = 2^Fit_struct_np2(ind(1,ipix)).True.s_.x(4);
    dF_comp(ipix,1) = Fit_struct(ind(1,ipix)).True.s_.x(1);
    dF_comp(ipix,2) = Fit_struct_np2(ind(1,ipix)).True.s_.x(1);
end
speed_comp = zeros(length(ind),2);
speed_comp(:,1) = TF_comp(:,1)./SF_comp(:,1);
speed_comp(:,2) = TF_comp(:,2)./SF_comp(:,2);

figure;
subplot(2,2,1)
semilogy(1:2, TF_comp,'-c');
hold on
errorbar(1:2, mean(TF_comp,1), std(TF_comp,[],1)./sqrt(length(ind)), '-sk', 'LineWidth', 3);
title('TF');
xlim([0 3])
ylim([0 20])
subplot(2,2,2)
semilogy(1:2, SF_comp,'-c');
hold on
errorbar(1:2, mean(SF_comp,1), std(SF_comp,[],1)./sqrt(length(ind)), '-sk', 'LineWidth', 3);
title('SF');
xlim([0 3])
subplot(2,2,3)
plot(1:2, dF_comp,'-c');
hold on
errorbar(1:2, mean(dF_comp,1), std(dF_comp,[],1)./sqrt(length(ind)), '-sk', 'LineWidth', 3);
title('dF');
xlim([0 3])
subplot(2,2,4)
semilogy(1:2, speed_comp,'-c');
hold on
errorbar(1:2, mean(speed_comp,1), std(speed_comp,[],1)./sqrt(length(ind)), '-sk', 'LineWidth', 3);
title('speed');
xlim([0 3])



all_TF = [];
all_TF_np2 = [];
all_SF = [];
all_SF_np2 = [];
all_dF = [];
all_dF_np2 =[];
all_TF_only = [];
all_SF_only = [];
all_dF_only = [];
for ipix = 1:n_pix
    if find(goodfit_ind==ipix)>0
        all_TF = [all_TF; 2^Fit_struct(ipix).True.s_.x(5)];
        all_SF = [all_SF; 2^Fit_struct(ipix).True.s_.x(4)];
        all_dF = [all_dF; Fit_struct(ipix).True.s_.x(1)];
        if find(setdiff(goodfit_ind, goodfit_ind_np2)==ipix)
            all_TF_only = [all_TF_only; 2^Fit_struct(ipix).True.s_.x(5)];
            all_SF_only = [all_SF_only; 2^Fit_struct(ipix).True.s_.x(4)];
            all_dF_only = [all_dF_only; Fit_struct(ipix).True.s_.x(1)];
        end
    end
    if find(goodfit_ind_np2==ipix)>0
        all_TF_np2 = [all_TF_np2; 2^Fit_struct_np2(ipix).True.s_.x(5)];
        all_SF_np2 = [all_SF_np2; 2^Fit_struct_np2(ipix).True.s_.x(4)];
        all_dF_np2 = [all_dF_np2; Fit_struct_np2(ipix).True.s_.x(1)];
    end
end
all_speed = all_TF./all_SF;
all_speed_np2 = all_TF_np2./all_SF_np2;
all_speed_only = all_TF_only./all_SF_only;

figure;
subplot(2,2,1)
[H, stats, xCDF, yCDF] = cdfplot_LG(all_TF);
[H_only, stats_only, xCDF_only, yCDF_only] = cdfplot_LG(all_TF_only);
[H_np2, stats_np2, xCDF_np2, yCDF_np2] = cdfplot_LG(all_TF_np2);
semilogx(xCDF, yCDF, 'k');
hold on
plot(xCDF_np2, yCDF_np2, 'r');
hold on
plot(xCDF_only, yCDF_only, 'b');
hold on
title('TF');
xlim([0 15])
subplot(2,2,2)
[H, stats, xCDF, yCDF] = cdfplot_LG(all_SF);
[H_only, stats_only, xCDF_only, yCDF_only] = cdfplot_LG(all_SF_only);
[H_np2, stats_np2, xCDF_np2, yCDF_np2] = cdfplot_LG(all_SF_np2);
semilogx(xCDF, yCDF, 'k');
hold on
plot(xCDF_np2, yCDF_np2, 'r');
hold on
plot(xCDF_only, yCDF_only, 'b');
title('SF');
xlim([0.02 .32])
subplot(2,2,3)
[H, stats, xCDF, yCDF] = cdfplot_LG(all_dF);
[H_only, stats_only, xCDF_only, yCDF_only] = cdfplot_LG(all_dF_only);
[H_np2, stats_np2, xCDF_np2, yCDF_np2] = cdfplot_LG(all_dF_np2);
plot(xCDF, yCDF, 'k');
hold on
plot(xCDF_np2, yCDF_np2, 'r');
hold on
plot(xCDF_only, yCDF_only, 'b');
title('dF');
xlim([0 1])
subplot(2,2,4)
[H, stats, xCDF, yCDF] = cdfplot_LG(all_speed);
[H_only, stats_only, xCDF_only, yCDF_only] = cdfplot_LG(all_speed_only);
[H_np2, stats_np2, xCDF_np2, yCDF_np2] = cdfplot_LG(all_speed_np2);
semilogx(xCDF, yCDF, 'k');
hold on
plot(xCDF_np2, yCDF_np2, 'r');
hold on
plot(xCDF_only, yCDF_only, 'b');
title('speed');
xlim([0 800])

figure;
subplot(10,10,1)
imagesq(reshape(Im_mat_nponly, [5 5])');
caxis([0 .8]);
start =2;
for ipix = 1:177
    if find(goodfit_ind==ipix)>0
        subplot(10,10,start)
        imagesq(reshape(Im_mat_np(ipix,:), [5 5])');
        caxis([0 .8]);
        start = start+1;
    end
end
colormap(gray)
    
anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120208';
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_exampleY13_lost_boutons_neuropilsub2']);
         print(gcf, '-dpdf', fn_out);
         

