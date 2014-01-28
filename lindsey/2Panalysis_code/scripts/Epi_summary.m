P = 1;
matrix = 'SF5xTF5';
image = 'intrinsic_10Hz';
inj = 'all';
nCond = 25;
if nCond == 9;
    oversamp = 33;
elseif nCond == 25;
    oversamp = 65;
end
nON=10;
nOFF=40;
post_win = [11 20];
sum_base = 'G:\users\lindsey\analysisLG\experiments';
base = 'G:\users\lindsey\analysisLG\active mice';
anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120203';
str_run = strvcat('allrun', 'running', 'norun');

list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
load(list_fn);
nexp = size(exp_list.mouse_mat,2);
mouse_list = [];
mice = exp_list.mouse_mat;

resp_avg_all = zeros(8, nCond+1, 3, nexp);
TC_avg_all = zeros(nON+nOFF, 8, nCond+1, 3, nexp);
% CM_mat_all = zeros(8, 2, 3, nexp);
H_ttest_all = zeros(8, 3, nexp);
Fit_all = zeros(8,6,3,nexp);
Fit_oversamp_all = zeros(8,oversamp,oversamp,3,nexp);
% F_all = zeros(8,2,nexp);

for iexp = 1:nexp
    mouse = char(exp_list.mouse_mat{iexp});
    date = char(exp_list.date_mat{iexp});
    userun = exp_list.runs_mat{iexp};
    count_prot = exp_list.prot_mat{iexp};
    run = exp_list.run_mat{iexp};
    blanks = exp_list.blanks_mat{iexp};
    
    mouse_list = [mouse_list '_' mouse];
    outDir = fullfile(base, mouse,date);
    
    fn_areas = fullfile(outDir,'analysis', [date '_' mouse '_run' num2str(userun) '_area_mask.mat']);
    load(fn_areas);
    fn_resp = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp_POST' num2str(post_win) '.mat']);
    load(fn_resp);
    fn_resp_norm = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp_norm_POST' num2str(post_win) '.mat']);
    load(fn_resp_norm);
    fn_var = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_vars.mat']);
    load(fn_var);
    fn_ttest = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_POST' num2str(post_win)  '_ttest.mat']);
    load(fn_ttest);
    fn_fit = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_POST' num2str(post_win) '_Fit_struct.mat']);
    load(fn_fit);
%     fn_F = fullfile(outDir,'analysis', [date '_' mouse '_run' num2str(userun) '_F.mat']);
%     load(fn_F);
    
    areas = size(area_list,1);

    a = 'V1'; b = 'PM'; c = 'LM'; d = 'AM'; e = 'AL'; f = 'RL'; g = 'A '; h = 'P ';
    all_areas = strvcat(a,b,c,d,e,f,g,h);
        col = ['g' 'c' 'k' 'b' 'r' 'm' 'y' 'y'];
    
    ind = zeros(1,length(all_areas));
    for iArea = 1:areas;
        if area_list(iArea,:) == 'V1'
            ind(:,1) = iArea;
        elseif area_list(iArea,:) == 'PM'
            ind(:,2) = iArea;
        elseif area_list(iArea,:) == 'LM'
            ind(:,3) = iArea;
        elseif area_list(iArea,:) == 'AM'
            ind(:,4) = iArea;
        elseif area_list(iArea,:) == 'AL'
            ind(:,5) = iArea;
        elseif area_list(iArea,:) == 'RL'
            ind(:,6) = iArea;
        elseif area_list(iArea,:) == 'A '
            ind(:,7) = iArea;
        elseif area_list(iArea,:) == 'P '
            ind(:,8) = iArea;   
        end
    end

    resp_avg_reorder = zeros(length(all_areas), nCond+1,3);
    TC_avg_reorder = zeros(nON+nOFF, length(all_areas), nCond+1,3);
%     CM_mat_reorder = zeros(length(all_areas), 2 ,3);
    H_ttest_reorder = zeros(length(all_areas), 3);
    Fit_reorder = zeros(length(all_areas),6,3);
    Fit_oversamp_reorder = zeros(length(all_areas),oversamp,oversamp,3);
%     F_reorder = zeros(length(all_areas),2);
    
    for iind = 1:length(all_areas);
        area = ind(:,iind);
        if area>0
            H_ttest_reorder(iind,:) = H_ttest(area,:);
%             F_reorder(iind,1) = area_F(:,area);
%             F_reorder(iind,2) = np_F(:,area);
            if blanks == 1;
                if run ==1;
                    resp_avg_reorder(iind,:,:) = resp_avg(area,1:nCond+1,:);
                    TC_avg_reorder(:,iind,:,:) = TC_avg(:,area,1:nCond+1,:);
                    for iRun = 1:3
%                         eval(['CM_mat_reorder(iind,:,iRun) = Fit_struct(area).True.s_' str_run(iRun,:) '.CM_data;']);                        
                        eval(['Fit_reorder(iind,:,iRun) = Fit_struct(area).True.s_' str_run(iRun,:) '.x;']);                        
                        eval(['Fit_oversamp_reorder(iind,:,:,iRun) = Fit_struct(area).True.s_' str_run(iRun,:) '.k2_plot_oversamp;']);                        
                    end
                else
                    resp_avg_reorder(iind,:,1) = resp_avg(area,1:nCond+1);
                    TC_avg_reorder(:,iind,:,1) = TC_avg(:,area,1:nCond+1);
%                     CM_mat_reorder(iind,:,1) = Fit_struct(area).True.s_allrun.CM_data';
                    Fit_reorder(iind,:,1) = Fit_struct(area).True.s_allrun.x';       
                    Fit_oversamp_reorder(iind,:,:,1) = Fit_struct(area).True.s_allrun.k2_plot_oversamp;
                end
            else
                if run == 1;
                    resp_avg_reorder(iind,1:nCond,:) = resp_avg(area,1:nCond,:);
                    resp_avg_reorder(iind,nCond+1,:) = NaN;
                    TC_avg_reorder(:,iind,1:nCond,:) = TC_avg(:,area,1:nCond,:);
                    TC_avg_reorder(:,iind,nCond+1,:) = NaN;
                    for iRun = 1:3
%                         eval(['CM_mat_reorder(iind,:,iRun) = Fit_struct(area).True.s_' str_run(iRun,:) '.CM_data;']);                        
                        eval(['Fit_reorder(iind,:,iRun) = Fit_struct(area).True.s_' str_run(iRun,:) '.x;']);                        
                        eval(['Fit_oversamp_reorder(iind,:,:,iRun) = Fit_struct(area).True.s_' str_run(iRun,:) '.k2_plot_oversamp;']);                        
                    end
                else
                    resp_avg_reorder(iind,1:nCond,1) = resp_avg(area,1:nCond);
                    resp_avg_reorder(iind,nCond+1,1) = NaN;
                    TC_avg_reorder(:,iind,1:nCond,1) = TC_avg(:,area,1:nCond);
                    TC_avg_reorder(:,iind,nCond+1,1) = NaN;
%                     CM_mat_reorder(iind,:,1) = Fit_struct(area).True.s_allrun.CM_data';
                    Fit_reorder(iind,:,1) = Fit_struct(area).True.s_allrun.x';       
                    Fit_oversamp_reorder(iind,:,:,1) = Fit_struct(area).True.s_allrun.k2_plot_oversamp;
                end
            end
        end
    end

    resp_avg_all(:,:,:,iexp) = resp_avg_reorder;
    TC_avg_all(:,:,:,:,iexp) = TC_avg_reorder;
%     CM_mat_all(:,:,:,iexp) = CM_mat_reorder;
    H_ttest_all(:,:,iexp) = H_ttest_reorder;
    Fit_all(:,:,:,iexp) = Fit_reorder;
    Fit_oversamp_all(:,:,:,:,iexp) = Fit_oversamp_reorder;
%     F_all(:,:,iexp) = F_reorder;
end

%make oversamp uvars
x = zeros(nCond,2);
SF_vec0 = flipud(uvar(:,2)); %flipped to have low to high SF in square  %flipud
TF_vec0 = uvar(:,1);
[tftf,sfsf] = meshgrid(TF_vec0,SF_vec0); 
grid2.sfsf = sfsf;
grid2.tftf = tftf;
x(:,1) = log2(grid2.sfsf(:));
x(:,2) = log2(grid2.tftf(:));
Nperfreq = sqrt(size(x,1));
xhigh = reshape(x(:,1),Nperfreq,Nperfreq);
xhigh2 = interp2(xhigh,4);
Nperfreq2 = (size(xhigh2,1));
yhigh = reshape(x(:,2),Nperfreq,Nperfreq);
yhigh2 = interp2(yhigh,4);
uvar_high = zeros(Nperfreq2,2);
uvar_high(:,2) = xhigh2(:,1);
uvar_high(:,1) = yhigh2(1,:)';

% averages
Fit_all_invlog = Fit_all;
Speed = zeros(length(all_areas),3,nexp);
for iArea = 1:length(all_areas)
    for iexp = 1:nexp
        run = exp_list.run_mat{iexp};
        if run == 1
            nRunning =3;
        else
            nRunning = 1;
        end
        for iRun = 1:nRunning
            if H_ttest_all(iArea,iRun,iexp)== 1;
                Fit_all_invlog(iArea,4:5,iRun,iexp) = 2.^(Fit_all(iArea,4:5,iRun,iexp));
                Speed(iArea,iRun,iexp) = Fit_all_invlog(iArea,5,iRun,iexp)./Fit_all_invlog(iArea,4,iRun,iexp);
            else
                Fit_all_invlog(iArea,4:5,iRun,iexp) = NaN;
                Speed(iArea,iRun,iexp) = NaN;
            end                
        end    
    end
end


% CM_mat_avg = zeros(2,length(all_areas),3);
% CM_mat_sem = zeros(2,length(all_areas),3);
resp_avg_norm = zeros(size(resp_avg_all));
resp_avg_avg = zeros(nCond+1,length(all_areas),3);
TC_avg_norm = zeros(size(TC_avg_all));
TC_avg_avg = zeros(nON+nOFF,nCond+1,length(all_areas),3);
TC_avg_sem = zeros(nON+nOFF,nCond+1,length(all_areas),3);
Fit_oversamp_norm = zeros(size(Fit_oversamp_all));
Fit_oversamp_avg= zeros(oversamp,oversamp,length(all_areas),3);
Fit_all_avg = zeros(7,length(all_areas),3);
Fit_all_sem = zeros(7,length(all_areas),3);
nexp_per_area = zeros(length(all_areas),3);

for iArea = 1:length(all_areas)
    for iexp = 1:nexp
        resp_avg_norm(iArea,:,:,iexp) = resp_avg_all(iArea,:,:,iexp)./max(max(resp_avg_all(iArea,:,:,iexp),[],2),[],3);
        TC_avg_norm(:,iArea, :,:,iexp) = TC_avg_all(:,iArea,:,:,iexp)./max(max(max(mean(TC_avg_all(post_win(1):post_win(2),:,:,:,iexp),1),[],4),[],3),[],2);
        Fit_oversamp_norm(iArea,:,:,:,iexp)=Fit_oversamp_all(iArea,:,:,:,iexp)./max(max(max(Fit_oversamp_all(iArea,:,:,:,iexp),[],4),[],3),[],2);
    end
end

%order of Fit_all_avg: 
%1st dim: dF; sigma SF; sigma TF; peak SF; peak TF; Xi; peak speed 
for iArea = 1:length(all_areas)
    for iRun = 1:3;
        ind = find(H_ttest_all(iArea,iRun,:)== 1);
        if ind>0
%             CM_mat_avg(:,iArea,iRun) = squeeze(nanmean(CM_mat_all(iArea,:,iRun,ind),4));
%             CM_mat_sem(:,iArea,iRun) = squeeze(nanstd(CM_mat_all(iArea,:,iRun,ind),[],4)/sqrt(size(ind,1)));
            resp_avg_avg(:,iArea,iRun) = squeeze(nanmean(resp_avg_norm(iArea,:,iRun,ind),4));
            TC_avg_avg(:,:,iArea,iRun) = squeeze(nanmean(TC_avg_norm(:,iArea,:,iRun,ind),5));
            TC_avg_sem(:,:,iArea,iRun) = nanstd(TC_avg_norm(:,iArea,:,iRun,ind),[],5)./sqrt(size(ind,1));
            Fit_all_avg(1:6,iArea,iRun) = squeeze(nanmean(Fit_all_invlog(iArea,:,iRun,ind),4));
            Fit_all_sem(1:6,iArea,iRun) = squeeze(nanstd(Fit_all_invlog(iArea,:,iRun,ind),[],4)/sqrt(size(ind,1)));
            Fit_all_avg(7,iArea,iRun) = nanmean(Speed(iArea,iRun,ind),3);
            Fit_all_sem(7,iArea,iRun) =  squeeze(nanstd(Speed(iArea,iRun,ind),[],3)/sqrt(size(ind,1)));
            Fit_oversamp_avg(:,:,iArea,iRun) = squeeze(nanmean(Fit_oversamp_norm(iArea,:,:,iRun,ind),5));
            nexp_per_area(iArea,iRun) = size(ind,1); 
        end
    end
end

% Log_fit_all_avg= zeros(3,6,length(all_areas),3);
% %order of Log_fit_all_avg: 
% %1st dim: SF; TF; speed
% %2nd dim: mu; sigma; lb_mu; ub_mu; lb_sigma; ub_sigma
% for iArea = 2:length(all_areas)
%     for iRun = 1:3;
%         ind = find(H_ttest_all(iArea,iRun,:)== 1);
%         if ind>0
%             [parmhat_SF parmci_SF] = lognfit(squeeze(Fit_all_invlog(iArea,4,iRun,ind)));
%             [parmhat_TF parmci_TF] = lognfit(squeeze(Fit_all_invlog(iArea,5,iRun,ind)));
%             [parmhat_speed parmci_speed] = lognfit(squeeze(Speed(iArea,iRun,ind)));
%             Log_fit_all_avg(1,:,iArea,iRun) =exp([parmhat_SF reshape(parmci_SF,[1 4])]);
%             Log_fit_all_avg(2,:,iArea,iRun) = exp([parmhat_TF reshape(parmci_TF,[1 4])]);
%             Log_fit_all_avg(3,:,iArea,iRun) = exp([parmhat_speed reshape(parmci_speed,[1 4])]);
%         end
%     end
% end

fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_fit_all_summary.mat']);
save(fn_out, 'all_areas', 'resp_avg_all', 'TC_avg_all', 'Fit_all', 'Fit_oversamp_all','H_ttest_all', 'mice', 'nexp_per_area');


[max_resp cond] = max(resp_avg_all,[],2);
max_dF = max(max(max(max_resp,[],4),[],3),[],1);
sf_vars = [reshape(grid2.sfsf',1,nCond) 0];
tf_vars = [reshape(grid2.tftf',1,nCond) 0];

sf_max = zeros(length(all_areas),3,nexp);
tf_max = zeros(length(all_areas),3,nexp);
speed_max = zeros(length(all_areas),3,nexp);
df_max = zeros(length(all_areas),3,nexp);

for iexp = 1:nexp
    if exp_list.run_mat{iexp}==1
        for iArea = 1:length(all_areas)
            for iRun = 2:3
                sf_max(iArea,iRun,iexp) = sf_vars(cond(iArea,:,iRun,iexp));
                tf_max(iArea,iRun,iexp) = tf_vars(cond(iArea,:,iRun,iexp));
                speed_max(iArea,iRun,iexp) = tf_max(iArea,iRun,iexp)/sf_max(iArea,iRun,iexp);
                df_max(iArea,iRun,iexp) = max_resp(iArea,:,iRun,iexp);
            end
        end
    else
        sf_max(:,:,iexp) = NaN;
        tf_max(:,:,iexp) = NaN;
        speed_max(:,:,iexp) = NaN;
        df_max(:,:,iexp) = NaN;
    end
end

sf_max_ratio = sf_max(:,2,:)./sf_max(:,3,:);
tf_max_ratio = tf_max(:,2,:)./tf_max(:,3,:);
speed_max_ratio = speed_max(:,2,:)./speed_max(:,3,:);
df_max_ratio = df_max(:,2,:)./df_max(:,3,:);

sf_max_avg = zeros(length(all_areas),1);
sf_max_sem = zeros(length(all_areas),1);
tf_max_avg = zeros(length(all_areas),1);
tf_max_sem = zeros(length(all_areas),1);
speed_max_avg = zeros(length(all_areas),1);
speed_max_sem = zeros(length(all_areas),1);
df_max_avg = zeros(length(all_areas),1);
df_max_sem = zeros(length(all_areas),1);

for iArea = 1:length(all_areas)
    ind = find(cond(iArea,:,3,:)<(nCond+1) & max_resp(iArea,:,3,:)>0);
    if size(ind)>0
    sf_max_avg(iArea,1) = squeeze(mean(sf_max_ratio(iArea,:,ind),3));
    sf_max_sem(iArea,1) = squeeze(std(sf_max_ratio(iArea,:,ind),[],3)./sqrt(size(sf_max_ratio(iArea,:,ind),3)));
    tf_max_avg(iArea,1) = squeeze(mean(tf_max_ratio(iArea,:,ind),3));
    tf_max_sem(iArea,1) = squeeze(std(tf_max_ratio(iArea,:,ind),[],3)./sqrt(size(tf_max_ratio(iArea,:,ind),3)));
    speed_max_avg(iArea,1) = squeeze(mean(speed_max_ratio(iArea,:,ind),3));
    speed_max_sem(iArea,1) = squeeze(std(speed_max_ratio(iArea,:,ind),[],3)./sqrt(size(speed_max_ratio(iArea,:,ind),3)));
    df_max_avg(iArea,1) = squeeze(mean(df_max_ratio(iArea,:,ind),3));
    df_max_sem(iArea,1) = squeeze(std(df_max_ratio(iArea,:,ind),[],3)./sqrt(size(df_max_ratio(iArea,:,ind),3)));
    end
end

%figures

%averaged SFxTF matrix
for iRun = 1
    figure
    start = 1;
    MAX = max(max(resp_avg_avg(:,:,iRun),[],2),[],1);
    MIN = 0;
    tit = [matrix ' ' inj ' ' image ' Resp matrix ' str_run(iRun,:)];
    suptitle(tit)
    for iArea = 1:length(all_areas)
        if nexp_per_area(iArea,iRun)>0
            subplot(ceil(sqrt(length(all_areas))), ceil(sqrt(length(all_areas))),start);      
            avg = reshape(resp_avg_avg(1:nCond,iArea,iRun),[sqrt(nCond) sqrt(nCond)]);
            imagesc(avg',[MIN MAX]); colormap('gray'); axis image;colorbar;
            n_exp= nexp_per_area(iArea,iRun);
            area = all_areas(iArea,:);
            title_txt = [area '  (' num2str(n_exp) ')'];
            title(title_txt);
            start = start+1;
        end
    end
    fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' image '_matrix.pdf']);
    print(gcf, '-dpdf', fn_out);
end

%averaged oversampled fits
MAX = 1;
MIN = 0;
for iRun = 1
    figure;
    tit = [matrix ' ' inj ' ' image ' Oversampled Fits ' str_run(iRun,:)];
    suptitle(tit)
    start = 1;
    for iArea = 1:length(all_areas)
        if nexp_per_area(iArea,iRun)>0
            h = subplot(ceil(sqrt(length(all_areas))), ceil(sqrt(length(all_areas))),start);      
            imagesc(Fit_oversamp_avg(:,:,iArea,iRun),[MIN MAX]); colormap('gray'); axis image;colorbar;
            n_exp= nexp_per_area(iArea,iRun);
            area = all_areas(iArea,:);
            title_txt = [area '  (' num2str(n_exp) ')'];
            title(title_txt);
            hold on
            for iexp = 1:nexp
                SF_fit = Fit_all(iArea,4,iRun,iexp);
                TF_fit = Fit_all(iArea,5,iRun,iexp);
                set(h,'XTick',[1:size(uvar_high,1)])
                set(h,'YTick',[1:size(uvar_high,1)])
                set(h,'YTickLabel',flipud(uvar_high(:,2)));set(h,'XTickLabel',(uvar_high(:,1)));
                xtxt = interp1(str2num(get(h,'XTickLabel')),get(h,'XTick')', TF_fit);
                ytxt = interp1(str2num(get(h,'YTickLabel')),flipud(get(h,'YTick')'),SF_fit);
                plot(xtxt,ytxt,'k.');
                axis off
                hold on
            end
        end
        start = start+1;
    end
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' image '_oversamp_fits.pdf']);
print(gcf, '-dpdf', fn_out);
end

%TCs by area
iRun = 1;
for iArea= 1:length(all_areas)
    figure; 
    tit = [matrix ' ' inj ' ' image ' Timecourse ' all_areas(iArea,:)];
    suptitle(tit)
    my = max(max(max(TC_avg_avg(:,:,find(nexp_per_area(:,iRun)>0),iRun),[],3),[],2),[],1);
    for iCond = 1:nCond;
        subplot(sqrt(nCond)+1,sqrt(nCond),iCond);
        if nexp_per_area(iArea,iRun)>0
        errorbar(1:(nOFF+nON), TC_avg_avg(:,iCond,iArea,iRun), TC_avg_sem(:,iCond,iArea,iRun), col(iArea));
        ylim([-.5 (my + 0.1*my)]);
        xlim([0 nON+nOFF]);
        hold on;
        end          
    end
    fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' image '_' all_areas(iArea,:) '_TCs.pdf']);
    print(gcf, '-dpdf', fn_out);
end

%summary fig
for iRun = 1;
figure;
tit = [matrix ' ' inj ' ' image ' Fit summary ' str_run(iRun,:)];
suptitle(tit)
h = subplot(4,2,1);
errorbar(1:length(all_areas), Fit_all_avg(4,:,iRun),  Fit_all_sem(4,:,iRun), 'k.')
hold on
bar(Fit_all_avg(4,:,iRun),'k');
ylim([0 0.35])
title('Peak SF')
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
h=subplot(4,2,2);
errorbar(1:length(all_areas),  Fit_all_avg(5,:,iRun), Fit_all_sem(5,:,iRun), 'k.')
hold on
bar(Fit_all_avg(5,:,iRun),'k');
ylim([0 20])
title('Peak TF')
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
h= subplot(4,2,3);
errorbar(1:length(all_areas),  Fit_all_avg(7,:,iRun), Fit_all_sem(7,:,iRun), 'k.')
hold on
bar(Fit_all_avg(7,:,iRun),'k');
ylim([0 500])
title('Peak Speed')
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
h = subplot(4,2,4);
errorbar(1:length(all_areas), Fit_all_avg(1,:,iRun), Fit_all_sem(1,:,iRun), 'k.')
hold on
bar(Fit_all_avg(1,:,iRun),'k');
title('Peak dF')
ylim([0 1])
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
h = subplot(4,2,5);
errorbar(1:length(all_areas), Fit_all_avg(2,:,iRun), Fit_all_sem(2,:,iRun), 'k.')
hold on
bar(Fit_all_avg(2,:,iRun),'k');
title('SF sigma')
ylim([0 6])
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
h = subplot(4,2,6);
errorbar(1:length(all_areas), Fit_all_avg(3,:,iRun), Fit_all_sem(3,:,iRun), 'k.')
hold on
bar(Fit_all_avg(3,:,iRun),'k');
title('TF sigma')
ylim([0 6])
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
h = subplot(4,2,7);
errorbar(1:length(all_areas), Fit_all_avg(6,:,iRun), Fit_all_sem(6,:,iRun), 'k.')
hold on
bar(Fit_all_avg(6,:,iRun),'k');
title('Xi')
ylim([-2 2])
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' image '_fit_summary.pdf']);
print(gcf, '-dpdf', fn_out);
end


%run/no run TCs
for iArea = 1:length(all_areas);
    my = max(max(max(TC_avg_avg(:,:,iArea,:),[],4),[],2),[],1);
    if nexp_per_area(iArea,2)>0
        figure;
        tit = [matrix ' ' inj ' ' image 'Timecourse ' all_areas(iArea,:) ' run/norun'];
        suptitle(tit)
        for iCond = 1:nCond;
            subplot(sqrt(nCond),sqrt(nCond),iCond);
            errorbar(1:20, TC_avg_avg(:,iCond,iArea,2), TC_avg_sem(:,iCond,iArea,2), col(iArea));
            hold on;
            errorbar(1:20, TC_avg_avg(:,iCond,iArea,3), TC_avg_sem(:,iCond,iArea,3), 'Color', [0.8 0.8 0.8]);
            ylim([-.5 (my + 0.1*my)]);
            xlim([0 25]);
        end
        fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_POST' num2str(post_win) mouse_list '_' all_areas(iArea,:) '_run_norun_TCs.ps']);
        print(gcf, '-depsc', fn_out);
    end
end

figure;
subplot(2,2,1)
for iexp = 1:nexp;
    if H_ttest_all(4,1,iexp) == 1;
        if H_ttest_all(2,1,iexp) == 1;
            plot(Speed(4,1,iexp), Speed(2,1,iexp),'*k');
            hold on
        end
    end
end
xlabel('LM');
ylabel('PM');
subplot(2,2,2)
for iexp = 1:nexp;
    if H_ttest_all(5,1,iexp) == 1;
        if H_ttest_all(2,1,iexp) == 1;
            plot(Speed(5,1,iexp), Speed(2,1,iexp),'*k');
            hold on
        end
    end
end
xlabel('AL');
ylabel('PM');
subplot(2,2,3)
for iexp = 1:nexp;
    if H_ttest_all(6,1,iexp) == 1;
        if H_ttest_all(2,1,iexp) == 1;
            plot(Speed(6,1,iexp), Speed(2,1,iexp),'*k');
            hold on
        end
    end
end
xlabel('RL');
ylabel('PM');
subplot(2,2,4)
for iexp = 1:nexp;
    if H_ttest_all(3,1,iexp) == 1;
        if H_ttest_all(2,1,iexp) == 1;
            plot(Speed(3,1,iexp), Speed(2,1,iexp),'*k');
            hold on
        end
    end
end
xlabel('AM');
ylabel('PM');
fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_POST' num2str(post_win) mouse_list '_PM_speed_corr.ps']);
print(gcf, '-depsc', fn_out);

ub_error = squeeze(Log_fit_all_avg(:,4,:,:)-Log_fit_all_avg(:,1,:,:));
for iRun = 1;
figure;
tit = [matrix ' ' inj ' ' image ' Fit summary log estimates ' str_run(iRun,:)];
suptitle(tit)
h = subplot(4,2,1);
errorbar(1:length(all_areas), Log_fit_all_avg(1,1,:,iRun),  ub_error(1,:,iRun), 'k.')
hold on
bar(squeeze(Log_fit_all_avg(1,1,:,iRun)),'k');
ylim([0 0.35])
title('Peak SF')
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
h=subplot(4,2,2);
errorbar(1:length(all_areas),  Log_fit_all_avg(2,1,:,iRun), ub_error(2,:,iRun), 'k.')
hold on
bar(squeeze(Log_fit_all_avg(2,1,:,iRun)),'k');
ylim([0 20])
title('Peak TF')
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
h= subplot(4,2,3);
errorbar(1:length(all_areas),  Log_fit_all_avg(3,1,:,iRun), ub_error(3,:,iRun), 'k.')
hold on
bar(squeeze(Log_fit_all_avg(3,1,:,iRun)),'k');
ylim([0 500])
title('Peak Speed')
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
end

figure;
count = 1;
for iArea1 = 2:7;
    for iArea2 = 2:7
        if iArea1 < iArea2
            subplot(4,4,count)
            for iexp = 1:nexp
                if H_ttest_all(iArea1,1,iexp) + H_ttest_all(iArea2,1,iexp) == 2;
                    semilogy([1:2],Speed([iArea1 iArea2],1,iexp),'-o');
                    hold on
                end
            end
            tit_area1 = all_areas(iArea1,:);
            tit_area2 = all_areas(iArea2,:);
            title([tit_area1 ' vs ' tit_area2])
            xlim([0 3])
            ylim([10^0 10^4])
            count = 1+count;
        end
    end
end

fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_POST' num2str(post_win) mouse_list '_' str_run(iRun,:) '_speed_area_comp.ps']);
print(gcf, '-depsc', fn_out);

figure;
col_exp = ['k' 'r' 'm' 'b' 'c' 'g' 'y' 'k' 'r' 'm' 'b' 'c' 'g' 'y'];
subplot(2,1,1)
for iexp = 1:nexp;
    h= semilogy([1:6], Speed(2:7, 1, iexp), '-ok');
    hold on;
end
xlim([0 7])
subplot(2,1,2)
for iexp = 1:nexp;
    h= semilogy([1:6], Speed(2:7, 1, iexp), ['-o' col_exp(iexp)]);
    hold on;
end
xlim([0 7])

fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_POST' num2str(post_win) mouse_list '_' str_run(iRun,:) '_speed_all_areas_black.ps']);
print(gcf, '-depsc', fn_out);

figure;
x = 0:1:1000;
y = x;
suptitle([matrix ' ' inj ' ' image ' run/norun scatter']);
for iArea = 1:length(all_areas)
    for iexp = 1:nexp
        if H_ttest_all(iArea,2,iexp) == 1;
            if H_ttest_all(iArea,3,iexp) == 1;
                subplot(2,2,1)
                scatter(squeeze(Fit_all_invlog(iArea,4,2,iexp)),  squeeze(Fit_all_invlog(iArea,4,3,iexp)), col(iArea));
                hold on
                ylim([0 0.35])
                xlim([0 0.35])
                title('Peak SF')                
                subplot(2,2,2)
                scatter(squeeze(Fit_all_invlog(iArea,5,2,iexp)),  squeeze(Fit_all_invlog(iArea,5,3,iexp)), col(iArea));                
                hold on
                ylim([0 15])
                xlim([0 15])
                title('Peak TF')  
            end
        end
    end    
end

for iArea = 1:length(all_areas)
    for iexp = 1:nexp
        if H_ttest_all(iArea,2,iexp) == 1;
            if H_ttest_all(iArea,3,iexp) == 1;
                subplot(2,2,3)
                scatter(squeeze(Speed(iArea,2,iexp)),  squeeze(Speed(iArea,3,iexp)), col(iArea));
                hold on
                ylim([0 500])
                xlim([0 500])
                title('Peak Speed')                
                subplot(2,2,4)
                scatter(squeeze(Fit_all_invlog(iArea,1,2,iexp)),  squeeze(Fit_all_invlog(iArea,1,3,iexp)), col(iArea));                
                hold on
                ylim([0 .03])
                xlim([0 .03])
                title('Peak dF')  
            end
        end
    end    
end
for sub= 1:4
    subplot(2,2,sub)
    xlabel('running')
    ylabel('stationary')
    plot(y,x,'-k');
end
fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_POST' num2str(post_win) mouse_list '_run_norun_Fit_scatter.ps']);
print(gcf, '-depsc', fn_out);

x = 0:0.1:8;
y = ones(size(x));
figure;
suptitle([matrix ' ' inj ' ' image ' Run/Norun Max Response'])
h = subplot(2,2,1);
errorbar(1:length(all_areas), sf_max_avg,  sf_max_sem, 'k.')
title('max SF')
hold on
bar(sf_max_avg,'k');
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
h = subplot(2,2,2);
errorbar(1:length(all_areas), tf_max_avg,  tf_max_sem, 'k.')
title('max TF')
hold on
bar(tf_max_avg,'k');
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
h= subplot(2,2,3);
errorbar(1:length(all_areas), speed_max_avg,  speed_max_sem, 'k.')
title('max speed')
hold on
bar(speed_max_avg,'k');
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
h= subplot(2,2,4);
errorbar(1:length(all_areas), df_max_avg,  df_max_sem, 'k.')
title('max dF')
hold on
bar(df_max_avg,'k');
set(h,'XTick', [1:length(all_areas)]);
set(h,'XTickLabel', all_areas);
for sub= 1:4
    subplot(2,2,sub)
    plot(x,y,'-k');
end

fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_POST' num2str(post_win) mouse_list '_run_norun_Max_response.ps']);
print(gcf, '-depsc', fn_out);

% statistics
Fit_all_Nan = Fit_all;
Fit_all_Nan(find(Fit_all ==0)) = NaN;

dF_allareas = squeeze(Fit_all_Nan(:,1,1,:))';
TF_allareas = squeeze(Fit_all_Nan(:,5,1,:))';
SF_allareas = squeeze(Fit_all_Nan(:,4,1,:))';
speed_allareas = squeeze(Speed(:,1,:))';

[p_dF anovatab_dF stats_dF] = anova1(dF_allareas,all_areas);
[p_TF anovatab_TF stats_TF] = anova1(TF_allareas,all_areas);
[p_SF anovatab_SF stats_SF] = anova1(SF_allareas,all_areas);
[p_speed anovatab_speed stats_speed] = anova1(speed_allareas,all_areas);

c_dF= multcompare(stats_dF);
c_SF= multcompare(stats_SF);
c_TF= multcompare(stats_TF);
c_speed= multcompare(stats_speed);

dF_diff = [];
for icomp = 1:size(c_dF,1);
    if c_dF(icomp,3)<0 & c_dF(icomp,5)<0;
        dF_diff = [dF_diff; c_dF(icomp,1:2)];
    elseif c_dF(icomp,3)>0 & c_dF(icomp,5)>0;
        dF_diff = [dF_diff; c_dF(icomp,1:2)];
    end
end
SF_diff = [];
for icomp = 1:size(c_SF,1);
    if c_SF(icomp,3)<0 & c_SF(icomp,5)<0;
        SF_diff = [SF_diff; c_SF(icomp,1:2)];
    elseif c_SF(icomp,3)>0 & c_SF(icomp,5)>0;
        SF_diff = [SF_diff; c_SF(icomp,1:2)];
    end
end
TF_diff = [];
for icomp = 1:size(c_TF,1);
    if c_TF(icomp,3)<0 & c_TF(icomp,5)<0;
        TF_diff = [TF_diff; c_TF(icomp,1:2)];
    elseif c_TF(icomp,3)>0 & c_TF(icomp,5)>0;
        TF_diff = [TF_diff; c_TF(icomp,1:2)];
    end
end
speed_diff = [];
for icomp = 1:size(c_speed,1);
    if c_speed(icomp,3)<0 & c_speed(icomp,5)<0;
        speed_diff = [speed_diff; c_speed(icomp,1:2)];
    elseif c_speed(icomp,3)>0 & c_speed(icomp,5)>0;
        speed_diff = [speed_diff; c_speed(icomp,1:2)];
    end
end

dF_diff_area = [];
for ipair = 1:size(dF_diff,1);
    dF_diff_area = [dF_diff_area; all_areas(dF_diff(ipair,1),:) ' '  all_areas(dF_diff(ipair,2),:)];
end
SF_diff_area = [];
for ipair = 1:size(SF_diff,1);
    SF_diff_area = [SF_diff_area; all_areas(SF_diff(ipair,1),:) ' '  all_areas(SF_diff(ipair,2),:)];
end
TF_diff_area = [];
for ipair = 1:size(TF_diff,1);
    TF_diff_area = [TF_diff_area; all_areas(TF_diff(ipair,1),:) ' '  all_areas(TF_diff(ipair,2),:)];
end
speed_diff_area = [];
for ipair = 1:size(speed_diff,1);
    speed_diff_area = [speed_diff_area; all_areas(speed_diff(ipair,1),:) ' '  all_areas(speed_diff(ipair,2),:)];
end

fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_POST' num2str(post_win) mouse_list '_stats.mat']);
save(fn_out, 'p_dF', 'anovatab_dF', 'stats_dF', 'p_TF', 'anovatab_TF', 'stats_TF', 'p_SF', 'anovatab_SF', 'stats_SF','c_dF', 'c_SF', 'c_TF', 'c_speed', 'dF_diff_area', 'SF_diff_area', 'TF_diff_area', 'speed_diff_area');
