areas = ['PM'; 'LM'; 'AL'];
for iArea = 1:3
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
image = areas(iArea,:);
nON = 12;
nOFF = 12;
nPlanes = 1;
begin = 1;
TFSFetc = [1:2];
pre_win = [7 12];
post_win = [13 24];

sum_base = 'G:\users\lindsey\analysisLG\experiments';
list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
load(list_fn);
nexp = size(exp_list.mouse_mat,2);

for iexp = 1:nexp;
    
    mouse = char(exp_list.mouse_mat{iexp});
    date = char(exp_list.date_mat{iexp});
    userun = exp_list.runs_mat{iexp};
    count_protocol = exp_list.prot_mat{iexp};
    run = exp_list.run_mat{iexp};
    blanks = exp_list.blanks_mat{iexp};
    dir = exp_list.dir_mat{iexp};
    
    base = 'G:\users\lindsey\analysisLG\active mice';    
    outDir = fullfile(base, mouse,date);
    
    fn_SNR = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_SNR.mat']);
    if exist(fn_SNR, 'file')== 0
    
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_sorted_npsub.tif']);
    stack_sorted = readtiff(fn_out);
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
    load(fn_out);
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
    load(fn_out);
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
    load(fn_out);
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
    stack_dF = readtiff(fn_out);
    
    siz = size(stack_sorted);
    
    ttest_mask_long = reshape(ttest_mask, [siz(1)*siz(2), 1]);
    ind = find(ttest_mask_long == 1);
    
    stack_sorted_long = reshape(stack_sorted, [siz(1)*siz(2) siz(3)]);
    clear stack_sorted
    
    stack_avg = mean(stack_sorted_long(ind,:),1);
    clear stack_sorted_long;
    
    siz = size(stack_dF);
    stack_dF_long = reshape(stack_dF, [siz(1)*siz(2) siz(3)]);
    [peak_dF peak_cond] = max(mean(stack_dF_long(ind,:),1),[],2);
    nRep = stim_reps(1,peak_cond);
    start_rep = sum(stim_reps(1,1:(peak_cond-1)));
    
    ncyc = size(stack_avg,2)/(nOFF+nON);
    stdev= zeros(nRep, 1);
    start = 1;
    for iRep = start_rep+1:nRep+start_rep;
        stdev(start,:) = std(stack_avg(:,pre_win(1)+((iRep-1)*(nON+nOFF)):(pre_win(2)-1+((iRep-1)*(nON+nOFF)))));
        start = 1+start;
    end
    noise = mean(stdev,1);
   
    resp_on_long = reshape(resp_on, [siz(1)*siz(2) ncyc]);
    resp_off_long = reshape(resp_off, [siz(1)*siz(2) ncyc]);
    resp_dF_peak = resp_on_long(ind,start_rep+1:nRep+start_rep)- resp_off_long(ind,start_rep+1:nRep+start_rep);
    resp_dF_avg = mean(mean(resp_dF_peak,1),2);
    SNR = resp_dF_avg./noise;
    
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_SNR.mat']);
    save(fn_out, 'SNR', 'resp_dF_avg', 'noise');
    end
end
end


