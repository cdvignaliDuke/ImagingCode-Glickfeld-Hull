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
    
    fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
    stack_dF = readtiff(fn_stack);

    fn_ttest = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
    load(fn_ttest);

    stack_max = max(stack_dF,[],3);
    f1 = fspecial('average');
    stack_max_sm = filter2(f1, stack_max);
    stack_max_ttest = ttest_mask.*stack_max;
    
    siz = size(stack_max);
    
    local_max = zeros(siz);
    for iy = 1:siz(1);
        for ix = 1:siz(2);
            if stack_max_ttest(iy,ix)>0
                sub = stack_max_ttest(iy-3:iy+3,ix-3:ix+3);
                if max(max(sub,[],1),[],2)==stack_max_ttest(iy,ix)
                    local_max(iy,ix) = 1;
                end
            end
        end
    end

    n_pix = sum(sum(local_max));

    fn_local = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
    save(fn_local, 'local_max', 'n_pix');

    fn_resp = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp.mat']);
    load(fn_resp);
    
    [i, j] = find(local_max ==1);
    
    resp_dF = zeros(n_pix, size(resp_on,3));
    for ipix = 1:n_pix
        sub_y = [i(ipix)-1:i(ipix)+1]; 
        sub_x = [j(ipix)-1:j(ipix)+1];
        roi_on = squeeze(mean(mean(resp_on(sub_y,sub_x,:),2),1));
        roi_off = squeeze(mean(mean(resp_off(sub_y,sub_x,:),2),1));
        resp_dF(ipix,:) = ((roi_on-roi_off)./roi_off)';
    end

    fn_resp = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp.mat']);
    save(fn_resp, 'resp_on', 'resp_off', 'resp_dF')
    
    
end
end