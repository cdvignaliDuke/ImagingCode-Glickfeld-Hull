areas = ['PM'; 'LM'; 'AL'];
for iArea = 3
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

for iexp = [13 14];
    
    mouse = char(exp_list.mouse_mat{iexp});
    date = char(exp_list.date_mat{iexp});
    userun = exp_list.runs_mat{iexp};
    count_protocol = exp_list.prot_mat{iexp};
    run = exp_list.run_mat{iexp};
    blanks = exp_list.blanks_mat{iexp};
    dir = exp_list.dir_mat{iexp};
    
    if dir ==1
        nCond = 25;
    elseif dir ==2
        nCond = 50;
    end
    
    base = 'G:\users\lindsey\analysisLG\active mice';
    outDir = fullfile(base, mouse, date);
   
    fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted.tif']);
    if exist(fn_stack,'file') == 0
        stack_sort
    elseif exist(fn_stack,'file') == 2
        stack_sorted = readtiff(fn_stack);
    end

    siz = size(stack_sorted);
    stack_sorted_np = zeros(size(stack_sorted));
    avg = mean(mean(mean(stack_sorted,3),2),1);
    for iframe = 1:siz(3)
        np_av = mean(mean(stack_sorted(:,:,iframe),2),1);
        stack_sorted_np(:,:,iframe) = stack_sorted(:,:,iframe)-np_av+avg;
    end
    
    fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted_npsub.tif']);
    writetiff(stack_sorted_np, fn_stack);
    
    clear stack_sorted
    
    fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
    load(fn);

    nReps = sum(stim_reps(1,:));
    if dir == 2
        add = 1;
        stim_reps_dir = zeros(1,26);
        for iCond = 1:nCond/2
        stim_reps_dir(1,iCond) = sum(stim_reps(1,add:add+1));
        add = add+2;
        end
        stim_reps_dir(1,end) = stim_reps(1,end);
        stim_reps = stim_reps_dir;
    end

    roi_stim = zeros(siz(1), siz(2), nOFF+nON,nReps);
    start = 1;
    rep = 1;
    for iCond = 1:length(stim_reps); 
        nRep = stim_reps(iCond);
        for iRep = 1:nRep;
            roi_stim(:,:,:, rep) = stack_sorted_np(:,:,start:start-1+nON+nOFF);
            start = start+nON+nOFF;
            rep = rep+1;
        end
    end
    resp_off = squeeze(mean(roi_stim(:,:,pre_win(1):pre_win(2),:),3));
    resp_on = squeeze(mean(roi_stim(:,:,post_win(1):post_win(2),:),3));

    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
    save(fn_out, 'resp_on', 'resp_off');

    clear stack_sorted_np
    clear roi_stim
    
    alphaB = .05./(25);
    Info_ttest_mat = zeros(siz(1),siz(2),25);

    b= 5;

    f1 = fspecial('average');
    siz = size(resp_on);
    resp_on_long = reshape(resp_on, [siz(1) siz(2)*siz(3)]);
    resp_off_long = reshape(resp_off, [siz(1) siz(2)*siz(3)]);
    
    clear resp_on
    clear resp_off
    
    resp_on_sm = reshape(filter2(f1,resp_on_long),[siz(1) siz(2) siz(3)]);
    resp_off_sm = reshape(filter2(f1,resp_off_long),[siz(1) siz(2) siz(3)]);

    clear resp_on_long
    clear resp_off_long
    
    for iy = b+1:240-b
        fprintf([num2str(iy) ' '])
        for ix = b+1:256-b
            start = 1;
            p_ttestB = zeros(1,1,25);
            for iCond = 1:(length(stim_reps)-1)
                nRep = stim_reps(1,iCond);
                [h_ttestB1,p_ttestB1] = ttest(resp_off_sm(iy,ix,start:start-1+nRep),resp_on_sm(iy,ix,start:start-1+nRep),alphaB,'left');
                p_ttestB(1,1,iCond) = p_ttestB1;
                start = start+nRep;
            end
        Info_ttest_mat(iy,ix,:) = p_ttestB;
        end
    end
    
    clear resp_on_sm
    clear resp_off_sm
    
    siz = size(Info_ttest_mat);
    Info_ttest_mat_long = reshape(Info_ttest_mat, [siz(1) siz(2)*siz(3)]);
    ttest_smooth = reshape(filter2(f1,Info_ttest_mat_long), [siz(1) siz(2) siz(3)]);
    ttest_mask = min(ttest_smooth,[],3) < alphaB;
    
    ttest_mask(1:5,1:end) = 0;
    ttest_mask(1:end, 1:5) = 0;
    ttest_mask(1:end, 251:end) = 0;
    ttest_mask(235:end,1:end) = 0;
    
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
    save(fn_out, 'ttest_mask', 'Info_ttest_mat');
end
end

