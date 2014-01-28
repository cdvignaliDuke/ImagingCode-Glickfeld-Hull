areas = ['RL'; 'AM'];
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
nON = 12;
nOFF = 12;
nPlanes = 1;
begin = 1;
TFSFetc = [1:2];
pre_win = [7 12];
post_win = [13 24];
for iArea = 1
    image = areas(iArea,:);
    
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    FOV_dF_all = zeros(26,nexp);
    SNR_all = zeros(1,nexp);
    ttest_mask_all = zeros(240,256,nexp);
    for iexp = 1:nexp;
    
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_protocol = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};

        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
        load(fn_out);
        
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_sorted.tif']);
        stack_sorted = readtiff(fn_out);
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
        load(fn_out);
        
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
        stack_dF = readtiff(fn_out);

        siz = size(stack_sorted);
        
        nReps = sum(stim_reps(1,:));
        
        roi_stim = zeros(siz(1), siz(2), nOFF+nON,nReps);
        start = 1;
        rep = 1;
        for iRep = 1:nReps;
            roi_stim(:,:,:, rep) = stack_sorted(:,:,start:start-1+nON+nOFF);
            start = start+nON+nOFF;
            rep = rep+1;
        end
        stim_off = squeeze(mean(roi_stim(:,:,pre_win(1):pre_win(2),:),3));
        stim_on = squeeze(mean(roi_stim(:,:,post_win(1):post_win(2),:),3));
        
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_allstim.mat']);
        save(fn_out, 'stim_on', 'stim_off');

        ttest_mask_long = reshape(ttest_mask, [siz(1)*siz(2), 1]);
        ind = find(ttest_mask_long == 1);
        ttest_mask_all(:,:,iexp) = ttest_mask;

        stack_sorted_long = reshape(stack_sorted, [siz(1)*siz(2) siz(3)]);
        clear stack_sorted
        stack_avg = mean(stack_sorted_long(ind,:),1);
        clear stack_sorted_long;

        if dirs ==1
            nCond = 25;
        elseif dirs ==2
            nCond = 50;
        end

        siz = size(stack_dF);
        stack_dF_long = reshape(stack_dF, [siz(1)*siz(2) siz(3)]);
        stack_dF_avg = mean(stack_dF_long(ind,:),1); 
        [peak_dF peak_cond] = max(stack_dF_avg,[],2);
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

        clear stack_avg
                
        resp_on_long = reshape(stim_on, [siz(1)*siz(2) ncyc]);
        resp_off_long = reshape(stim_off, [siz(1)*siz(2) ncyc]);
        
        clear stim on
        clear stimoff
        
        resp_dF_peak = resp_on_long(ind,start_rep+1:nRep+start_rep)- resp_off_long(ind,start_rep+1:nRep+start_rep);
        resp_dF_avg = mean(mean(resp_dF_peak,1),2);
        SNR = resp_dF_avg./noise;
        
        clear resp_on_long
        clear resp_off_long
        
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_SNR.mat']);
        save(fn_out, 'SNR', 'resp_dF_avg', 'noise');

        stack_dF_avg_temp = zeros(1,26);
        if nCond == 50;
            start = 1;
            for iCond = 1:25
                stack_dF_avg_temp(:,iCond) = mean(stack_dF_avg(:,start:start+1),2);
                start = start+2;
            end
            stack_dF_avg_temp(:,end) = stack_dF_avg(:,end);
            stack_dF_avg = stack_dF_avg_temp;
        end
        FOV_dF_all(:,iexp) = stack_dF_avg';
        SNR_all(:,iexp) = SNR;
    end
    [SNR_sort ind_order] = sort(SNR_all);
    ind = find(SNR_all>5);
    FOV_dF_norm = zeros(size(FOV_dF_all));
    for iexp = 1:nexp
        FOV_dF_norm(:,iexp) = FOV_dF_all(:,iexp)./max(FOV_dF_all(:,iexp),[],1);
    end
    figure;
    start =1;
    for iexp = 1:nexp
        subplot(ceil(sqrt((nexp+1)*2)),ceil(sqrt(nexp*2)),start)
        imagesq(ttest_mask_all(:,:,ind_order(iexp)));
        colormap(gray);
        title([char(exp_list.date_mat{ind_order(iexp)}) ' ' char(exp_list.mouse_mat{ind_order(iexp)})]);
        subplot(ceil(sqrt(nexp*2)),ceil(sqrt(nexp*2)),start+1)
        imagesq(reshape(FOV_dF_norm(1:25,ind_order(iexp)),[5 5])');
        colormap(gray);
        title(num2str(SNR_sort(iexp)));
        start = start+2;
    end
    subplot(ceil(sqrt((nexp+1)*2)),ceil(sqrt(nexp*2)),start)
    imagesq(reshape(mean(FOV_dF_norm(1:25,ind),2),[5 5])');
    colormap(gray);
    title('Norm avg SNR>5');
    fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_resp_avg.pdf']);
    print(gcf, '-dpdf', fn_out);
    fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_resp_avg.mat']);
    save(fn_out, 'SNR_all', 'ttest_mask_all', 'FOV_dF_all');
end
clear all


