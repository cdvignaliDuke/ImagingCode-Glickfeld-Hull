areas = ['PM'; 'AL'; 'V1'];
P = 2;
matrix = 'SF5xTF5';
inj = 'LM';
nON = 12;
nOFF = 12;
nPlanes = 1;
begin = 1;
TFSFetc = [1:2];
pre_win = [7 12];
post_win = [13 24];
for iArea = 1:3
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
        
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_SNR.mat']);
        load(fn_out);
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
        load(fn_out); 
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
        stack_dF = readtiff(fn_out);

        siz = size(stack_dF);

        ttest_mask_long = reshape(ttest_mask, [siz(1)*siz(2), 1]);
        ind = find(ttest_mask_long == 1);
        
        stack_dF_long = reshape(stack_dF, [siz(1)*siz(2) siz(3)]);
        stack_dF_avg = mean(stack_dF_long(ind,:),1); 
        stack_dF_avg_temp = zeros(1,26);
        
        if dirs ==1
            nCond = 25;
        elseif dirs ==2
            nCond = 50;
        end
        
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
        ttest_mask_all(:,:,iexp) = ttest_mask;
    end
    [SNR_sort ind_order] = sort(SNR_all);
    ind = find(SNR_all>5);
    FOV_dF_norm = zeros(size(FOV_dF_all));
    for iexp = 1:nexp
        FOV_dF_norm(:,iexp) = FOV_dF_all(:,iexp)./max(FOV_dF_all(:,iexp),[],1);
    end
    FOV_dF = mean(FOV_dF_norm,2);
    figure;
    start =1;
    for iexp = 1:nexp
        subplot(4,4,start)
        imagesq(ttest_mask_all(:,:,ind_order(iexp)));
        colormap(gray);
        title([char(exp_list.date_mat{ind_order(iexp)}) ' ' char(exp_list.mouse_mat{ind_order(iexp)})]);
        subplot(4,4,start+1)
        imagesq(reshape(FOV_dF_norm(1:25,ind_order(iexp)),[5 5])');
        colormap(gray);
        title(num2str(SNR_sort(iexp)));
        start = start+2;
    end
    subplot(4,4,start)
    imagesq(reshape(FOV_dF(1:25,:),[5 5])');
    title('Avg Norm SNR>5');
    suptitle([matrix ' ' num2str(P) 'P ' inj ' ' image]);
    
    anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120120';
    fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_' image '_resp_avg.pdf']);
    print(gcf, '-dpdf', fn_out);
    
    fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_resp_avg.mat']);
    save(fn_out, 'SNR_all', 'ttest_mask_all', 'FOV_dF_all');
end
