anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120209';
for iArea = 1:3
    P = 2;
    matrix = 'SF5xTF5';
    inj = 'V1';
    image = areas(iArea,:);
    nON = 12;
    nOFF = 12;
    nPlanes = 1;

    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    
    neuropil_dFoverF = zeros(nexp,25);
    figure;
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
        fn_neuropil = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_neuropil.mat']);
        load(fn_neuropil)
        
        subplot(4,4,iexp)
        imagesq(npmask_all);
        title(mouse);
        neuropil_dFoverF(iexp, :) = Im_mat_nponly;
    end
    suptitle([areas(iArea,:) 'neuropil masks']);
    fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_' areas(iArea,:) '_neuropil_masks.pdf' ]);
         print(gcf, '-dpdf', fn_out);
    figure;
    for iexp = 1:nexp
        subplot(4,4,iexp)
        imagesq(reshape(neuropil_dFoverF(iexp,:),5,5)');
        caxis([0 .5])
        max_dF = max(neuropil_dFoverF(iexp,:),[],2);
        title(num2str(max_dF));
    end
    colormap(gray)
    suptitle([areas(iArea,:) 'neuropil response']);
    fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_' areas(iArea,:) '_neuropil_response.pdf' ]);
         print(gcf, '-dpdf', fn_out);
end

for iArea = 1:3;
    P = 2;
    matrix = 'SF5xTF5';
    inj = 'V1';
    image = areas(iArea,:);
    nON = 12;
    nOFF = 12;
    nPlanes = 1;
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    figure;
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
        fn_stack = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
        stack_dF = readtiff(fn_stack);
        fn_neuropil = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_neuropil.mat']);
        load(fn_neuropil);

        stack_max = max(stack_dF,[],3);
        stack_max_norm = stack_max./max(max(stack_max,[],1),[],2);
        overlay = stack_max_norm+npmask_all;
        residual_np = overlay;
        residual_np(find(overlay<=1)) = 1;
        
        subplot(4,4,iexp)
        imagesq(residual_np);
        caxis([1 2])
    end
    suptitle([areas(iArea,:) 'neuropil response']);
end
