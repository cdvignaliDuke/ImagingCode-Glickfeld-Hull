P = 2;
matrix = 'SF5xTF5';
inj = 'LM';
image = 'AM';

sum_base = 'G:\users\lindsey\analysisLG\experiments';
list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
load(list_fn);
nexp = size(exp_list.mouse_mat,2);


for iexp = 1:nexp;
    figure;
    beg = 1;
    
    mouse = char(exp_list.mouse_mat{iexp});
    date = char(exp_list.date_mat{iexp});
    userun = exp_list.runs_mat{iexp};
    count_protocol = exp_list.prot_mat{iexp};
    run = exp_list.run_mat{iexp};
    blanks = exp_list.blanks_mat{iexp};
    dir = exp_list.dir_mat{iexp};
    
    base = 'G:\users\lindsey\analysisLG\active mice';
    outDir = fullfile(base, mouse, date);

    if dir ==1;
        nCond = 25;
    elseif dir ==2;
        nCond = 50;
    end
    if length(userun)==1
        fn_avg = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun(end)) '_dF_F_Nstim' num2str(nCond+1) '_POST_1_5.tif']);
    else
        fn_avg = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun(end)) '_dec_reg_dF_F_Nstim' num2str(nCond+1) '_POST_1_5_MEAN.tif']);
    end
    sfxtf_avg =  readtiff(fn_avg);
    
    siz = size(sfxtf_avg);
    
    r = zeros(siz(1), siz(2));
    b= 5;
    for iy = b+1:240-b
        fprintf('.');
        for ix = b+1:256-b
                sub = sfxtf_avg(iy-1:iy+1,ix-1:ix+1,:);
                r(iy,ix)=mean(triu2vec((corrcoef(reshape(sub,[3*3,siz(3)])')),1));
        end;
    end;
    
    sfxtf_long = reshape(sfxtf_avg,[siz(1)*siz(2) siz(3)]);
    corr = [.4:.1:.8];
    axon_mask = zeros(siz(1),siz(2),5);
    resp_avg = zeros(26,5);
    for icorr =1:5    
        axon_mask_temp = zeros(size(r));
        resp_avg_temp  = zeros(26,1);
        axon_mask_ind = find(r>corr(icorr));
        axon_mask_temp(axon_mask_ind) = 1;
        axon_mask(:,:,icorr) = reshape(axon_mask_temp,[siz(1) siz(2)]);
        subplot(5,2,beg)
        imagesq(axon_mask_temp);
        if icorr ==1
        title([mouse ' ' date ' ' image])
        end
        start = 1;
        if dir == 2
            for iCond = 1:25
                resp_avg_temp(iCond,:) = mean(mean(sfxtf_long(axon_mask_ind,start:start+1),2),1);
                start = start+2;
            end
            resp_avg_temp(26,:) = mean(sfxtf_long(axon_mask_ind,end),1);
        elseif dir ==1
            resp_avg_temp = mean(sfxtf_long(axon_mask_ind,:),1)';
        end
        resp_avg(:,icorr) = resp_avg_temp;
        resp_avg_sq = reshape(resp_avg_temp(1:25,:), [5 5])';
        subplot(5,2,beg+1)
        imagesq(resp_avg_sq);colormap(gray); colorbar
        beg = beg+2;
    end
    
    
    fn_out = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_axon_mask_resp_fig.ps']);
    print(gcf, '-depsc', fn_out);
    
    fn_out = fullfile('\\Zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging',[matrix '_' num2str(P) 'P_' image], [date '_' mouse '_run' num2str(userun) '_axon_mask_resp_fig.ps']);
    print(gcf, '-depsc', fn_out);
    
    fn_out = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_axon_mask.mat']);
    save(fn_out, 'r', 'axon_mask');

    fn_out = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp_avg.mat']);
    save(fn_out, 'resp_avg');
end
