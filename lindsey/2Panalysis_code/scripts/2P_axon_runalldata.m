P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
image = 'PM';

sum_base = 'G:\users\lindsey\analysisLG\experiments';
list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
load(list_fn);
nexp = size(exp_list.mouse_mat,2);
beg = 1;
figure;
for iexp = 1:nexp
    base = 'G:\users\lindsey\analysisLG\active mice';
    outDir = fullfile(base, mouse, date);

    if dir ==1;
        nCond = 25;
    elseif dir ==2;
        nCond = 50;
    end

    fn_avg = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun(end)) '_dec_reg_dF_F_Nstim' num2str(nCond+1) '_POST_1_5_MEAN.tif']);
    sfxtf_avg =  readtiff(fn_avg);

    siz = size(sfxtf_avg);
    b= 5;
    r = zeros(siz(1),siz(2));
    for iy = b+1:240-b
        fprintf('.');
        for ix = b+1:256-b
            sub = sfxtf_avg(iy-1:iy+1,ix-1:ix+1,:);
            sub_avg = mean(mean(sub,2),1);
            r(iy,ix)=mean(triu2vec((corrcoef(reshape(sub,[3*3,siz(3)])')),1));
        end;
    end;
    axon_mask = zeros(size(r));
    
    axon_mask_ind = find(r>.7);
    axon_mask(axon_mask_ind) = 1;
    subplot(ceil(sqrt(nexp))*2, ceil(sqrt(nexp))*2,beg)
    imagesq(axon_mask);
    
    axon_avg = mean(sfxtf(axon_mask_ind,:),1)';

    resp_avg = zeros(26,1);
    if dir ==2
        start = 1;
        for iCond = 1:25
            resp_avg(iCond,:) = mean(axon_avg(start:start+1,:),1);
            start = start+2;
        end
        resp_avg(26,:) = axon_avg(end,:);
    elseif dir ==1
        resp_avg = axon_avg;
    end
    
    
    resp_avg_sq = reshape(resp_avg(1:25,:), [5 5])';
    subplot(ceil(sqrt(nexp))*2, ceil(sqrt(nexp))*2,beg+1)
    imagesq(resp_avg_sq)
    beg = beg+2;
    fn_out = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp_avg.m']);
    save(fn_out, 'resp_avg');
end







