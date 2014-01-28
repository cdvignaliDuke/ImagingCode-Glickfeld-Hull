P = 2;
matrix = 'SF5xTF5';
image = 'AL';
inj = 'V1';
nCond = 25;
nON=12;
nOFF=12;

sum_base = 'G:\users\lindsey\analysisLG\experiments';
base = 'G:\users\lindsey\analysisLG\active mice';

list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
load(list_fn);
nexp = size(exp_list.mouse_mat,2);
mouse_list = [];
mice = exp_list.mouse_mat;

resp_all = zeros(27,nexp);
ttest_mask_all = zeros(240, 256, nexp);

for iexp = 1:nexp
    mouse = char(exp_list.mouse_mat{iexp});
    date = char(exp_list.date_mat{iexp});
    userun = exp_list.runs_mat{iexp};
    count_prot = exp_list.prot_mat{iexp};
    run = exp_list.run_mat{iexp};
    blanks = exp_list.blanks_mat{iexp};
    dir = exp_list.dir_mat{iexp};
    
    base = 'G:\users\lindsey\analysisLG\active mice';    
    outDir = fullfile(base, mouse,date);

    fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
    stack = readtiff(fn_stack);
    siz = size(stack);
    stack_dF = zeros(siz(1), siz(2), 26);
    
    if dir == 2
        start = 1;
        for iCond = 1:25
            stack_dF(:,:,iCond) = mean(stack(:,:,start:start+1),3);
            start= start+1;
        end
        stack_dF(:,:,end) = stack(:,:,end);
    else
        stack_dF = stack;
    end

    fn_ttest = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
    load(fn_ttest);
    
    ttest_mask_all(:,:,iexp) = ttest_mask;
    ind = find(ttest_mask == 1);
    
    stack_dF_long = reshape(stack_dF, [siz(1)*siz(2) 26]);
    resp_avg = mean(stack_dF_long(ind,:),1);
    
    resp_all(1:26,iexp) = resp_avg';
    
    fn_SNR = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_SNR.mat']);
    load(fn_SNR)
    resp_all(27,iexp) = SNR;
end

figure;
start = 1;
[y,ind] = sort(resp_all(27,:),2);
for iexp = ind;
    resp_avg_sq = reshape(resp_all(1:25,iexp), [5 5])';
    mouse = char(exp_list.mouse_mat{iexp});
    date = char(exp_list.date_mat{iexp});
    SNR = num2str(resp_all(27,iexp));
    max_dF= num2str(max(resp_all(1:25,iexp),[],1));
    subplot(6,6,start);
    imagesq(ttest_mask_all(:,:,iexp));
    title([mouse ' ' date])
    subplot(6,6,start+1);    
    imagesq(resp_avg_sq);
    colormap(gray);
    title(['SNR=' SNR(:,1:4) ' dF=' max_dF(:,1:3)]);
    start = start+2;
end
    
ind = find(resp_all(27,:)>5);
resp_all_norm = zeros(25,nexp);
for iexp = 1:nexp
    resp_all_norm(:,iexp) = resp_all(1:25,iexp)./max(resp_all(1:25,iexp),[],1);
end
resp_norm_avg = mean(resp_all_norm(:,ind),2);
resp_norm_sq = reshape(resp_norm_avg, [5 5])';
subplot(6,6,start)
imagesq(resp_norm_sq);
title('Norm Avg SNR>5');

fn_out = '\\Zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\120112\all_ttestmask_SNR_avg_AL.pdf';
        print(gcf, '-dpdf', fn_out);


        