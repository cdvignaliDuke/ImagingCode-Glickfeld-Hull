
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
image = 'AL';

if matrix == 'SF5xTF5'
    SF_vec0 = [.32; .16; .08; .04; .02];
    TF_vec0 = [1 2 4 8 15];
end

[tftf,sfsf] = meshgrid(TF_vec0,SF_vec0); 
grid2.sfsf = sfsf;
grid2.tftf = tftf;
x(:,1) = log2(grid2.sfsf(:));
x(:,2) = log2(grid2.tftf(:));

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

    stack = readtiff(fn_avg);
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
    fn_mask = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_axon_mask.mat']);
    load(fn_mask)
    stack_dF_long = reshape(stack_dF, [siz(1)*siz(2) 26]);
    r_long = reshape(r,[siz(1)*siz(2),1]);
    r_sorted= sort(r_long);
    CM_data = zeros(3,2);
    percent = [1 10 100];
    resp_avg = zeros(3,25);
    
    for iper = 1:3;
        if percent(iper) < 100
            rank = (siz(1)*siz(2))-ceil(((percent(iper)/100)*(siz(1)*siz(2))-1));
            thresh = r_sorted(rank);
            ind = find(r>thresh);
            resp_avg_temp =mean(stack_dF_long(ind,1:25),1);
        else
            resp_avg_temp = mean(stack_dF_long(:,1:25),1);
        end
        CM = zeros(1,2);
        CM(1,1) = sum(x(:,1).*resp_avg_temp')/sum(resp_avg_temp');
        CM(1,2) = sum(x(:,2).*resp_avg_temp')/sum(resp_avg_temp');
        CM_data(iper,:) = 2.^CM;
        resp_avg(iper,:) = resp_avg_temp;
    end
    fn_CM = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_CM_data.mat']);
    save(fn_CM, 'CM_data');
    fn_resp = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp_avg.mat']);
    save(fn_resp, 'resp_avg');
end