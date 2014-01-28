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
    
    if dir ==1
        nCond = 25;
    elseif dir ==2
        nCond = 50;
    end
    
    base = 'G:\users\lindsey\analysisLG\active mice';
    outDir = fullfile(base, mouse, date);
   
    fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
    if exist(fn_stack,'file') == 0
       if length(userun) == 1
            seqfile =[date '_' mouse '_run' num2str(userun) '_Seqposition.mat'];
            load(fullfile(outDir,'analysis',seqfile));
            Big_Seqposition = Seqposition;
       else
            seqfile = [date '_' mouse '_run' num2str(userun) '_Big_Seqposition.mat'];
            load(fullfile(outDir,'analysis',seqfile));
       end
        fn_np = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted_npsub.tif']);
        stack_sorted_np = readtiff(fn_np);
       stack_dF_np = zeros(siz(1), siz(2), nCond+1);
        start = 0;
        for iCond = 1:Cond;
            nRep = length(Big_Seqposition(iCond).ind);
            rep_dF = zeros(siz(1), siz(2),nRep);
            for iRep = 1:nRep
                rep_base = mean(stack_sorted_np(:,:,start+pre_win(1):start+pre_win(2)),3);
                rep_resp = mean(stack_sorted_np(:,:,start+post_win(1):start+post_win(2)),3);
                rep_dF(:,:,iRep) = (rep_resp-rep_base)./rep_base;
                start = ((nOFF+nON)/nPlanes)+start;
            end
            stack_dF_np(:,:,iCond) = mean(rep_dF,3);
        end
        fn_out = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
        writetiff(stack_dF_np, fn_out);
    end
end
end