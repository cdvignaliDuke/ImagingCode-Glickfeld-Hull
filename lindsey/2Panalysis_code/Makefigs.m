dir = 'G:\users\lindsey\analysisLG\active mice\LG25\110220\analysis\';
run = '110220_LG25_PM_run3';
analysis = '_dec_reg_tc_Nstim33';

fn_out = fullfile(dir,sprintf('%s_Seqposition.mat',run));
save(fn_out, 'Seqposition');

%make dF/F of resorted data

stack = readtiff([dir run analysis '.tif']);
for iCond = 1:Nstimtypes0000;
    baseline(1+((iCond-1)*round((Noff1)/2)):iCond*round((Noff1)/2)) = Ntot1*(iCond-1)+round((Noff1)/2):Ntot1*(iCond-1)+Noff1;
end
base = mean(stack(:,:,baseline),3);
stack_d = double(stack);
stack_dfoverf = zeros(size(stack));
avg= mean(stack(:,:,:),3);
for time = 1:size(stack,3);
    stack_dfoverf(:,:,time) = (stack_d(:,:,time)-avg)./avg;
end
fn_out = fullfile(dir,sprintf('%s%s_dFoverF.tif',run,analysis));
writetiff(stack_dfoverf,fn_out);

%make average of each stim
[a b c] = size(stack);
stack_collapse = zeros(a,b,Nstimtypes0000);
for iCond = 1:Nstimtypes0000;
    stack_collapse(:,:,iCond) = mean(stack_dfoverf(:,:,((iCond-1)*Ntot1)+(Noff1+1):((iCond-1)*Ntot1)+(Noff1+1)+(round((Noff1)/2))),3);
end
fn_out2 = fullfile(dir,sprintf('%s%s_dFoverF_avg.tif',run,analysis));
writetiff(stack_collapse,fn_out2);

%make plots
stack_collapse_trim = stack_collapse(5:a-5,5:b-5, :);
SF = .16;
TF = 2;
SF_tot = 8;
TF_tot = 8;

%TF plot
SF_ind = find(SFvec==SF);
TF_all = unique(TFvec);
figure;
for iTF = 2:TF_tot+1
    SFTF_comb_ind = SF_ind(find(TFvec(SF_ind) == TF_all(iTF)));
    subplot(2, 4,iTF-1);
    imagesq(stack_collapse_trim(:,:,SFTF_comb_ind(1)));
    colormap gray;
    title([num2str(TF_all(iTF)) ' Hz']);
    hold on;
end

%SF plot
figure;
TF_ind = find(TFvec==TF);
SF_all = unique(SFvec);
for iSF = 1:SF_tot
    TFSF_comb_ind = TF_ind(find(SFvec(TF_ind) == SF_all(iSF)));
    subplot(2, 4,iSF);
    imagesq(stack_collapse_trim(:,:,TFSF_comb_ind(1)));
    colormap gray;
    title([num2str(SF_all(iSF)) ' cyc/deg']);
    hold on;
end


