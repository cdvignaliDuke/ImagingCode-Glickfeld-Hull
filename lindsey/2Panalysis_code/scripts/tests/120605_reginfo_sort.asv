reg_base = '\\zmey\storlab\users\Lindsey\DataShare\Layer5\All Reg_info';

date = '120524';
mouse = 'CM63';
userun = [1:4];
nON = 5;
nOFF = 5;
pre_win = [3 5];
post_win = [6 10];
nCond = 25;
P = 2;
nPlanes = 1;
Nshuf = 500;
base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse, date);

run_str = [];
for iRun = 1:length(userun)
    run_str = [run_str num2str(userun(iRun)) '_'];
end

fn_reg = fullfile(reg_base, [date '_' mouse '_' run_str 'RegInfo.mat']);
load(fn_reg);

all_TC = [];
for iRun = 1:length(userun);
    sub_TC = RegInfo.Runs(iRun).tc_VBnpsub;
    seqfile =[date '_' mouse '_run' num2str(userun(iRun)) '_Seqposition.mat'];
    load(fullfile(outDir,'analysis',seqfile));
    ntrials = 0;
    for iCond = 1:nCond+1
        ntrials = length(Seqposition(iCond).ind)+ntrials;
    end
    all_TC = [all_TC; sub_TC(1:ntrials*(nON+nOFF),:)];
end

nCells = size(all_TC,2);

if length(userun) == 1
    seqfile =[date '_' mouse '_run' num2str(userun) '_Seqposition.mat'];
    load(fullfile(outDir,'analysis',seqfile));
    Big_Seqposition = Seqposition;
else
    seqfile = [date '_' mouse '_run' num2str(userun) '_Big_Seqposition.mat'];
    load(fullfile(outDir,'analysis',seqfile));
end

begin = 1;
start = 1;
tc_sorted = zeros(size(all_TC));
for iCond = 1:nCond;
    nRep = length(Big_Seqposition(iCond).ind);
    for iRep = 1:nRep;
        ind = Big_Seqposition(iCond).ind(iRep);
        if begin-1+(nOFF+nON)+((ind-1)*(nOFF+nON)) > size(all_TC,1);
            tc_sorted(start:start-1+((nOFF+nON)/nPlanes)-begin+1,:) = all_TC(begin+((ind-1)*((nOFF+nON)/nPlanes)):((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)),:);
            tc_sorted(start+((nOFF+nON)/nPlanes)-begin+1:start+((nOFF+nON)/nPlanes)-1,:) = all_TC(1:begin-1,:);
        else
        tc_sorted(start:start-1+((nOFF+nON)/nPlanes),:) = all_TC(begin+((ind-1)*((nOFF+nON)/nPlanes)):begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)),:);
        end
        start = start+((nOFF+nON)/nPlanes);
    end
end

%resort blanks
nblanks = length(Big_Seqposition(end).ind);
for iblank = 1:nblanks;
    ind = Big_Seqposition(end).ind(iblank);
        if begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)) > size(all_TC,1);
            tc_sorted(start:start-1+((nOFF+nON)/nPlanes)-begin+1,:) = all_TC(begin+((ind-1)*((nOFF+nON)/nPlanes)):((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)),:);
            tc_sorted(start+((nOFF+nON)/nPlanes)-begin+1:start+((nOFF+nON)/nPlanes)-1,:) = all_TC(1:begin-1,:);
        else
        tc_sorted(start:start-1+((nOFF+nON)/nPlanes),:) = all_TC(begin+((ind-1)*((nOFF+nON)/nPlanes)):begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)),:);
        end
    start = start+((nOFF+nON)/nPlanes);
end
if start<size(tc_sorted,3);
    tc_sorted(start:end,:)=[];
end

fn_tc = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_tcsorted.mat']);
save(fn_tc, 'tc_sorted');

stim_reps = zeros(1,nCond+1);
for iCond = 1:nCond;
    stim_reps(1,iCond) = length(Big_Seqposition(iCond).ind);
    stim_reps(1,end) = length(Big_Seqposition(end).ind);
end

nReps = sum(stim_reps);
fn2 = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
save(fn2, 'stim_reps');


cell_resp = zeros(nOFF+nON,nCells,nReps);
cell_resp_avg = zeros(nOFF+nON,nCells,length(stim_reps));
start = 1;
rep = 1;
for iCond = 1:length(stim_reps); 
    nRep = stim_reps(iCond);
    for iRep = 1:nRep;
        cell_resp(:,:,rep) = tc_sorted(start:start-1+nON+nOFF,:);
        start = start+nON+nOFF;
        rep = rep+1;
    end
    cell_resp_avg(:,:,iCond) = mean(cell_resp(:,:,rep-10:rep-1),3);
end
resp_off = squeeze(mean(cell_resp(pre_win(1):pre_win(2),:,:),1));
resp_on = squeeze(mean(cell_resp(post_win(1):post_win(2),:,:),1));
resp_dF = bsxfun(@minus, resp_on, resp_off);
resp_dFoverF = bsxfun(@rdivide, resp_dF, resp_on);

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
        save(fn_out, 'resp_off', 'resp_on', 'resp_dF', 'resp_dFoverF');

% cell_resp_avg_norm = bsxfun(@minus, cell_resp_avg,mean(cell_resp_avg,1));
% 
% cell_resp_avg_on = squeeze(mean(cell_resp_avg(post_win(1):post_win(2),:,:),1));
% cell_resp_avg_off = squeeze(mean(cell_resp_avg(pre_win(1):pre_win(2),:,:),1));
% cell_dF = bsxfun(@minus,cell_resp_avg_on, cell_resp_avg_off);
% cell_dFoverF = bsxfun(@rdivide,cell_dF, cell_resp_avg_off);
% 
alphaB = .05./(25);
Info_ttest_mat = zeros(nCells-1,25);

for iCell = 1:nCells
    start = 1;
    p_ttestB = zeros(1,25);
    for iCond = 1:(length(stim_reps)-1)
        nRep = stim_reps(1,iCond);
        [h_ttestB1,p_ttestB1] = ttest(resp_off(iCell,start:start-1+nRep),resp_on(iCell,start:start-1+nRep),alphaB,'left');
        p_ttestB(1,iCond) = p_ttestB1;
        start = start+nRep;
    end
    Info_ttest_mat(iCell,:) = p_ttestB;
end

Cell_ttest = min(Info_ttest_mat,[],2) < alphaB;
cell_ind = find(Cell_ttest==1);
nCells = length(Cell_ttest_pass);

SF_vec0 = [.32 .16 .08 .04 .02]; %flipped to have low to high SF in square  %flipud
TF_vec0 = [1 2 4 8 15];

[tftf,sfsf] = meshgrid(TF_vec0,SF_vec0); 
grid2.sfsf = sfsf;
grid2.tftf = tftf;

dSF = median(diff(log2(SF_vec0)));
dTF = median(diff(log2(TF_vec0)));
SF_vec00 = log2(SF_vec0(1)):(dSF/10):log2(SF_vec0(end));
TF_vec00 = log2(TF_vec0(1)):(dTF/10):log2(TF_vec0(end));
[sfsf00,tftf00]=meshgrid(SF_vec00,TF_vec00);
grid2.sfsf00 = sfsf00;
grid2.tftf00 = tftf00;

Ind_struct = [];
start = 1;
for iCond = 1:nCond
    nRep = stim_reps(iCond);
    Ind_struct(iCond).all_trials = [start:start-1+nRep];
    start = start+nRep;
end

Fit_struct = [];
for count_shuf = 0:Nshuf
    fprintf('.')
    Im_mat_USE = zeros(nCells, 25);
    Im_mat_std = zeros(nCells, 25);
    dF_mat = zeros(nCells, 25);
    for iCond = 1:25        
        ind_all = Ind_struct(iCond).all_trials;
        if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
            ind_all_1 = ind_all(randsample(length(ind_all),length(ind_all),1));
        else
            ind_all_1 = ind_all;        
        end
        Im_mat_USE(:,iCond) = mean(resp_dFoverF(cell_ind,ind_all_1),2);
        dF_mat(:,iCond) = mean(resp_dF(cell_ind,ind_all_1),2);
    end

    start = 1;
    for iCell = 1:nCells;
        a = Im_mat_USE(iCell,:);
        if max(a,[],2) > 0     
            if min(a,[],2) < 0
                min_resp = min(a,[],2);
                b = a-min_resp;
            else
                b= a;
            end
            c = reshape(b',length(SF_vec0),length(TF_vec0));
            %b2 = b( ind_SFuse(:,1),ind_TFuse(:,1));
            data = c';
%             ind0 = find(data<0);
%             data(ind0) = NaN;    
            if count_shuf == 0
                PLOTIT_FIT = 0;
                SAVEALLDATA = 1;
                Fit_2Dellipse_LG_soma
                eval(['Fit_struct(iCell).True.s_',' = s;']);
            else
                SAVEALLDATA = 0;
                PLOTIT_FIT = 0;
                Fit_2Dellipse_LG_soma
                eval(['Fit_struct(iCell).Shuf(count_shuf).s_',' = s;']);
            end
        end               
    end
end

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct.mat']);   
        save(fn_out, 'Fit_struct')


if Nshuf>1;
    for iCell = 1:nCells
        if ~isempty(Fit_struct(iCell).True)                
            eval(['tmp = Fit_struct(iCell).True.s_.x;']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.SFhicut_50];']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.TFhicut_50];']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.SFhicut_10];']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.TFhicut_10];']);
            fit_true_vec(iCell,:) = tmp;
        end
    end

    for count_shuf = 1:Nshuf
        for iCell = 1:nCells
            if ~isempty(Fit_struct(iCell).Shuf)
                eval(['tmp = Fit_struct(iCell).Shuf(count_shuf).s_.x;']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.SFhicut_50];']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.TFhicut_50];']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.SFhicut_10];']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.TFhicut_10];']);
                %fit is: %A sigma_SF sigma_TF sf0 tf0 xi
                fit_shuf_vec(iCell,:,count_shuf) = tmp;
            end
        end
    end

    Npars = size(fit_shuf_vec,2);
    lbub_fits = zeros(nCells,Npars,5);
    alpha_bound = .025;
    for iCell = 1:nCells
        for count2 = 1:Npars
            tmp = squeeze(fit_shuf_vec(iCell,count2,:));
            [i,j] = sort(tmp);
            ind_shuf_lb = ceil(Nshuf*alpha_bound);
            ind_shuf_ub = ceil(Nshuf*(1-alpha_bound));
            lbub_fits(iCell,count2,1) = i(ind_shuf_lb);
            lbub_fits(iCell,count2,2) = i(ind_shuf_ub);
            lbub_fits(iCell,count2,3) = mean(i); 
            lbub_fits(iCell,count2,5) = std(i);
        end
        %now take means from truedata fit:
        lbub_fits(iCell,:,4) = fit_true_vec(iCell,:);
    end
end

lbub_diff = lbub_fits(:,:,2)-lbub_fits(:,:,1);

goodfit_ind = [];
for iCell = 1:nCells
    if lbub_diff(iCell,4)<2 
        if lbub_diff(iCell,5)<2
            goodfit_ind = [goodfit_ind iCell];
        end
    end
end

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);   
save(fn_out, 'lbub_fits', 'lbub_diff', 'goodfit_ind')

anal_base = '\\Zmey\storlab\users\Lindsey\DataShare\Layer5\';    
figure;
start = 1;
fig = 1;
for iCell = 1:nCells
    if start==37
        fn_out = fullfile(anal_base, mouse, [date '_' mouse '_' num2str(userun) '_data_fits_figure' num2str(fig) '.pdf']);
        print(gcf, '-dpdf', fn_out);figure;
        start = 1;
        fig = fig+1;
    end
    subplot(6,6,start)
    imagesq(reshape(Fit_struct(iCell).True.s_.orig,5,5)');
    title(num2str(max(max(Fit_struct(iCell).True.s_.orig,[],1),[],2)));
    subplot(6,6,start+1)
    imagesq(Fit_struct(iCell).True.s_.k2b_plot);
    colormap(gray)
    if find(goodfit_ind==iCell);
        title('**')
    end
    %title(num2str(max(cell_dFoverF(iCell,:),[],2)));
    start = start+2;
end

