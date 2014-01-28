date = '120525';
mouse = 'CM63';
userun = [1:4];
nON = 5;
nOFF = 5;
pre_win = [3 5];
post_win = [6 10];
base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse, date);

nCells = size(RegInfo.All.tc_USE,2);
fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
        load(fn);
nReps = sum(stim_reps);

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
nPlanes = 1;
nCond = 25;
tc_sorted = zeros(size(RegInfo.All.tc_USE));
for iCond = 1:nCond;
    nRep = length(Big_Seqposition(iCond).ind);
    for iRep = 1:nRep;
        ind = Big_Seqposition(iCond).ind(iRep);
        if begin-1+(nOFF+nON)+((ind-1)*(nOFF+nON)) > size(RegInfo.All.tc_USE,1);
            tc_sorted(start:start-1+((nOFF+nON)/nPlanes)-begin+1,:) = RegInfo.All.tc_USE(begin+((ind-1)*((nOFF+nON)/nPlanes)):((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)),:);
            tc_sorted(start+((nOFF+nON)/nPlanes)-begin+1:start+((nOFF+nON)/nPlanes)-1,:) = RegInfo.All.tc_USE(1:begin-1,:);
        else
        tc_sorted(start:start-1+((nOFF+nON)/nPlanes),:) = RegInfo.All.tc_USE(begin+((ind-1)*((nOFF+nON)/nPlanes)):begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)),:);
        end
        start = start+((nOFF+nON)/nPlanes);
    end
end

%resort blanks
nblanks = length(Big_Seqposition(end).ind);
for iblank = 1:nblanks;
    ind = Big_Seqposition(end).ind(iblank);
        if begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)) > size(RegInfo.All.tc_USE,1);
            tc_sorted(start:start-1+((nOFF+nON)/nPlanes)-begin+1,:) = RegInfo.All.tc_USE(begin+((ind-1)*((nOFF+nON)/nPlanes)):((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)),:);
            tc_sorted(start+((nOFF+nON)/nPlanes)-begin+1:start+((nOFF+nON)/nPlanes)-1,:) = RegInfo.All.tc_USE(1:begin-1,:);
        else
        tc_sorted(start:start-1+((nOFF+nON)/nPlanes),:) = RegInfo.All.tc_USE(begin+((ind-1)*((nOFF+nON)/nPlanes)):begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)),:);
        end
    start = start+((nOFF+nON)/nPlanes);
end
if start<size(tc_sorted,3);
    tc_sorted(start:end,:)=[];
end

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


cell_resp_avg_norm = bsxfun(@minus, cell_resp_avg,mean(cell_resp_avg,1));

cell_resp_avg_on = squeeze(mean(cell_resp_avg(post_win(1):post_win(2),:,:),1));
cell_resp_avg_off = squeeze(mean(cell_resp_avg(pre_win(1):pre_win(2),:,:),1));
cell_dF = bsxfun(@minus,cell_resp_avg_on, cell_resp_avg_off);
cell_dFoverF = bsxfun(@rdivide,cell_dF, cell_resp_avg_off);

for iCell = 1:nCells
    figure;
    subplot(1,2,1)
    resp_sq = reshape(cell_dFoverF(iCell,1:nCond),[5 5])';
    imagesq(resp_sq);
    colormap(gray)
    if Cell_ttest(iCell,:) == 1
        title('**')
    end
    subplot(1,2,2)
    plot(tc_sorted(:,iCell));
end

figure;
for iCell = 1:nCells
    subplot(5,5,iCell)
    resp_sq = reshape(cell_dFoverF(iCell,1:nCond),[5 5])';
    imagesq(resp_sq);
    colormap(gray)
    title(num2str(max(cell_dFoverF(iCell,:),[],2)));
end

tc_sorted_norm = bsxfun(@minus,tc_sorted,mean(tc_sorted,1));

figure;
for iCell = 1:9
    subplot(3,3,iCell)
    tcoffsetplot(squeeze(cell_resp_avg_norm(:,iCell,:)));
end

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