run_str = ['runs-' ImgFolder(1,:)];
if nrun>1
    run_str = [run_str '-' ImgFolder(nrun,:)];
end

clear temp
for irun = 1:nrun
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    temp(irun) = input;
end
input = concatenateDataBlocks(temp);
clear temp
ntrials = sum(input.trialsSinceReset,2);

load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
    
sz_mat = celleqel2mat_padded(input.tGratingDiameterDeg);
szs = unique(sz_mat);
nStim = length(szs);
Stims = szs;
con_mat = celleqel2mat_padded(input.tGratingContrast);
cons = unique(con_mat);

nOn = input.nScansOn;
nOff = input.nScansOff;
nCells = size(npSub_tc,2);
data_mat = zeros(nOn+nOff, nCells, ntrials);
for itrial = 1:ntrials
    data_mat(:,:,itrial) = npSub_tc(1+((itrial-1).*(nOn+nOff)):(itrial.*(nOn+nOff)),:);
end
data_f = mean(data_mat(nOff/2:nOff,:,:),1);
data_df = bsxfun(@minus, data_mat, data_f);
data_dfof = bsxfun(@rdivide, data_df, data_f);

clear data_mat data_f data_df

figure;
if nCells<37
    [n, n2] = subplotn(nCells);
else
    [n, n2] = subplotn(36);
end
tt= (1-nOff:nOn)*(1000./frame_rate);
tuning_mat = zeros(nStim,2,nCells);
tc_mat = zeros(nOn+nOff,nStim,nCells);
start = 1;
f = 1;
cmap = flipud(gray(nStim+1));
Ind_struct = [];
for iCell = 1:nCells
    if start >36
        for i = 1:36
            subplot(n,n2,i)
            ylim([min(min(min(tc_mat(:,:,1:start-1),[],1),[],2),[],3) max(max(max(tc_mat(:,:,1:start-1),[],1),[],2),[],3)])
            vline(nOff)
        end
        suptitle(['Stims: ' num2str(chop(Stims,2))])
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs' num2str(f) '.pdf']), '-dpdf')
        figure;
        f = 1+f;
        start = 1;
    end
    subplot(n, n2, start)
    for iStim = 1:nStim
        ind1 = find(sz_mat == Stims(iStim));
        ind2 = find(con_mat == cons(end));
        ind = intersect(ind1,ind2);
        plot(tt',squeeze(mean(data_dfof(:,iCell,ind),3)), 'Color', cmap(iStim+1,:))
        hold on
        tc_mat(:,iStim,iCell) = squeeze(mean(data_dfof(:,iCell,ind),3));
        tuning_mat(iStim,1,iCell) = mean(mean(data_dfof(nOff+1:nOn+nOff,iCell,ind),3),1);
        tuning_mat(iStim,2,iCell) = std(mean(data_dfof(nOff+1:nOn+nOff,iCell,ind),1),[],3)./sqrt(length(ind));
        Ind_struct(iStim).all_trials = ind;
    end
    start = start+1;
end
for i = 1:start-1
    subplot(n,n2,i)
    ylim([min(min(min(tc_mat,[],1),[],2),[],3) max(max(max(tc_mat,[],1),[],2),[],3)])
end
suptitle([date ' ' mouse ' ' run_str ': ' img_area])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs' num2str(f) '.pdf']), '-dpdf')

figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning' num2str(f) '.pdf']), '-dpdf')
        figure;
        f = 1+f;
        start = 1;
    end
    subplot(n, n2, start)
    errorbar(Stims, tuning_mat(:,1,iCell), tuning_mat(:,2,iCell), '-ok');
    start = start + 1;
end
suptitle([date ' ' mouse ' ' run_str ': ' img_area])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning'  num2str(f) '.pdf']), '-dpdf')
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TuningSummary.mat']), 'tuning_mat', 'tc_mat', 'Ind_struct')

%% size tuning summary
[max_val max_ind] = max(squeeze(tuning_mat(:,1,:)),[],1);
max_sz = szs(max_ind);
max_resp = squeeze(tuning_mat(end,1,:))';
suppIx = max_resp./max_val;
figure;
subplot(2,2,1)
hist(max_sz); xlabel('Size (deg)'); title('Peak size')
subplot(2,2,2)
hist(suppIx); xlabel('MaxSizeResp/MaxResp'); title('Suppression Index')
subplot(2,2,3)
tuning_avg = squeeze(mean(tuning_mat(:,1,:),3));
tuning_sem = squeeze(std(tuning_mat(:,1,:),[],3)./sqrt(size(tuning_mat,3)));
errorbar(szs,tuning_avg, tuning_sem, '-o');
ylabel('dF/F'); xlabel('Size (deg)'); title('Average tuning curve')
subplot(2,2,4)
tuning_mat_norm = bsxfun(@rdivide, squeeze(tuning_mat(:,1,:)),squeeze(max(tuning_mat(:,1,:),[],1))');
tuning_norm_avg = squeeze(mean(tuning_mat_norm,2));
tuning_norm_sem = squeeze(std(tuning_mat_norm,[],2)./sqrt(size(tuning_mat_norm,2)));
errorbar(szs,tuning_norm_avg, tuning_norm_sem, '-o');
ylabel('Normalized dF/F'); xlabel('Size (deg)'); title('Normalized average tuning curve')
suptitle([date ' ' mouse ' ' run_str ': ' img_area])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SizeHist'  num2str(f) '.pdf']), '-dpdf')

%%
figure;
for isz = 1:length(szs)
    ellipse(szs(isz)./2, szs(isz)./2, 0, 0, 0,cmap(isz+1,:));
    hold on
end
axis square
legend(num2str(chop(szs,3)'))