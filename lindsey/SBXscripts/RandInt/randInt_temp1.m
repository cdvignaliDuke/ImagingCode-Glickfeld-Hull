h_Bg90 = zeros(nCells,1);
p_Bg90 = zeros(nCells,1);
h_Bl90 = zeros(nCells,1);
p_Bl90 = zeros(nCells,1);
data_dfof_trial = bsxfun(@minus, data_dfof,mean(data_dfof(19:24,:,:,:)));
nCells = size(data_dfof,2);
base_resp_avg = zeros(nCells,1);
targ_resp_avg = zeros(nCells,1);
for i = 1:length(good_resp_ind)
    iC = good_resp_ind(i);
    base_ind = find(baseDir == 0);
    targ_ind = intersect(base_ind, find(targetDelta == 90));
    base_resp = squeeze(mean(data_dfof(29:34,iC,1,base_ind),1)-mean(data_dfof(19:24,iC,1,base_ind),1));
    targ_resp = squeeze(mean(data_dfof(29:34,iC,6,targ_ind),1)-mean(data_dfof(19:24,iC,6,targ_ind),1));
    [h_Bg90(iC,:) p_Bg90(iC,:)] = ttest2(targ_resp, base_resp, 'tail','left');
    [h_Bl90(iC,:) p_Bl90(iC,:)] = ttest2(targ_resp, base_resp, 'tail','right');
    base_resp_avg(iC,:) = mean(base_resp,1);
    targ_resp_avg(iC,:) = mean(targ_resp,1);
end

base_resp_rect = base_resp_avg;
base_resp_rect(base_resp_avg<0) = 0;
targ_resp_rect = targ_resp_avg;
targ_resp_rect(targ_resp_avg<0) = 0;
diff_rat = base_resp_rect-targ_resp_rect./(base_resp_rect+targ_resp_rect);
ind0 = find(h_Bg90);
ind90 = find(h_Bl90);
indint = intersect(find((h_Bg90+h_Bl90)==0), good_resp_ind);
ind_str{1} = 'ind0'; 
ind_str{2} = 'indint';
ind_str{3} = 'ind90';

figure; 
for i = 1:3
    subplot(3,1,i)
    hist((max_dir(eval(ind_str{i}),:)-1)*30,[-30:30:180])
    title([ind_str{i} '; n = ' num2str(length(eval(ind_str{i}))) ' cells'])
    xlabel('Ori (deg)')
end
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriHist_bysig.pdf']),'-dpdf')

figure;
for i = 1:3
    subplot(3,1,i)
    for idir = 1:ndir
        ind = find(baseDir == dirs(idir));
        resp = squeeze(mean((mean(data_dfof(29:34,eval(ind_str{i}),1,ind),1)-mean(data_dfof(19:24,eval(ind_str{i}),1,ind),1)),4));
        errorbar(dirs(idir), mean(resp,2), std(resp,[],2)./length(eval(ind_str{i})),'o')
        hold on
    end
    xlim([-30 180])
    xlabel('Ori (deg)')
    title([ind_str{i} '; n = ' num2str(length(eval(ind_str{i}))) ' cells'])
end
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning_bysig.pdf']),'-dpdf')

figure
for i = 1:3
    subplot(3,1,i)
    hist(diff_rat(eval(ind_str{i})), [-1:.1:1])
    title([ind_str{i} '; n = ' num2str(length(eval(ind_str{i}))) ' cells'])
    xlim([-1.1 1.1])
end
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_rat_bysig.pdf']),'-dpdf')

ind0 = find(diff_rat>0);
ind90 = find(diff_rat<-.5);
indint = intersect((find(diff_rat>-.5)),(find(diff_rat<0)));