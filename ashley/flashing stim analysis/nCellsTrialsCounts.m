% cells
ncells_expt{1}(i) = cellfun(@(x) size(x,2),hitVal_short(i));
ncells_expt{2}(i) = cellfun(@(x) size(x,2),missVal_short(i));
ncells_expt{3}(i) = cellfun(@(x) size(x,2),allVal_short(i));

ncells_expt{4}(i) = cellfun(@(x) size(x,2),hitVal_long(i));
ncells_expt{5}(i) = cellfun(@(x) size(x,2),missVal_long(i));
ncells_expt{6}(i) = cellfun(@(x) size(x,2),allVal_long(i));

ncells_expt{7}(i) = cellfun(@(x) size(x,2),hitVal_all(i));
ncells_expt{8}(i) = cellfun(@(x) size(x,2),missVal_all(i));
ncells_expt{9}(i) = cellfun(@(x) size(x,2),allVal_all(i));

ncells_expt{10}(i) = cellfun(@(x) size(x,2),hitInv_short(i));
ncells_expt{11}(i) = cellfun(@(x) size(x,2),missInv_short(i));
ncells_expt{12}(i) = cellfun(@(x) size(x,2),allInv_short(i));

ncells_expt{13}(i) = cellfun(@(x) size(x,2),hitInv_long(i));
ncells_expt{14}(i) = cellfun(@(x) size(x,2),missInv_long(i));
ncells_expt{15}(i) = cellfun(@(x) size(x,2),allInv_long(i));

ncells_expt{16}(i) = cellfun(@(x) size(x,2),hitInv_all(i));
ncells_expt{17}(i) = cellfun(@(x) size(x,2),missInv_all(i));
ncells_expt{18}(i) = cellfun(@(x) size(x,2),allInv_all(i));

% trials
ntrials_expt{1}(i) = sum(cellfun(@length,hitVal_short_ind));
ntrials_expt{2}(i) = sum(cellfun(@length,missVal_short_ind));
ntrials_expt{3}(i) = sum(cellfun(@length,allVal_short_ind));

ntrials_expt{4}(i) = sum(cellfun(@length,hitVal_long_ind));
ntrials_expt{5}(i) = sum(cellfun(@length,missVal_long_ind));
ntrials_expt{6}(i) = sum(cellfun(@length,allVal_long_ind));

ntrials_expt{7}(i) = sum(cellfun(@length,hitVal_all_ind));
ntrials_expt{8}(i) = sum(cellfun(@length,missVal_all_ind));
ntrials_expt{9}(i) = sum(cellfun(@length,allVal_all_ind));

ntrials_expt{10}(i) = sum(cellfun(@length,hitInv_short_ind));
ntrials_expt{11}(i) = sum(cellfun(@length,missInv_short_ind));
ntrials_expt{12}(i) = sum(cellfun(@length,allInv_short_ind));

ntrials_expt{13}(i) = sum(cellfun(@length,hitInv_long_ind));
ntrials_expt{14}(i) = sum(cellfun(@length,missInv_long_ind));
ntrials_expt{15}(i) = sum(cellfun(@length,allInv_long_ind));

ntrials_expt{16}(i) = sum(cellfun(@length,hitInv_all_ind));
ntrials_expt{17}(i) = sum(cellfun(@length,missInv_all_ind));
ntrials_expt{18}(i) = sum(cellfun(@length,allInv_all_ind));

