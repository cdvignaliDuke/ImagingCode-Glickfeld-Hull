fn = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_avg.tif']);
stack_dF = readtiff(fn);
sz = size(stack_dF);

stack_dF_slow = zeros(sz(1),sz(2),(nON+nOFF)*3);
stack_dF_med = zeros(sz(1),sz(2),(nON+nOFF)*13);
stack_dF_fast = zeros(sz(1),sz(2),(nON+nOFF)*3);

TF_mat = repmat(uvar(:,1),length(uvar),1);
SF_vec = flipud(uvar(:,2));
SF_mat = zeros(size(TF_mat));
start = 1;
for ivar = 1:length(SF_vec)
    SF_mat(start:start+length(uvar)-1) = repmat(SF_vec(ivar,1),1,5);
    start = start+length(uvar);
end

speed_mat = TF_mat./SF_mat;

[speed_sort ind] = sort(speed_mat);

start = 1;
for iCond = 1:6;
    use = ind(iCond);
    stack_dF_slow(:,:,start:start+nON+nOFF-1) = stack_dF(:,:,1+((use-1)*(nON+nOFF)):((use-1)*(nON+nOFF))+nON+nOFF);
    start = start+nON+nOFF;
end
start=1;
for iCond = 7:19;
    use = ind(iCond);
    stack_dF_med(:,:,start:start+nON+nOFF-1) = stack_dF(:,:,1+((use-1)*(nON+nOFF)):((use-1)*(nON+nOFF))+nON+nOFF);
    start = start+nON+nOFF;
end
start=1;
for iCond = 20:25;
    use = ind(iCond);
    stack_dF_fast(:,:,start:start+nON+nOFF-1) = stack_dF(:,:,1+((use-1)*(nON+nOFF)):((use-1)*(nON+nOFF))+nON+nOFF);
    start = start+nON+nOFF;
end

stack_dF_slow_avg = zeros(sz(1),sz(2),nON+nOFF);
stack_dF_med_avg = zeros(sz(1),sz(2),nON+nOFF);
stack_dF_fast_avg = zeros(sz(1),sz(2),nON+nOFF);
for iFrame = 1:nON+nOFF;
    stack_dF_slow_avg(:,:,iFrame) = mean(stack_dF_slow(:,:,iFrame:nON+nOFF:end),3);
    stack_dF_med_avg(:,:,iFrame) = mean(stack_dF_med(:,:,iFrame:nON+nOFF:end),3);
    stack_dF_fast_avg(:,:,iFrame) = mean(stack_dF_fast(:,:,iFrame:nON+nOFF:end),3);
end

stack_dF_slow_med_fast = cat(3, stack_dF_slow_avg, stack_dF_fast_avg, stack_dF_fast_avg);


fn_out = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_slow_med_fast.tif']);
writetiff(stack_dF_slow_med_fast, fn_out);
