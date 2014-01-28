mouse = 'CM10';
date = '120517';
userun = 1:4;

base = 'G:\users\lindsey\analysisLG\active mice';    
outDir = fullfile(base, mouse,date);

fn_local = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
load(fn_local);
        
fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
load(fn_out);

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct.mat']);
load(fn_out);

start = 1;
fig = 1;
for iCell = 1:n_pix
    if start == 37
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_bouton_tuning_fig' num2str(fig) '.pdf']);
        print(gcf, '-dpdf', fn_out);
        figure;
        fig = fig+1;
        start = 1;
    end
    subplot(6,6,start)
    imagesq(Fit_struct(iCell).True.s_.orig);
    title(num2str(max(max(Fit_struct(iCell).True.s_.orig,[],2),[],1)));
    subplot(6,6,start+1)
    imagesq(Fit_struct(iCell).True.s_.k2b_plot);
    if find(goodfit_ind == iCell)
        title('**');
    end
    colormap(gray)
    start = start+2;
end
fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_bouton_tuning_fig' num2str(fig) '.pdf']);
        print(gcf, '-dpdf', fn_out);

figure;
imagesq(local_max)
fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_local_maxima.pdf']);
        print(gcf, '-dpdf', fn_out);
        
good_local = local_max;
for iCell = 1:n_pix
    if find(goodfit_ind == iCell)
        good_local(Fit_struct(iCell).True.s_.ypos, Fit_struct(iCell).True.s_.xpos)= 5;
    end
end


fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_good_local_maxima.pdf']);
        print(gcf, '-dpdf', fn_out);
        
figure;
for iCell = 1:n_pix
    if find(goodfit_ind== iCell)
        plot(2^Fit_struct(iCell).True.s_.x(5),2^Fit_struct(iCell).True.s_.x(4), 'o');
        hold on;
    end
end
fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_SFTF_scatter.pdf']);
        print(gcf, '-dpdf', fn_out);

speed_local = zeros(size(local_max));
for iCell = 1:n_pix
    if find(goodfit_ind == iCell)
        speed_local(Fit_struct(iCell).True.s_.ypos, Fit_struct(iCell).True.s_.xpos)= (2^Fit_struct(iCell).True.s_.x(5)./2^Fit_struct(iCell).True.s_.x(4))*10;
    end
end
figure; imagesq(log2(speed_local));
fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_speed_local_maxima.pdf']);
        print(gcf, '-dpdf', fn_out);
        
fn_sorted = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_sorted.tif']);
stack = readtiff(fn_sorted);

fn_reps = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
load(fn_reps);

siz = size(stack);
stack_tc = zeros(siz(1), siz(2), (nON+nOFF)*26);
begin = 0;
start = 0;
for iCond = 1:26
    nReps = stim_reps(1,iCond);
    stack_tc_temp = zeros(siz(1), siz(2), nON+nOFF, nReps);
    for iRep = 1:nReps
        stack_tc_temp(:,:,:,iRep) = stack(:,:,start+1:start+nON+nOFF);
        start = start+nON+nOFF;
    end
    stack_tc(:,:,begin+1:begin+nON+nOFF) = squeeze(mean(stack_tc_temp,4));
    begin = begin+nON+nOFF;
end

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_avgTC.tif']);
writetiff(stack_tc, fn_out);

