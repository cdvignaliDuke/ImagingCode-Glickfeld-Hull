fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 'S1';

col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

%intrinsic ret
fn = 'G:\users\lindsey\analysisLG\active mice\M14\111026_M14_run3_dF_F_Nstim6_POST_12_19_MEAN.tif';
ret = readtiff(fn);

for ipos = 1:3
figure;
imagesq(ret(:,:,ipos));
colormap(gray)
colorbar
caxis([-0.01 0.01]);
fn_out = fullfile(fig_base, ['Figure' fig], [mouse '_111026_ret_pos' num2str(ipos) '.ps']);
print(gcf, '-depsc2', fn_out);
end

%widefield image
fn = 'G:\users\lindsey\analysisLG\active mice\M14\111128\111128_M14_run2_FOV.tif';
islands = readtiff(fn);

figure;
imagesq(islands);
colormap(gray);
caxis([200 3500]);
colorbar;

fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_' mouse '_111128_FOV.ps']);
        print(gcf, '-depsc2', fn_out);
        
%GCAMP ret
mouse = 'M14';
userun = '2';
date = '111112';
stim = '6pos';

base = 'G:\users\lindsey\analysisLG\active mice';    
outDir = fullfile(base, mouse,date);
FOV = readtiff(fullfile(outDir,[date '_' mouse '_run1_FOV.tif']));

fn_data = fullfile(outDir,[date '_' mouse '_run2_' stim '.tif']);
ret_stack = readtiff(fn_data);
[max_stack ind_max2] = max(ret_stack,[],3); 
max_stack = max(ret_stack,[],3)-1;
max_stack(find(max_stack<0))=0;

cmap = colormap(prism);
cmap = cmap([4 5 6 2 1 3],:);
cmap(2,:) = [0 1 1];

sz = size(ind_max2);
tmp3 = zeros(sz(1),sz(2),3);

for count = 1:6
    [i,j] = find(ind_max2 == count);
    for count2 = 1:length(i);
        tmp3(i(count2),j(count2),:) = cmap(count,:);
    end
end
tmp4 = rgb2hsv(tmp3);

tmp5 = double(FOV) ./ max(max(double(FOV))); %zeros(size(ind_max2));

tmp4(:,:,3) = tmp5;
tmp4(:,:,2) = 1;
tmp4(:,:,1) = tmp4(:,:,1);

tmp6 = hsv2rgb(tmp4);

figure; image(tmp6)
fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_' mouse '_' date '_6pos_ret_baselineF.ps']);
        print(gcf, '-depsc2', fn_out);

tmp5 = double(max_stack) ./ max(max(double(max_stack))); %zeros(size(ind_max2));

tmp4(:,:,3) = tmp5;

tmp6 = hsv2rgb(tmp4);

figure; image(tmp6)
fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_' mouse '_' date '_6pos_ret_maxdF.ps']);
        print(gcf, '-depsc2', fn_out);

ind_max3 = [1 2 3; 4 5 6];

sz = size(ind_max3);
tmp2 = zeros(sz(1),sz(2),3);
tmp3 = zeros(sz(1),sz(2),3);

for count = 1:6
 [i,j] = find(ind_max3 == count);
 for count2 = 1:length(i);
     tmp3(i(count2),j(count2),:) = cmap(count,:);
 end
end

figure; image(tmp3)

fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_' mouse '_' date '_6pos_colormap.ps']);
        print(gcf, '-depsc2', fn_out);
       
%boutons

mouse = 'M14';
date = '111127';
userun = [1:4];
iArea = 1;
iexp = 12;

base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse,date);

fn_stim = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_allstim.mat']);
load(fn_stim);
FOV = mean(stim_off,3);

fn_df =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
stack = readtiff(fn_df);

siz = size(stack);

dF_nodir = zeros(siz(1), siz(2), 25);
start = 1;
for iCond = 1:25
    dF_nodir(:,:,iCond) = mean(stack(:,:,start:start+1),3);
    start = start+2;
    subplot(5,5,iCond)
    imagesq(dF_nodir(:,:,iCond));
end
dF_max = max(dF_nodir,[],3);

ind = [];
n = all_fits(iArea).expt(iexp).n(1);
for iCell = 1:n
    pos = all_fits(iArea).expt(iexp).bouton(iCell).pos;
    if 134<pos(1) & pos(1)<211
        if 99<pos(2) & pos(2)<176
            ind = [ind iCell];
        end
    end
end

x = [101 175];
y = [136 210];

figure
subplot(2,2,1)
imagesq(dF_max);
colormap(gray);
caxis([0 .5]);
hold on
n = length(ind);
xlim(x);
ylim(y);
for iCell = 1:n
    pos = all_fits(iArea).expt(iexp).bouton(ind(iCell)).pos; 
    scatter(pos(2), pos(1), 100, [0.7 0.7 0.7])
    hold on
end
colorbar
subplot(2,2,2)
imagesq(FOV);
colormap(gray);
caxis([50 200]);
hold on
n = length(ind);
xlim(x);
ylim(y);
for iCell = 1:n
    pos = all_fits(iArea).expt(iexp).bouton(ind(iCell)).pos; 
    scatter(pos(2), pos(1), 100, [0.7 0.7 0.7])
    hold on
end
colorbar

fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_M14_FOV_boutons.ps']);
        print(gcf, '-depsc2', fn_out);

        
