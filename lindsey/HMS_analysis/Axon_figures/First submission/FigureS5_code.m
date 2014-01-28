fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 'S5';

P = 2;
matrix = 'SF5xTF5';
inj = 'LM';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

areas = strvcat('PM', 'AL', 'V1');
col = strvcat('c', 'r', 'k');


%intrinsic ret
fn = 'G:\users\lindsey\analysisLG\active mice\Y28\111012_Y28_run3_dF_F_Nstim6_POST_12_19_MEAN.tif';
ret = readtiff(fn);

for ipos = 2:3
figure;
imagesq(ret(:,:,ipos));
colormap(gray)
colorbar
caxis([-0.01 0.01]);
fn_out = fullfile(fig_base, ['Figure' fig], [mouse '_111012_ret_pos' num2str(ipos) '.ps']);
print(gcf, '-depsc2', fn_out);
end



%GCAMP ret
mouse = 'Y28';
userun = '1';
date = '111108';
stim = '6pos';

base = 'G:\users\lindsey\analysisLG\active mice';    
outDir = fullfile(base, mouse,date);
FOV = readtiff(fullfile(outDir,[date '_' mouse '_run1_FOV.tif']));

fn_data = fullfile(outDir,[date '_' mouse '_run1_' stim '.tif']);
ret_stack_1 = readtiff(fn_data);
[max_stack1 ind_max1] = max(ret_stack_1,[],3); 

fn_data = fullfile(outDir,[date '_' mouse '_run2_' stim '.tif']);
ret_stack_2 = readtiff(fn_data);
[max_stack2 ind_max2] = max(ret_stack_2,[],3); 
max_stack = max(ret_stack_2,[],3)-1;
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

%widefield islands

fn = 'G:\users\lindsey\analysisLG\active mice\Y28\111117\111117_Y28_run1_FOV.tif';
islands = readtiff(fn);

figure;
imagesq(islands);
colormap(gray);
caxis([200 1800]);
colorbar;

fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_' mouse '_111117_FOV.ps']);
        print(gcf, '-depsc2', fn_out);
        

