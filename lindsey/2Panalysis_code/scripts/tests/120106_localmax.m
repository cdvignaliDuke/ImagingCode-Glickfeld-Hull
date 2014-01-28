mouse = 'Y13';
date = '110509';
userun = [1:2];

base = 'G:\users\lindsey\analysisLG\active mice';    
outDir = fullfile(base, mouse,date);

fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
stack_dF = readtiff(fn_stack);

fn_ttest = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
load(fn_ttest);

stack_max = max(stack_dF,[],3);
f1 = fspecial('average');
stack_max_sm = filter2(f1, stack_max);
stack_max_ttest = ttest_mask.*stack_max;

figure; imagesq(overlay);
siz = size(stack_max);

local_max_sm = zeros(size(stack_max));
for iy = 1:siz(1);
    for ix = 1:siz(2);
        if stack_max_ttest(iy,ix)>0
            sub = stack_max_ttest(iy-3:iy+3,ix-3:ix+3);
            if max(max(sub,[],1),[],2)==stack_max_ttest(iy,ix)
                local_max_sm(iy,ix) = 1;
            end
        end
    end
end

n_pix = sum(sum(local_max));

fn_local = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
save(fn_local, 'local_max');

fn_out = fullfile('\\Zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\120109', [date '_' mouse '_run' num2str(userun) '_Raw_vs_Smooth_Overlay.pdf']);
print(gcf, '-dpdf', fn_out);

local_max_resp = zeros(26, n_pix);
[i, j] = find(local_max ==1);
for ipix = 1:n_pix
    sub_y = [i(ipix)-1:i(ipix)+1]; 
    sub_x = [j(ipix)-1:j(ipix)+1];
    local_max_resp(:,ipix) = squeeze(mean(mean(stack_dF(sub_y, sub_x,:),2),1));
end

local_max_resp_blank = zeros(25,n_pix);
for iCond = 1:25;
local_max_resp_blank(iCond,:) = local_max_resp(iCond,:)- local_max_resp(26,:);
end

figure;
start = 1;
for iresp = 1:n_pix
    if start >100
        figure;
        start = 1;
    end
subplot(10,10,start)
local_max_sq = reshape(local_max_resp_blank(1:25,iresp), [5 5])';
max_resp = max(local_max_resp_blank(:,iresp),[],1);
ax = [0 max_resp];
imagesq(local_max_sq);
colormap(gray)
caxis(ax);
max_dF = num2str(max_resp*100);
% title(max_dF(:,1:3));
start = start+1;
end

local_max_spot = zeros(size(local_max));

figure;
start = 1;
for ipix = 1:n_pix
    if start >100
        figure;
        start = 1;
    end
subplot(10,10,start)
sub_y = [i(ipix)-1:i(ipix)+1]; 
sub_x = [j(ipix)-1:j(ipix)+1];
local_max_spot(sub_y, sub_x) = 1;
imagesq(local_max_spot);
colormap(gray)
start = start+1;
end

fn_resp = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp.mat']);
load(fn_resp);

resp_dF = zeros(n_pix, size(resp_on,3));
for ipix = 1:n_pix
    sub_y = [i(ipix)-1:i(ipix)+1]; 
    sub_x = [j(ipix)-1:j(ipix)+1];
    roi_on = squeeze(mean(mean(resp_on(sub_y,sub_x,:),2),1));
    roi_off = squeeze(mean(mean(resp_off(sub_y,sub_x,:),2),1));
    resp_dF(ipix,:) = ((roi_on-roi_off)./roi_off)';
end
Nshuf = 0;

fn_resp = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp.mat']);
save(fn_resp, 'resp_on', 'resp_off', 'resp_dF')