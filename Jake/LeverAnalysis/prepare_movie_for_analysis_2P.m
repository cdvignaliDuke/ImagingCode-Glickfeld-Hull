clear;
%  go over all the movies in a session
% 1. select an ROI
% 2. stack all the movies to one tiff file

base_dir =  '\\CRASH.dhe.duke.edu\data\home\jake\';
SubNum = '925';
date = '150704';
run = '_000_000';
time = '1937';
mouse = 'img25';
image_dest  = fullfile(base_dir, [date '_' mouse], mouse);
BIN_SIZE =1;

out_base = 'Z:\home\lindsey\Analysis\2P\Jake';
run_name = [date '_' mouse '_run' run(length(run)-2:end)];
out_path = fullfile(out_base,run_name);
% load MWorks file
 behav_dir = [base_dir '\Data\2P_imaging\behavior'];
 cd(behav_dir);
 mworks = ['data-' 'i' SubNum '-' date '-' time '.mat']; 
 load (mworks);

% find sbx info files
data_dir = [base_dir '\Data\2P_imaging\' date '_' mouse '\' mouse];
cd(data_dir);
fName = [mouse run];
imgMatFile = [fName '.mat'];
load(imgMatFile);

[frame_times frame_count input] = get_frame_time_by_counters(input, info);

dest =  fullfile(out_path,run_name);
save([dest '_frame_times.mat'],  'frame_times', 'input');

%load and register dataset
nframes = info.config.frames;
tic
data = sbxread(fName,0,nframes);
toc
data = squeeze(data);

%remove negative data (by addition)
data_sub = data-min(min(min(data,[],1),[],2),[],3);
clear data

% register
data_avg = mean(data_sub(:,:,90:190),3);
figure; imagesq(data_avg); colormap(gray)

data_std = std(double(data_sub(:,:,1:1000)),[],3);
% use first file to calculate ROI
[ROI_x, ROI_y] = get_2P_ROI(data_std); % get the ROI -  must be a rectangle   

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

dest =  fullfile(out_path,run_name);
save([dest '_data_reg.mat'],  'data_reg');

img = data_reg(ROI_x,ROI_y,:);
writetiff(img,[dest '_ROI.tif']);

img_down = stackGroupProject(img,10);

%prep for pca
global stack
stack = single(img_down);
defaultopts = {'nComp',300,'BorderWidth',4};
options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);
[ny,nx,nt]=size(stack);
roi = imborder([ny,nx],options.BorderWidth,0); 
fprintf('Masking edges... ');
stack= bsxfun(@times,stack,single(roi));
fprintf('Done\n');
% compute thin svd using randomized algorithm
pcs = stackFastPCA(1,options.nComp);
% save extracted components 
fprintf('Saving principal components\n');
save([dest '_pca_usv.mat'],'-struct','pcs');

%visualize pca components
nt = size(pcs.v,1);
figure;
sm = stackFilter(pcs.U,1.5);
ax=[];
for pc = 1:25;                   % in order to visualize additional PCs simply alter the range (e.g. 26:50) Then subtract the appropriate amount from pc in the next line
    ax(pc)=subplot(5,5,pc);
    imagesc(sm(:,:,[pc]));
end;
colormap gray;

%compute independent components
PCuse = [1:100];
mu = 0;
nIC = 32;
termtol = 1e-6;
maxrounds = 400;
mixedsig = pcs.v';
mixedfilters = pcs.U;
CovEvals = diag(pcs.s).^2;
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, ...
    mixedfilters, CovEvals, PCuse, mu, nIC,[],termtol,maxrounds);

dt = 1/frGetFrameRate;
tt = [0:nt-1]/frGetFrameRate;


%% TC amd ROI code
cs = permute(ica_filters,[2,3,1]);
sm = stackFilter(cs,1.5);
figure;
ind = 1;
sel = [1:32];    
for ic = sel
    subplot(8,4,ind);                 %change here too
    imstretch(sm(:,:,ic),[.5 .99],1.5);
    ind = ind+1;
    text(.8,.1,num2str(ic),'fontsize',12,'color','w','fontweight','bold','unit','norm');
end;

save([dest '_ICs.mat'],'sm');
%% Start here!!
load([dest '_ICs.mat']);
img = readtiff([dest '_ROI.tif']);
img_down = stackGroupProject(img,10);
%segment from ICs
nIC = 32;
sel = [1:nIC];  
mask_cell = zeros(size(sm));
for ic = sel
    bwimgcell = imCellEditInteractive(sm(:,:,ic),[]);
    mask_cell(:,:,ic) = bwlabel(bwimgcell);
    close all
end

sz = size(img);
mask_cell_temp = reshape(mask_cell,[sz(1)*sz(2) nIC]);
for ic = sel
    ind = find(mask_cell_temp(:,ic));
    mask_cell_temp(ind,ic)=1;
end
mask_cell_temp = reshape(mask_cell_temp,[sz(1)*sz(2) nIC]);

data_tc = zeros(size(img_down,3), nIC);
for ic = sel;
    if sum(mask_cell_temp(:,ic),1)>0
        data_tc(:,ic) = stackGetTimeCourses(img_down, reshape(mask_cell_temp(:,ic), [sz(1) sz(2)]));
    end
end
data_corr = corrcoef(data_tc);


mask_all = zeros(1,sz(1)*sz(2));
count = 0;
for ic = 1:nIC
    ind_new = find(mask_cell_temp(:,ic))';
    if length(ind_new)>1
        ind_old = find(mask_all);
        overlap = ismember(ind_old,ind_new);
        ind_both = find(overlap);
        if length(ind_both)>1
            ic_match = unique(mask_all(ind_old(ind_both)));
            for im = 1:length(ic_match)
                if data_corr(ic, ic_match(im))> 0.8
                    count = count+1;
                    mask_all(ind_new) = ic_match(im);
                else
                    mask_all(ind_new) = ic;
                    mask_all(ind_old(ind_both)) = 0;
                end
            end
        else
             mask_all(ind_new) = ic;
        end
    end
end
figure; imagesc(reshape(mask_all,[sz(1) sz(2)]))

start = 1;
mask_final = zeros(size(mask_all));
for ic = 1:max(mask_all,[],2)
    ind = find(mask_all==ic);
    if length(ind)>0
        mask_final(ind)=start;
        start= start+1;
    end
end

figure; imagesc(reshape(mask_final,[sz(1) sz(2)]))
print([dest '_mask_final.eps'], '-depsc');
print([dest '_mask_final.pdf'], '-dpdf');

data_tc = stackGetTimeCourses(img, reshape(mask_final,[sz(1) sz(2)]));
save([dest '_ROI_TCs.mat'],'data_tc', 'mask_final');
