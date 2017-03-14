clear all
close all
rc = behavConstsAV;
awFSAVdatasets_V1
iexp = 22;
doRegMovie = 0;

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
dirFolder = expt(iexp).dirtuning;
down = 10;
expTime = expt(iexp).dirtuning_time;

fntun = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate,dirFolder);
%% Load dir tuning dataset/watch movie
fName = [dirFolder '_000_000'];

if ~doRegMovie
    
    [input] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,dirFolder,fName);    

elseif doRegMovie
    [input, data] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,dirFolder,fName);    


    % watch registered direction tuning movie
    data_down = stackGroupProject(data,down);
    clear data

    data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
    data = data_sub;
    clear data_sub

    %load previously registered parameters
    load(fullfile(fntun,'reg_img'))
    tic
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    toc

    figure;
    colormap gray
    for iplot = 1:size(data_reg,3)
        imagesc(data_reg(:,:,iplot))
        drawnow
    end
    
     % watch registered behavior movie
     expFolder = expt(iexp).runs(1,:);
     fName_exp = [expFolder '_000_000'];
     
    [input, data] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,expFolder,fName_exp,1000);    
    
    data_sub = data-min(min(min(data,[],1),[],2),[],3);
    data = data_sub;
    clear data_sub

    %load previously registered parameters
    load(fullfile(fntun,'reg_img'))
    tic
    [out, data_reg]=stackRegister(data,data_avg);
    toc
    fnout = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate,expFolder)
%     save(fullfile(fnout,'reg_out.mat'),'out')

    figure;
    colormap gray
    for iplot = 1:1000
        imagesc(data(:,:,iplot))
        drawnow
    end
    
    % after downsampling
    data_down = stackGroupProject(data,down);
    tic
        [out, data_down_reg]=stackRegister(data_down,data_avg);
    toc
    
    figure;
    colormap gray
    for iplot = 1:1000
        imagesc(data_down_reg(:,:,iplot))
        drawnow
    end
end    

%% load all max dF/F images, crop if necessary

load(fullfile(fntun,'dir_max_images.mat'))

fnbx = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate);
load(fullfile(fnbx,'bx_max_images.mat'))

tun_maxDFF = max(dFF_dirmax,[],3);
bx_maxDFF = max(cat(3,start_max,long_max,tar_max),[],3);
figure;colormap gray; setFigParams4Print('portrait')
subplot(2,1,1)
imagesc(tun_maxDFF)
title('direction tuning')
subplot(2,1,2)
imagesc(bx_maxDFF)
title('behavior')
print(fullfile(fnbx,'max_dFF_used4cellSelect'),'-dpdf','-fillpage')

ypix = size(dFF_dirmax,1);
xpix = size(dFF_dirmax,2);
nstim = length(unique(cell2mat(input.tGratingDirectionDeg)));
%% ****select cells
% iterate through each max dF/F for each direction tuning stimulus
mask_cell_dir = zeros(ypix,xpix,nstim);
for istim = 1:nstim
    img_temp = dFF_dirmax(:,:,istim);
    if istim > 1
        img_temp(mask_cell_dir(:,:,istim-1)>0) = 0;
    end
    bwout = imCellEditInteractive(img_temp); % select cells
    if istim > 1
        mask_temp = bwlabel(bwout);
        prev_mask = mask_cell_dir(:,:,istim-1);
        nc = length(unique(prev_mask(:)))-1;
        mask_temp(mask_temp >0) = mask_temp(mask_temp > 0)+nc;
        overlap = intersect(find(mask_temp(:) > 0), find(prev_mask(:) > 0)); % make sure there's no overlapping cells selected
        if ~isempty(overlap) % fix overlap if there is some
            overlap_cell = unique(mask_temp(overlap));
            for icell = 1:length(overlap_cell)
                mask_temp(mask_temp == overlap_cell(icell)) = 0;                
            end
            cells_temp = unique(mask_temp(mask_temp > 0));
            for icell = 1:length(cells_temp)
                mask_temp(mask_temp == cells_temp(icell)) = nc+icell;
            end
        end
        mask_temp = mask_temp +prev_mask; % combine cell masks
        mask_cell_dir(:,:,istim) = mask_temp;
    else
        mask_cell_dir(:,:,istim) = bwlabel(bwout);
    end
end

figure;
for iplot = 1:nstim
    subplot(3,3,iplot)
    imagesc(mask_cell_dir(:,:,iplot))
end

tun_mask = mask_cell_dir(:,:,nstim);

figure; colormap gray; 
imagesc(tun_mask);
title([num2str(length(unique(tun_mask(:)))-1) ' cells'])

% try with max across bx
% ceil_tun_max = max(tun_maxDFF(:));

img_temp = bx_maxDFF;
% img_temp(1,1) = ceil_tun_max;
img_temp(tun_mask > 0) = 0;
bwout = imCellEditInteractive(img_temp); % select cells
bx_mask = bwlabel(bwout);
nc = length(unique(tun_mask(:)))-1;
bx_mask(bx_mask >0) = bx_mask(bx_mask > 0)+nc;
overlap = intersect(find(bx_mask(:) > 0), find(tun_mask(:) > 0)); % make sure there's no overlapping cells selected
if ~isempty(overlap) % fix overlap if there is some
    overlap_cell = unique(bx_mask(overlap));
    for icell = 1:length(overlap_cell)
        mask_temp(bx_mask == overlap_cell(icell)) = 0;                
    end
    cells_temp = unique(bx_mask(bx_mask > 0));
    for icell = 1:length(cells_temp)
        bx_mask(bx_mask == cells_temp(icell)) = nc+icell;
    end
end

mask_cell = tun_mask+bx_mask;


figure; colormap gray; setFigParams4Print('portrait')
imagesc(mask_cell);
title({[num2str(length(unique(mask_cell(:)))-1) ' cells with behavior'];[mouse '-' expDate]})
print(fullfile(fnbx,'final_mask'),'-dpdf')

save(fullfile(fnbx,'final_mask.mat'),'mask_cell');
