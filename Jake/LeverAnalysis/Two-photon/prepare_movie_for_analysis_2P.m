%  prepares movie for timecourse analysis
% 1. select an image ROI
% 2. register images
% 3. extract PCA/ICA components
% 4. segment dendrite ROIs
% 5. extract timecourses
% 6. subtract neuropil

%% load data and associated info files
% % load MWorks file
 behav_dir = [base_dir 'behavior'];
 cd(behav_dir);
 mworks = ['data-' 'i' subNum '-' date '-' time '.mat']; 
 load (mworks);

% find sbx info files
data_dir = [base_dir date '_' mouse '\' mouse];
cd(data_dir);
fName = [mouse run];
imgMatFile = [fName '.mat'];
load(imgMatFile);

[frame_times frame_count input] = get_frame_time_by_counters(input, info);

dest =  fullfile(out_path,run_name);
dest1 = fullfile(out_path,run1_name);
save([dest '_frame_times.mat'],  'frame_times', 'input');

%load and register dataset
nframes = info.config.frames;
disp('loading sbx image file');
tic
data = sbxread(fName,0,nframes);
toc
data = squeeze(data);

%remove negative data (by addition)
data_sub = data-min(min(min(data,[],1),[],2),[],3);
clear data

%% 1. select ROI 
data_std = std(double(data_sub(:,:,1:1000)),[],3);
% use first file to calculate ROI
if irun == 1
    [ROI_x, ROI_y] = get_2P_ROI(data_std); % get the ROI -  must be a rectangle 
    save([dest '_ROI_xy.mat'],  'ROI_x', 'ROI_y');
else
    load([dest1 '_ROI_xy.mat'])
end

%% 2. register images
%add in an if statement to use the first movie to register the second one.
%Save the data_avg
data_avg = mean(data_sub(:,:,90:190),3);
figure; imagesq(data_avg); colormap(gray)
[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

save([dest '_data_reg.mat'],  'data_reg');

img = data_reg(ROI_x,ROI_y,:);
%writetiff(img,[dest '_ROI.tif']);   %commented out to save time
clear data_reg

%% 3. PCA and ICA
if irun == 1
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
    save([dest '_pca_usv.mat'],'pcs');


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

    PCuse = [1:40];
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

    cs = permute(ica_filters,[2,3,1]);
    sm = stackFilter(cs,1.5);
    figure;
    ind = 1;
    sel = [1:32];    
    for ic = sel
        subplot(8,4,ind);                 %change here too
        %imstretch(sm(:,:,ic),[.5 .99],1.5);
        imagesc(sm(:,:,ic));
        ind = ind+1;
        text(.8,.1,num2str(ic),'fontsize',12,'color','w','fontweight','bold','unit','norm');
    end;

    save([dest '_ICs.mat'],'sm', 'ica_sig');

    %% 4. segment ROIs from ICs

    nIC = 32;
    sel = [1:nIC];  
    mask_cell = zeros(size(sm));
    for ic = sel
        bwimgcell = imCellEditInteractive(sm(:,:,ic),[]);
        mask_cell(:,:,ic) = bwlabel(bwimgcell);
        close all
    end

    %consolidates all ROIs within IC into single ROI
    thresh = 0.8; %correlation threshold for calling two dendrites one thing
    sz = size(img);
    mask_cell_temp = zeros(sz(1)*sz(2), nIC);
    for ic = sel
        if length(unique(reshape(mask_cell(:,:,ic),[1 sz(1)*sz(2)])))>2
            data_tc_temp = stackGetTimeCourses(img_down,mask_cell(:,:,ic));
            data_corr_temp = corrcoef(data_tc_temp);
            ind_rem = 1:length(unique(reshape(mask_cell(:,:,ic),[1 sz(1)*sz(2)])))-1;
            for i = 1:length(unique(reshape(mask_cell(:,:,ic),[1 sz(1)*sz(2)])))-1
                ind = ind_rem(find(data_corr_temp(min(ind_rem,[],2),ind_rem)>thresh));
                if length(ind)>1
                    for ii = ind
                        if i == 1
                            mask_cell_temp(find(mask_cell(:,:,ic)== ii),ic) = 1;
                        else
                            mask_cell_temp(find(mask_cell(:,:,ic)== ii),nIC) = 1;
                        end
                    end
                else
                    if i == 1
                        mask_cell_temp(find(mask_cell(:,:,ic)== i),ic) = 1;
                    else
                        mask_cell_temp(find(mask_cell(:,:,ic)== ind),nIC) = 1;
                    end  
                end
                ind_rem = ind_rem(~ismember(ind_rem,ind));
                if ~isempty(ind_rem)
                    cat(3, mask_cell_temp, zeros(size(mask_cell_temp(:,:,1))));
                    nIC = nIC+1;
                else
                    break
                end
            end
        else
            mask_cell_temp(find(mask_cell(:,:,ic)),ic) = 1;
        end
    end

    %get preliminary timecourses for segregating and grouping ROIs
    data_tc = zeros(size(img_down,3), nIC);
    for ic = sel;
        if sum(mask_cell_temp(:,ic),1)>0
            data_tc(:,ic) = stackGetTimeCourses(img_down, reshape(mask_cell_temp(:,ic), [sz(1) sz(2)]));
        end
    end
    data_corr = corrcoef(data_tc);
    figure; imagesc(data_corr)

    %consolidate timecourses that are highly correlated
    [i j] = find(and(data_corr>thresh,data_corr<1));
    n = size(i,1);
    if n>1
        for ii = 1:(n/2)
            ind = find(mask_cell_temp(:,i(ii)));
            mask_cell_temp(ind,j(ii)) = 1;
            mask_cell_temp(ind,i(ii)) = 0;
        end
    end

    %finds overlapping pixels of ROIs and based on correlations decides whether
    %to group them or to split them- if splitting, then overlapping pixels are
    %eliminated from both ROIs

    mask_overlap = zeros(1,sz(1)*sz(2));
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
                    if data_corr(ic, ic_match(im))> thresh
                        count = count+1;
                        mask_all(ind_new) = ic_match(im);
                    else
                        mask_all(ind_new) = ic;
                        mask_all(ind_old(ind_both)) = 0;
                        mask_overlap(ind_old(ind_both)) = 1;
                    end
                end
            else
                 mask_all(ind_new) = ic;
            end
        end
    end
    figure; imagesc(reshape(mask_all,[sz(1) sz(2)]))


    % removes ICs smaller than 200 pixels, renumbers ROIs so in continuous ascending order
    start = 1;
    mask_final = zeros(size(mask_all));
    for ic = 1:max(mask_all,[],2)
        ind = find(mask_all==ic);
        if length(ind)<200
             mask_overlap(find(mask_all==ic)) = 1;
             mask_all(ind) = 0;
        end
        ind = find(mask_all==ic);
        if length(ind)>0
            mask_final(ind)=start;
            start= start+1;
        end
    end

    figure; imagesc(reshape(mask_final,[sz(1) sz(2)]))

    print([dest '_mask_final.eps'], '-depsc');
    print([dest '_mask_final.pdf'], '-dpdf');
else
    load([dest1 '_ROI_TCs.mat'])
end

%% 5. Extract timecourses
data_tc = stackGetTimeCourses(img, reshape(mask_final,[sz(1) sz(2)]));

sz = size(mask_cell);
save([dest '_ROI_TCs.mat'],'data_tc', 'mask_final', 'sz', 'mask_cell', 'mask_cell_temp', 'data_corr', 'mask_overlap');


%% 6. Neuropil subtraction
%create masks, get timecourses
%COMMENTED OUT 3/16/16 because the variable created here are not being
%used. Also concerned that there is not enough space between dendrites in
%XY for this to be meaningful. Additionally reducing neuropil contamination
%from the Z axis is less meaningful in this prep since the PC's own
%dendrite extends out of the plain of focus in Z

% nCells = max(unique(mask_final),[],2);
% npTC = zeros(size(data_tc));
% 
% buf = 4;
% np = 6;
% neuropil = squeeze(imCellNeuropil(mask_final,buf,np));
% np_mask = zeros(sz(1),sz(2),nCells);
% for i = 1:nCells
% %     ind_both = and(neuropil(:,i),mask_overlap');
% %     neuropil(ind_both,i) = 0;
%     np_mask(:,:,i) = reshape(neuropil(:,i),[sz(1) sz(2)]);
%     npTC(:,i) = stackGetTimeCourses(img,np_mask(:,:,i));
% end
% 
% %get weights by maximizing skew
% ii= 0.01:0.01:1;
% x = zeros(length(ii), nCells);
% tc_avg = tsmovavg(data_tc,'s',1,1);
% np_avg = tsmovavg(npTC,'s',1,1);
% for i = 1:100
%     x(i,:) = skewness(tc_avg-tcRemoveDC(np_avg.*ii(i)));
% end
% [max_skew ind] =  max(x,[],1);
% np_w = 0.01*ind;
% npSubTC = data_tc-bsxfun(@times,tcRemoveDC(npTC),np_w);
% save([dest '_npSubTCs.mat'],'npSubTC',  'neuropil', 'np_w');    
