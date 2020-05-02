%% get path names
% date is day2
date = '200108'; 
ImgFolder = strvcat('003');
time = strvcat('1158');
mouse = 'i1316';
% set align to 1 if aligning any day to day 1. set to 0 if you are just
% working with day 1
alignToRef = 1;
% ref_date is day1
ref_date = '200106';
ref_run = strvcat('003');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
ref_str = catRunName(ref_run, size(ref_run,1));
% Type in your own stuff here:
gl_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\jerry\2P_Imaging_Grace';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\jerry\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';

%% load and register
data = [];
clear temp
trial_n = [];
offset = 0;
%     loading data_dfof_max and data_reg_avg from day 1 (_reg)
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']))
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    data_dfof_max_reg = data_dfof_max;
    data_avg_reg = data_reg_avg;
%     loading data_dfof_max and data_reg_avg from day 1 (_ref)
    load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_stimActFOV.mat']))
    load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']))
    data_dfof_max_ref = data_dfof_max;
    data_avg_ref = data_reg_avg;

%% Thresholding    
%     defining day 1 and day 2 data_reg_avg and taking out the brightest
%     points of the image
    reg = data_avg_reg;
    ref = data_avg_ref;
    reg(find(reg>2000)) = 0;
    reg = (reg./max(max(abs(reg))));
    ref(find(ref>2000)) = 0;
    ref = (ref./max(max(abs(ref)))); 
    sz_target  = size(reg); % i moved this from the select points section below (line 52)
    figure(1);clf;imshowpair(reg, ref, 'montage'); %displays reg and ref images side by side

%% Transforming
%     selecting the matching landmarks on day 1 and day 2 data_reg_avg and
%     shifting the position of day 2 data_reg_avg and data_dfof_max accordingly
    %% Select Points
    sz_target  = size(reg);
    [input_points, base_points] = cpselect(double(reg),double(ref),'Wait', true);
    save('day2Inputs.mat','input_points'); % saving input points
    save('day2Base.mat','base_points'); % saving base points
    %input_points = load('day2Inputs.mat'); %load presaved input coordinates
    %base_points = loadi('day2Base.mat'); %load presaved base coordinates
    
    %% maketform + imtransform (original)
    mytform = maketform('affine',input_points(1:3,:), base_points(1:3,:));
    reg2ref = imtransform(double(reg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    reg2ref_dfof = imtransform(double(data_dfof_max_reg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    [rgb_reg2ref, rgb_reg2ref_dfof, yellowTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rOriginal.png'); %save figure 1
    saveas(2, 'dfofOriginal.png'); %save figure 2
    
    %% fitgeotrans affine
    fitGeoTForm = fitgeotrans(input_points(:,:), base_points(:,:),'affine'); %creating transformation
    reg2ref = imwarp(double(reg),fitGeoTForm, 'OutputView', imref2d(size(reg))); %applying transformation on reg
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),fitGeoTForm, 'OutputView', imref2d(size(reg))); %applying transformation on data_dfof_max_reg
    [rgb_reg2ref, rgb_reg2ref_dfof, yellowTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rFitGeoAffine.png'); %saving figure 1
    saveas(2, 'dfofFitGeoAffine.png'); %saving figure 2
    
    %% fitgeotrans nonreflective similarity
    fitGeoTForm = fitgeotrans(input_points(:,:), base_points(:,:),'nonreflectivesimilarity'); %creating transformation
    reg2ref = imwarp(double(reg),fitGeoTForm, 'OutputView', imref2d(size(reg))); %applying transformation on reg
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),fitGeoTForm, 'OutputView', imref2d(size(reg))); %applying transformation on data_dfof_max_reg
    [rgb_reg2ref, rgb_reg2ref_dfof, yellowTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target); %performing create RBG section below (line 114)
    saveas(1, 'r2rFitGeoNonReflSim.png'); %saving figure 1
    saveas(2, 'dfofFitGeoNonReflSim.png'); %saving figure 2
   
    %% imregcorr similarity (nonreflective)
    figure(1); clf; plot(input_points(:,1), input_points(:,2), '.'); %creating plot of coordinates of input_points
    figure(2); clf; plot(base_points(:,1), base_points(:,2), '.'); %creating plot of coordinates of base_points
    saveas(1, 'inputPlot.png'); %saving plot of input_points
    saveas(2, 'basePlot.png'); %saving plot of base_points
    regIm = imread('inputPlot.png'); %loading image of input_points
    refIm = imread('basePlot.png'); %loading image of base_points
    imRegTForm = imregcorr(regIm, refIm,'similarity'); %creating transformation with imregcorr on images of input_points and base_points
    reg2ref = imwarp(double(reg),imRegTForm, 'OutputView', imref2d(size(reg))); %applying transformation on reg
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegTForm, 'OutputView', imref2d(size(reg))); %applying transformation on data_dfof_max_reg
    [rgb_reg2ref, rgb_reg2ref_dfof, yellowTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(3, 'r2rImregcorrSim.png'); %saving figure 3
    saveas(4, 'dfofImregcorrSim.png'); %saving figure 4
        
    %% imregtform Similarity
    [optimizer, metric] = imregconfig('Monomodal'); %creating optimizer and metric for transformation
    imRegisTForm = imregtform(double(reg), double(ref), 'similarity', optimizer, metric); %creating transformation with imregtform on images of input_points and base_points
    reg2ref = imwarp(double(reg),imRegisTForm, 'OutputView', imref2d(size(reg))); %applying transformation on reg
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegisTForm, 'OutputView', imref2d(size(reg))); %applying transformation on data_dfof_max_reg
    [rgb_reg2ref, rgb_reg2ref_dfof, yellowTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegisSim.png'); %saving figure 1
    saveas(2, 'dfofImRegisSim.png'); %saving figure 2
    
    %% imregdemons
    [D, reg2ref] = imregdemons(double(reg), double(ref),300); %D is the displacement field generated, reg2ref is the transformed image
    %[D2, reg2ref_dfof] = imregdemons(double(data_dfof_max_reg), double(data_dfof_max_ref)); % D2 is displacement field for data_dfof_max_ref, reg2ref_dfof is transformed image
    reg2ref_dfof = imwarp(double(data_dfof_max_reg), D); %transforming data_dfof_max_reg with the displacement field generated through reg and ref
    [rgb_reg2ref, rgb_reg2ref_dfof, yellowTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImregdemons.png'); %saving figure 1
    saveas(2, 'dfofImregdemons.png'); %saving figure 2
   
%% Creating RGB Images, this section is in the function createRGB
%     Red/green images created here. rgb_ref2ref is using data_reg_avg, and rgb_reg2ref_dfof is using
%     data_dfof_max
    rgb_reg2ref = zeros(sz_target(1), sz_target(2), 3);
    rgb_reg2ref_dfof = zeros(sz_target(1), sz_target(2), 3);
    rgb_reg2ref(:,:,1) = ref;
    rgb_reg2ref(:,:,2) = reg2ref;
    rgb_reg2ref_dfof(:,:,1) = data_dfof_max_ref;
    rgb_reg2ref_dfof(:,:,2) = reg2ref_dfof;
    figure; imagesc(rgb_reg2ref); title(['day 2 on day 1, data reg rgb overlay'])
%     mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays']))
%     print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_registered_reg2ref_overlay.pdf']), '-dpdf','-bestfit')
    figure; imagesc(rgb_reg2ref_dfof); title(['day 2 on day 1, dfof rgb overlay'])
%     print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_registered_dfof_overlay.pdf']), '-dpdf','-bestfit')

%% Plotting
%     the first figure I showed
    figure; 
    subplot(2,3,1); imagesc(ref); axis off; axis equal; title('day 1 data avg')
    subplot(2,3,2); imagesc(reg); axis off; axis equal; title('day 2 data avg')
    subplot(2,3,3); imagesc(reg2ref); axis off; axis equal; title('registered day 2 data avg')
    subplot(2,3,4); imagesc(rgb_reg2ref_dfof(:,:,1)); axis off; axis equal; title('day 1 data dfof')
    subplot(2,3,5); imagesc(data_dfof_max_reg); axis off; axis equal; title('day 2 data dfof')
    subplot(2,3,6); imagesc(rgb_reg2ref_dfof(:,:,2)); axis off; axis equal; title('registered day 2 data dfof')
%     print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_rgb_individuals.pdf']), '-dpdf','-bestfit')
