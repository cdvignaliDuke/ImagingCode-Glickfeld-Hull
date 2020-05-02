%% get path names
% date is day2
date = '200110'; 
ImgFolder = strvcat('003');
time = strvcat('1409');
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
    sz_target  = size(reg);
    figure(1);clf;imshowpair(reg, ref, 'montage');
%     BW1 = edge(ref,'Canny', [0.01 0.15]); BW2 = edge(reg,'Canny', [0.01 0.1]);
%     figure(2);clf;imshowpair(BW1, BW2, 'montage');
%     refTL = ref(1:512/2, 1:796/2);
%     refTR = ref(1:512/2, 796/2:end);
%     refBL = ref(512/2:end, 1:796/2);
%     refBR = ref(512/2:end, 796/2:end);

%% Transforming
%     selecting the matching landmarks on day 1 and day 2 data_reg_avg and
%     shifting the position of day 2 data_reg_avg and data_dfof_max accordingly
    %% Select Points
    sz_target  = size(reg);
    [input_points, base_points] = cpselect(double(reg),double(ref),'Wait', true);
    save('day4Inputs.mat','input_points');
    save('day4Base.mat','base_points');
    
    %% maketform + imtransform (original)
    mytform = maketform('affine',input_points(1:3,:), base_points(1:3,:));
    reg2ref = imtransform(double(reg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    reg2ref_dfof = imtransform(double(data_dfof_max_reg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rOriginal.png');
    saveas(2, 'dfofOriginal.png');
    
    %% fitgeotrans affine
    fitGeoTForm = fitgeotrans(input_points(:,:), base_points(:,:),'affine');
    %fitGeoOutput = affineOutputView(size(reg), fitGeoTForm);
    reg2ref = imwarp(double(reg),fitGeoTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),fitGeoTForm, 'OutputView', imref2d(size(reg)));
    sz_target = size(ref);
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rFitGeoAffine.png');
    saveas(2, 'dfofFitGeoAffine.png');
    
    %% fitgeotrans nonreflective similarity
    fitGeoTForm = fitgeotrans(input_points(:,:), base_points(:,:),'nonreflectivesimilarity');
    %fitGeoOutput = affineOutputView(size(reg), fitGeoTForm);
    reg2ref = imwarp(double(reg),fitGeoTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),fitGeoTForm, 'OutputView', imref2d(size(reg)));
    sz_target = size(ref);
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rFitGeoNonReflSim.png');
    saveas(2, 'dfofFitGeoNonReflSim.png');
    
    %% fitgeotrans piecewise linear
    fitGeoTForm = fitgeotrans(input_points(:,:), base_points(:,:),'pwl');
    %fitGeoOutput = affineOutputView(size(reg), fitGeoTForm);
    reg2ref = imwarp(double(reg),fitGeoTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),fitGeoTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rFitGeoPWL.png');
    saveas(2, 'dfofFitGeoPWL.png');

    %% fitgeotrans local weighted mean
    fitGeoTForm = fitgeotrans(input_points(:,:), base_points(:,:),'lwm',10);
    %fitGeoOutput = affineOutputView(size(reg), fitGeoTForm);
    reg2ref = imwarp(double(reg),fitGeoTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),fitGeoTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rFitGeoLWM.png');
    saveas(2, 'dfofFitGeoLWM.png');
    
    %% fitgeotrans projective
    fitGeoTForm = fitgeotrans(input_points(:,:), base_points(:,:),'projective');
    %fitGeoOutput = affineOutputView(size(reg), fitGeoTForm);
    reg2ref = imwarp(double(reg),fitGeoTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),fitGeoTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rFitGeoProj.png');
    saveas(2, 'dfofFitGeoProj.png');
   
    %% imregcorr similarity (nonreflective)
    imRegTForm = imregcorr(input_points(:,:), base_points(:,:),'similarity');
    reg2ref = imwarp(double(reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegCorrNonReflSim.png');
    saveas(2, 'dfofImRegCorrNonReflSim.png');

    %% imregcorr rigid
    imRegTForm = imregcorr(input_points(:,:), base_points(:,:),'rigid');
    reg2ref = imwarp(double(reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegCorrRigid.png');
    saveas(2, 'dfofImRegCorrRigid.png');
    
    %% imregcorr translation 
    imRegTForm = imregcorr(input_points(:,:), base_points(:,:),'translation');
    reg2ref = imwarp(double(reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegCorrTransl.png');
    saveas(2, 'dfofImRegCorrTransl.png');
    
    %% imregister translation (intensity based image registration)
    [optimizer, metric] = imregconfig('Monomodal');
    reg2ref = imregister(double(reg), double(ref), 'translation', optimizer, metric);
    reg2ref_dfof = imregister(double(data_dfof_max_reg), double(data_dfof_max_ref), 'translation', optimizer, metric);
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegisterTransl.png');
    saveas(2, 'dfofImRegisterTransl.png');

    %% imregister rigid (intensity based image registration)
    [optimizer, metric] = imregconfig('Monomodal');
    reg2ref = imregister(double(reg), double(ref), 'rigid', optimizer, metric);
    reg2ref_dfof = imregister(double(data_dfof_max_reg), double(data_dfof_max_ref), 'rigid', optimizer, metric);
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegisterRigid.png');
    saveas(2, 'dfofImRegisterRigid.png');

    %% imregister similarity (intensity based image registration)
    [optimizer, metric] = imregconfig('Monomodal');
    reg2ref = imregister(double(reg), double(ref), 'similarity', optimizer, metric);
    reg2ref_dfof = imregister(double(data_dfof_max_reg), double(data_dfof_max_ref), 'similarity', optimizer, metric);
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegisterSimilarity.png');
    saveas(2, 'dfofImRegisterSimilarity.png');

    %% imregister affine (intensity based image registration)
    [optimizer, metric] = imregconfig('Monomodal');
    reg2ref = imregister(double(reg), double(ref), 'affine', optimizer, metric);
    reg2ref_dfof = imregister(double(data_dfof_max_reg), double(data_dfof_max_ref), 'affine', optimizer, metric);
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegisterAffine.png');
    saveas(2, 'dfofImRegisterAffine.png');
    
    %% imregtform translation
    [optimizer, metric] = imregconfig('Monomodal');
    imRegisTForm = imregtform(double(reg), double(ref), 'translation', optimizer, metric);
    reg2ref = imwarp(double(reg),imRegisTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegisTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegisTransl.png');
    saveas(2, 'dfofImRegisTransl.png');
    
    %% imregtform rigid
    [optimizer, metric] = imregconfig('Monomodal');
    imRegisTForm = imregtform(double(reg), double(ref), 'rigid', optimizer, metric);
    reg2ref = imwarp(double(reg),imRegisTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegisTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegisRigid.png');
    saveas(2, 'dfofImRegisRigid.png');

    %% imregtform Similarity
    [optimizer, metric] = imregconfig('Monomodal');
    imRegisTForm = imregtform(double(reg), double(ref), 'similarity', optimizer, metric);
    reg2ref = imwarp(double(reg),imRegisTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegisTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegisSim.png');
    saveas(2, 'dfofImRegisSim.png');

    %% imregtform affine
    [optimizer, metric] = imregconfig('Monomodal');
    imRegisTForm = imregtform(double(reg), double(ref), 'affine', optimizer, metric);
    reg2ref = imwarp(double(reg),imRegisTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegisTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegisAffine.png');
    saveas(2, 'dfofImRegisAffine.png');

    %% imregmtb
    [reg2ref, regShift] = imregmtb(double(reg), double(ref));
    [reg2ref_dfof, dfofShift] = imregmtb(double(data_dfof_max_reg), double(data_dfof_max_ref));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImregmtb.png');
    saveas(2, 'dfofImregmtb.png');
%     saveas(1, 'r2rImregmtb1500.png');
%     saveas(2, 'dfofImregmtb1500.png');
    
    %% imregdemons
    [~, reg2ref] = imregdemons(double(reg), double(ref),300);
    [~, reg2ref_dfof] = imregdemons(double(data_dfof_max_reg), double(data_dfof_max_ref));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
%     saveas(1, 'r2rImregdemons.png');
%     saveas(2, 'dfofImregdemons.png');    
    saveas(1, 'r2rImregdemons1500.png');
    saveas(2, 'dfofImregdemons1500.png');     
    
    
    
   
    
%% Images   
    %% imregcorr similarity (nonreflective) image
    imRegTForm = imregcorr(double(reg), double(ref),'similarity');
    reg2ref = imwarp(double(reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegCorrNonReflSimImage.png');
    saveas(2, 'dfofImRegCorrNonReflSimImage.png');

    %% imregcorr rigid image
    imRegTForm = imregcorr(double(reg), double(ref),'rigid');
    reg2ref = imwarp(double(reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegCorrRigidImage.png');
    saveas(2, 'dfofImRegCorrRigidImage.png');
    
    %% imregcorr translation image
    imRegTForm = imregcorr(double(reg), double(ref),'translation');
    reg2ref = imwarp(double(reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(1, 'r2rImRegCorrTranslImage.png');
    saveas(2, 'dfofImRegCorrTranslImage.png');    
    
%% imregcorr testing    
    %% imregcorr similarity (nonreflective) imagePlot 1
    figure(1); clf; plot(input_points, '.');
    figure(2); clf; plot(base_points, '.');
    saveas(1, 'inputPlot.png');
    saveas(2, 'basePlot.png');
    regIm = imread('inputPlot.png');
    refIm = imread('basePlot.png');
    imRegTForm = imregcorr(regIm, refIm,'similarity');
    reg2ref = imwarp(double(reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(3, 'r2rImRegCorrNonReflSimPointsImage1.png');
    saveas(4, 'dfofImRegCorrNonReflSimPointsImage1.png');    
    
    
    %% imregcorr similarity (nonreflective) imagePlot 2
    figure(1); clf; plot(input_points(:,1), input_points(:,2), 'o');
    %set(gca, 'ydir', 'reverse');
    figure(2); clf; plot(base_points(:,1), base_points(:,2), 'o');
    %set(gca, 'ydir', 'reverse');
    saveas(1, 'inputPlot.png');
    saveas(2, 'basePlot.png');
    regIm = imread('inputPlot.png');
    refIm = imread('basePlot.png');
    imRegTForm = imregcorr(regIm, refIm,'similarity');
    reg2ref = imwarp(double(reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    reg2ref_dfof = imwarp(double(data_dfof_max_reg),imRegTForm, 'OutputView', imref2d(size(reg)));
    [rgb_reg2ref, rgb_reg2ref_dfof, greenTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target);
    saveas(3, 'r2rImRegCorrNonReflSimPointsImage2.png');
    saveas(4, 'dfofImRegCorrNonReflSimPointsImage2.png');    
    
%% Creating RGB Images
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
    greenpixels = rgb_reg2ref_dfof(:,:,1) < 0.12 & rgb_reg2ref_dfof(:,:,2) >= 0.24 & rgb_reg2ref_dfof(:,:,3) < 0.12; sum(greenpixels(:))

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

%% Saving
%    each frame within day 2 data_reg is shifted according to how the
%    landmarks match up
    for i = 1:nframes
        data_reg(:,:,i) = imtransform(double(data_reg(:,:,i)),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
        if rem(i,50) == 0
            fprintf([num2str(i) '/n'])
        end
    end

    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_transform.mat']), 'input_points', 'base_points', 'mytform', 'data_dfof_max_ref', 'ref', 'reg2ref', 'reg2ref_dfof');
    load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'mask_cell', 'mask_np')

%% neuropil subtraction
data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_tc_down = stackGetTimeCourses(stackGroupProject(data_reg,5), mask_cell);
nCells = size(data_tc,2);
%np_tc = stackGetTimeCourses(data_reg,mask_np);
clear np_tc np_tc_down
sz = size(data_reg);
down = 5;
data_reg_down  = stackGroupProject(data_reg,down);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s/n']) 
end
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
