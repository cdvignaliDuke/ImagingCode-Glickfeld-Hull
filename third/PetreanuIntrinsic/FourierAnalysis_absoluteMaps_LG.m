close all
clear all

%%
animalName = 'i1303';
date = '191115';
time = '0956';
saveName = 'i1303';
cond = {'VR', 'HD', 'VL', 'HU'};
stim_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
database_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\Widefield_images';
fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Widefield_imaging';

fprintf(['Create ' animalName ' ' date ' absolute maps \n' ])
%% stimulus data

stimdata_fn = ['data-' animalName '-' date '-' time '.mat'];
load(fullfile(stim_fn,stimdata_fn))

% parameters
scalePxPerMM=62.5; %example for our average zoom

elevationBorders = double([-input.monitorYDeg input.monitorYDeg]);
azimuthBorders = double([-input.monitorXDeg input.monitorXDeg]);
azimuthOffset=0;
elevationOffset=0;
save_images=1;
save_maps=1;

fCamera = double(input.frameRateHz);
cITIStart = celleqel2mat_padded(input.cITIStart);
cStimOn = celleqel2mat_padded(input.cStimOn);

roiDiskRadius_um=50;
radius=roiDiskRadius_um/(1000/scalePxPerMM);

colormapRes = 128;
anglesMap = jet(colormapRes);

fprintf(['Loaded stimulus data \n' ])
%% Import section

artificial_delay=pi/3;
% filt_size =50;
% filt_sigma = 1;
%h = fspecial('gaussian', filt_size, filt_sigma);
h=fspecial('disk',radius);

% Load files - output of FourierAnalysis_relativeMaps

vessels = load(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_relative_' cond{1} '.mat']),'vessels'); vessels = rot90(vessels.vessels,-1);

phaseRight = load(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_relative_' cond{1} '.mat']),'phase');
phaseRight = fliplr(phaseRight.phase);
phaseLeft = load(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_relative_' cond{3} '.mat']),'phase');
phaseLeft = fliplr(phaseLeft.phase);
phaseDown = load(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_relative_' cond{2} '.mat']),'phase');
phaseDown = fliplr(phaseDown.phase);
phaseUp = load(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_relative_' cond{4} '.mat']),'phase');
phaseUp = fliplr(phaseUp.phase);

amplitudeRight = load(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_relative_' cond{1} '.mat']),'amplitude');
amplitudeRight = fliplr(amplitudeRight.amplitude);
amplitudeLeft = load(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_relative_' cond{3} '.mat']),'amplitude');
amplitudeLeft = fliplr(amplitudeLeft.amplitude);
amplitudeDown = load(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_relative_' cond{2} '.mat']),'amplitude');
amplitudeDown = fliplr(amplitudeDown.amplitude);
amplitudeUp = load(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_relative_' cond{4} '.mat']),'amplitude');
amplitudeUp = fliplr(amplitudeUp.amplitude);

x_px=1:size(phaseRight,2);
y_px=1:size(phaseRight,1);
[X_mat, Y_mat]=meshgrid(x_px,y_px);
phaseRightFilt=zeros(size(phaseRight));
phaseLeftFilt=zeros(size(phaseRight));
phaseDownFilt=zeros(size(phaseRight));
phaseUpFilt=zeros(size(phaseRight));

for i=1:size(phaseRight,1)
    for j=1:size(phaseRight,2)
        mask=(Y_mat-i).^2+(X_mat-j).^2<=radius^2;
        
        phaseR_values=phaseRight(mask);
        phaseL_values=phaseLeft(mask);
        phaseRightFilt(i,j)=circ_mean(phaseR_values);
        phaseLeftFilt(i,j)=circ_mean(phaseL_values);
        
        phaseD_values=phaseDown(mask);
        phaseU_values=phaseUp(mask);
        phaseDownFilt(i,j)=circ_mean(phaseD_values);
        phaseUpFilt(i,j)=circ_mean(phaseU_values);
        
    end
end

cycDur = cITIStart(2)-cStimOn(1);
framesPerBin = cycDur+fCamera;
stimPeriod = framesPerBin./fCamera;
T_right = stimPeriod;
cycDur = cITIStart(3)-cStimOn(2);
framesPerBin = cycDur+fCamera;
stimPeriod = framesPerBin./fCamera;
T_down = stimPeriod;
azimuthRange = double(input.monitorXDeg.*2);
elevationRange = double(input.monitorYDeg.*2);

fprintf(['Imported image data \n' ])

%% correct for delays and make phase and sign maps

% Introduce delay
phaseRight=phaseRight+artificial_delay;
phaseLeft=phaseLeft+artificial_delay;
phaseDown=phaseDown+artificial_delay;
phaseUp=phaseUp+artificial_delay;

phaseRightFilt=phaseRightFilt+artificial_delay;
phaseLeftFilt=phaseLeftFilt+artificial_delay;
phaseDownFilt=phaseDownFilt+artificial_delay;
phaseUpFilt=phaseUpFilt+artificial_delay;

fprintf(['Introduced delay \n' ])

% Unwrap
phaseRight(phaseRight>pi)=phaseRight(phaseRight>pi)-2*pi;
phaseLeft(phaseLeft>pi)=phaseLeft(phaseLeft>pi)-2*pi;
phaseDown(phaseDown>pi)=phaseDown(phaseDown>pi)-2*pi;
phaseUp(phaseUp>pi)=phaseUp(phaseUp>pi)-2*pi;

phaseRightFilt(phaseRightFilt>pi)=phaseRightFilt(phaseRightFilt>pi)-2*pi;
phaseLeftFilt(phaseLeftFilt>pi)=phaseLeftFilt(phaseLeftFilt>pi)-2*pi;
phaseDownFilt(phaseDownFilt>pi)=phaseDownFilt(phaseDownFilt>pi)-2*pi;
phaseUpFilt(phaseUpFilt>pi)=phaseUpFilt(phaseUpFilt>pi)-2*pi;

fprintf(['Unwrapped \n' ])

% Calculates angles
phaseAzimuth = (phaseLeft-phaseRight)/2;
phaseElevation = (phaseDown-phaseUp)/2;
angleAzimuth=(phaseAzimuth-(-pi))/(2*pi)*azimuthRange+azimuthBorders(1);
angleElevation=(phaseElevation-(-pi))/(2*pi)*elevationRange+elevationBorders(1);

phaseAzimuthFilt = (phaseLeftFilt-phaseRightFilt)/2;
phaseElevationFilt = (phaseDownFilt-phaseUpFilt)/2;
angleAzimuthFilt=(phaseAzimuthFilt-(-pi))/(2*pi)*azimuthRange+azimuthBorders(1);
angleElevationFilt=(phaseElevationFilt-(-pi))/(2*pi)*elevationRange+elevationBorders(1);

fprintf(['Calculated angles \n' ])
% Calculates amplitude
amplitudeAzimuth = (amplitudeRight+amplitudeLeft)/2;
amplitudeElevation = (amplitudeDown+amplitudeUp)/2;

amplitudeAzimuthFilt = imfilter(amplitudeAzimuth,h);
amplitudeAzimuthFilt= amplitudeAzimuthFilt/max(amplitudeAzimuthFilt(:));
amplitudeElevationFilt = imfilter(amplitudeElevation,h);
amplitudeElevationFilt= amplitudeElevationFilt/max(amplitudeElevationFilt(:));

fprintf(['Calculated amplitude \n' ])
% Calculates delay
phaseAzimuthDelay = (phaseRight+phaseLeft)/2-2*artificial_delay;
phaseElevationDelay = (phaseDown+phaseUp)/2;
delayAzimuth = (phaseAzimuthDelay-(-pi))/(2*pi)*T_right;
delayElevation = (phaseElevationDelay-(-pi))/(2*pi)*T_right;

phaseAzimuthFiltDelay = (phaseRightFilt+phaseLeftFilt)/2-2*artificial_delay;
phaseElevationFiltDelay = (phaseDownFilt+phaseUpFilt)/2;
delayAzimuthFilt = (phaseAzimuthFiltDelay-(-pi))/(2*pi)*T_right;
delayElevationFilt = (phaseElevationFiltDelay-(-pi))/(2*pi)*T_right;

fprintf(['Calculated delay \n' ])

% Calculates sign map - Gabriela 2019-04-08
[dhdx dhdy] = gradient(angleElevationFilt);
[dvdx dvdy] = gradient(angleAzimuthFilt);

graddir_hor = atan2(dhdy,dhdx);
graddir_vert = atan2(dvdy,dvdx);

vdiff = exp(1i*graddir_hor) .* exp(-1i*graddir_vert); %Should be vert-hor, but the gradient in Matlab for y is opposite.
VFS = sin(angle(vdiff)); %Visual field sign map
id = find(isnan(VFS));
VFS(id) = 0;

hh = fspecial('gaussian',size(VFS),3);
hh = hh/sum(hh(:));
VFS = ifft2(fft2(VFS).*abs(fft2(hh)));

fprintf(['Calculated sign map \n' ])
%% Output section
lambda=10;
ampThrs=0.3;

azimuthFigure = figure('Position',[100 100 800 800]);
azimuthHandle = imagesc(angleAzimuth); colormap(anglesMap)
set(azimuthHandle,'AlphaData',(1./(1+exp(-lambda*(amplitudeAzimuthFilt-ampThrs)))));
%set(azimuthHandle,'AlphaData',(1-exp(-amplitudeAzimuthFilt*lambda))/(1-exp(-lambda)));
axis image
colorbar
caxis(azimuthBorders)
title([animalName, ' - Azimuth'])

elevationFigure = figure('Position',[100 100 800 800]);
elevationHandle = imagesc(angleElevation); colormap(anglesMap)
set(elevationHandle,'AlphaData',(1./(1+exp(-lambda*(amplitudeElevationFilt-ampThrs)))));
alphamap('increase',1)
axis image
colorbar
caxis(elevationBorders)
title([animalName, ' - Elevation'])

azimuthFiltFigure=figure('Position',[100 100 800 800]);
azimuthFiltHandle = imagesc(angleAzimuthFilt); colormap(anglesMap)
set(azimuthFiltHandle,'AlphaData',(1./(1+exp(-lambda*(amplitudeAzimuthFilt-ampThrs)))));
axis image
colorbar
caxis(azimuthBorders)
title([animalName, ' - Azimuth'])

elevationFiltFigure=figure('Position',[100 100 800 800]);
elevationFiltHandle = imagesc(angleElevationFilt); colormap(anglesMap)
set(elevationFiltHandle,'AlphaData',(1./(1+exp(-lambda*(amplitudeElevationFilt-ampThrs)))));
alphamap('increase',1)
axis image
colorbar
caxis(elevationBorders)
title([animalName, ' - Elevation'])

vesselsFigure = figure('Position',[100 100 800 800]);
imagesc(vessels); colormap('gray')
axis image
title([animalName, ' - Vessels'])

VFSFigure=figure('Position',[100 100 800 800]);
VFSHandle = imagesc(VFS,[-1 1]); colormap(anglesMap)
alphamap('increase',1)
axis image
colorbar
title([animalName, ' - VSF'])


if save_images
    figure(azimuthFigure)
    saveas(gcf,fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_Azimuth']),'fig');
    %screen2png(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_Azimuth']));
    
    figure(azimuthFiltFigure)
    saveas(gcf,fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_AzimuthFilt']),'fig');
    %screen2png(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_AzimuthFilt']));
    
    figure(elevationFigure)
    saveas(gcf,fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_Elevation']),'fig');
    %screen2png(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_Elevation']));
    
    figure(elevationFiltFigure)
    saveas(gcf,fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_ElevationFilt']),'fig');
    %screen2png(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_ElevationFilt']));
    
    figure(vesselsFigure)
    saveas(gcf,fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_Vessels']),'fig');
    %screen2png(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_Vessels']));
    
    figure(VFSFigure)
    saveas(gcf,fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_VFS']),'fig');
    %screen2png(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_VFS']));
end

if save_maps
    save(fullfile(fn_out, saveName, [date '_' animalName], [date '_' animalName '_maps']), 'angleAzimuth','angleElevation',...
        'angleAzimuthFilt','angleElevationFilt','amplitudeAzimuthFilt','amplitudeElevationFilt','vessels','VFS');
end