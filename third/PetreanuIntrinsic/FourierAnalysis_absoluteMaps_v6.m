
animalName = 'RD10161';
saveName = 'RD10161';

% parameters
scalePxPerMM=62.5; %example for our average zoom

elevationBorders = [-40 40];
azimuthBorders = [-60 60];
azimuthOffset=0;
elevationOffset=0;
save_images=1;
save_maps=1;

roiDiskRadius_um=50;
radius=roiDiskRadius_um/(1000/scalePxPerMM);

colormapRes = 128;
anglesMap = jet(colormapRes);

%%


% parameters


roiDiskRadius_um=50;
radius=roiDiskRadius_um/(1000/scalePxPerMM);

colormapRes = 128;
anglesMap = jet(colormapRes);

%% Import section
artificial_delay=pi/3;
% filt_size =50;
% filt_sigma = 1;
%h = fspecial('gaussian', filt_size, filt_sigma);
h=fspecial('disk',radius);

% Load files - output of FourierAnalysis_relativeMaps
[verticalRightFile pathAnalysis]=uigetfile('.mat','Select Vertical Right File');
cd(pathAnalysis);
[verticalLeftFile]=uigetfile('.mat','Select Vertical Left File');
[horizontalDownFile]=uigetfile('.mat','Select Horizontal Down File');
[horizontalUpFile]=uigetfile('.mat','Select Horizontal Up File');
cd(pathAnalysis);

vessels = load([pathAnalysis verticalRightFile],'vessels'); vessels = rot90(vessels.vessels,-1);

stimDataRight = load([pathAnalysis verticalRightFile],'stimData'); 
stimDataRight = stimDataRight.stimData;
stimDataDown = load([pathAnalysis horizontalDownFile],'stimData'); 
stimDataDown = stimDataDown.stimData;
phaseRight = load([pathAnalysis verticalRightFile],'phase');
phaseRight = fliplr(phaseRight.phase);
phaseLeft = load([pathAnalysis verticalLeftFile],'phase');
phaseLeft = fliplr(phaseLeft.phase);
phaseDown = load([pathAnalysis horizontalDownFile],'phase');
phaseDown = fliplr(phaseDown.phase);
phaseUp = load([pathAnalysis horizontalUpFile],'phase');
phaseUp = fliplr(phaseUp.phase);

amplitudeRight = load([pathAnalysis verticalRightFile],'amplitude');
amplitudeRight = fliplr(amplitudeRight.amplitude);
amplitudeLeft = load([pathAnalysis verticalLeftFile],'amplitude');
amplitudeLeft = fliplr(amplitudeLeft.amplitude);
amplitudeDown = load([pathAnalysis horizontalDownFile],'amplitude');
amplitudeDown = fliplr(amplitudeDown.amplitude);
amplitudeUp = load([pathAnalysis horizontalUpFile],'amplitude');
amplitudeUp = fliplr(amplitudeUp.amplitude);

x_px=1:size(phaseRight,2);
y_px=1:size(phaseRight,1);
[X_mat, Y_mat]=meshgrid(x_px,y_px);
phaseRightFilt=zeros(size(phaseRight));
phaseLeftFilt=zeros(size(phaseRight));
phaseDownFilt=zeros(size(phaseRight));
phaseUpFilt=zeros(size(phaseRight));

for i=1:size(phaseRight,1)
    i
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

T = stimDataRight.stimPeriod;
azimuthRange = stimDataRight.barPos(end)-stimDataRight.barPos(1);
elevationRange = stimDataDown.barPos(end)-stimDataDown.barPos(1);


%% Introduce delay
phaseRight=phaseRight+artificial_delay;
phaseLeft=phaseLeft+artificial_delay;
phaseDown=phaseDown+artificial_delay;
phaseUp=phaseUp+artificial_delay;

phaseRightFilt=phaseRightFilt+artificial_delay;
phaseLeftFilt=phaseLeftFilt+artificial_delay;
phaseDownFilt=phaseDownFilt+artificial_delay;
phaseUpFilt=phaseUpFilt+artificial_delay;

%% Unwrap
phaseRight(phaseRight>pi)=phaseRight(phaseRight>pi)-2*pi;
phaseLeft(phaseLeft>pi)=phaseLeft(phaseLeft>pi)-2*pi;
phaseDown(phaseDown>pi)=phaseDown(phaseDown>pi)-2*pi;
phaseUp(phaseUp>pi)=phaseUp(phaseUp>pi)-2*pi;

phaseRightFilt(phaseRightFilt>pi)=phaseRightFilt(phaseRightFilt>pi)-2*pi;
phaseLeftFilt(phaseLeftFilt>pi)=phaseLeftFilt(phaseLeftFilt>pi)-2*pi;
phaseDownFilt(phaseDownFilt>pi)=phaseDownFilt(phaseDownFilt>pi)-2*pi;
phaseUpFilt(phaseUpFilt>pi)=phaseUpFilt(phaseUpFilt>pi)-2*pi;

%% Calculates angles
phaseAzimuth = (phaseLeft-phaseRight)/2;
phaseElevation = (phaseDown-phaseUp)/2;
angleAzimuth=(phaseAzimuth-(-pi))/(2*pi)*azimuthRange+stimDataRight.barPos(1);
angleElevation=(phaseElevation-(-pi))/(2*pi)*elevationRange+stimDataDown.barPos(1);

phaseAzimuthFilt = (phaseLeftFilt-phaseRightFilt)/2;
phaseElevationFilt = (phaseDownFilt-phaseUpFilt)/2;
angleAzimuthFilt=(phaseAzimuthFilt-(-pi))/(2*pi)*azimuthRange+stimDataRight.barPos(1);
angleElevationFilt=(phaseElevationFilt-(-pi))/(2*pi)*elevationRange+stimDataDown.barPos(1);

%% Calculates amplitude
amplitudeAzimuth = (amplitudeRight+amplitudeLeft)/2;
amplitudeElevation = (amplitudeDown+amplitudeUp)/2;

amplitudeAzimuthFilt = imfilter(amplitudeAzimuth,h);
amplitudeAzimuthFilt= amplitudeAzimuthFilt/max(amplitudeAzimuthFilt(:));
amplitudeElevationFilt = imfilter(amplitudeElevation,h);
amplitudeElevationFilt= amplitudeElevationFilt/max(amplitudeElevationFilt(:));


%% Calculates delay
phaseAzimuthDelay = (phaseRight+phaseLeft)/2-2*artificial_delay;
phaseElevationDelay = (phaseDown+phaseUp)/2;
delayAzimuth = (phaseAzimuthDelay-(-pi))/(2*pi)*T;
delayElevation = (phaseElevationDelay-(-pi))/(2*pi)*T;

phaseAzimuthFiltDelay = (phaseRightFilt+phaseLeftFilt)/2-2*artificial_delay;
phaseElevationFiltDelay = (phaseDownFilt+phaseUpFilt)/2;
delayAzimuthFilt = (phaseAzimuthFiltDelay-(-pi))/(2*pi)*T;
delayElevationFilt = (phaseElevationFiltDelay-(-pi))/(2*pi)*T;

%
% angleAzimuth=imfilter(angleAzimuth,h);
% angleElevation=imfilter(angleElevation,h);

%% Calculates sign map - Gabriela 2019-04-08
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
    folder_name = uigetdir(pathAnalysis,'Choose saving directory');
    figure(azimuthFigure)
    saveas(gcf,[folder_name filesep saveName '_Azimuth'],'fig');
    screen2png([folder_name filesep saveName '_Azimuth']);
    
    figure(azimuthFiltFigure)
    saveas(gcf,[folder_name filesep saveName '_AzimuthFilt'],'fig');
    screen2png([folder_name filesep saveName '_AzimuthFilt']);
    
    figure(elevationFigure)
    saveas(gcf,[folder_name filesep saveName '_Elevation'],'fig');
    screen2png([folder_name filesep saveName '_Elevation']);
    
    figure(elevationFiltFigure)
    saveas(gcf,[folder_name filesep saveName '_ElevationFilt'],'fig');
    screen2png([folder_name filesep saveName '_ElevationFilt']);
    
    figure(vesselsFigure)
    saveas(gcf,[folder_name filesep saveName '_Vessels'],'fig');
    screen2png([folder_name filesep saveName '_Vessels']);
    
    figure(VFSFigure)
    saveas(gcf,[folder_name filesep saveName '_VFS'],'fig');
    screen2png([folder_name filesep saveName '_VFS']);
end

if save_maps
    folder_name = uigetdir(pathAnalysis,'Choose saving directory');
    save([folder_name filesep saveName '_maps'], 'angleAzimuth','angleElevation',...
        'angleAzimuthFilt','angleElevationFilt','amplitudeAzimuthFilt','amplitudeElevationFilt','vessels','VFS');
end