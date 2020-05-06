%% align 2p images to WF
%% load csv to define exp

clear all; clc;
%expfile = '\\CRASH.dhe.duke.edu\data\home\kevin\Code\Ai9x_experiment_list.txt';
expfile = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Code\Ai9x_experiment_list.txt';
fID = fopen(expfile);
head = textscan(fID,'%s%s%s%s%s',1,'delimiter',',');
head = vertcat(head{:});
temp = textscan(fID,'%s%s%s%s%s','delimiter',',','HeaderLines',1);
temp = horzcat(temp{:});
expdata = cell2table(temp,'VariableNames',head);
nExp = size(expdata,1);
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);

fprintf(['Receptive-field-size visual-area comparison analysis - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

%% load 2P images
mean_imgs = cell(4,1);
figure(1);clf;
for iExp = 1:4
    data = [];
    clear info;
    date = expdata.date{iExp};
    mouse = expdata.mouse{iExp};
    run_str = expdata.run_str{iExp};
    
    CD = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Data\2P\' date '_' mouse '\003' ];
    cd(CD);
    
    imgMatFile = '003_000_000.mat';
    load(fullfile(CD,imgMatFile));
    
    %nframes = info.config.frames;
    nframes = 1000;
    fprintf('\nReading exp %i/%i, run 003 - %i frames \n',iExp,nExp,nframes)
    fprintf('Mouse: %s, date: %s\n',mouse,date)
    tic
    %data = sbxread('003_000_000',0,nframes);
    data = sbxread(fullfile(CD,'003_000_000'),0,nframes);
    toc
    data = squeeze(data); % squeeze because only 1 pmt channel
    
    fprintf('Loaded data, saving mean image\n')
    mean_imgs{iExp} = mean(data,3);
    
    figure(1)
    subplot(2,2,iExp)
    imagesc(mean_imgs{iExp})
    title(['Exp ' num2str(iExp) ' - Mouse: ' mouse ', date: ' date])
    axis image
end
%% load WF image
switch mouse
    case 'i838'
        WFroot = 'N:\home\kevin\Data\Widefield\KM38_180305';
        WFimage = 'KM38_180305_45deg_el[15 0 -15]_az[-20 0 20]_FOV.tif';
    case 'i840'
        WFroot = 'N:\home\kevin\Data\Widefield\KM40_180305';
        WFimage = 'KM40_180305_45deg_el[15 0 -15]_az[-20 0 20]_FOV.tif';
    case 'i842'
        WFroot = 'N:\home\kevin\Data\Widefield\KM42_180301';
        WFimage = 'KM42_180301_45deg_el[15 0 -15]_az[-20 0 20]_FOV.tif';
    case 'i843'
        WFroot = 'N:\home\kevin\Data\Widefield\KM43_180302';
        WFimage = 'KM43_180302_45deg_el[15 0 -15]_az[-20 0 20]_FOV.tif';
    case 'i844'
        WFroot = 'N:\home\kevin\Data\Widefield\KM44_180302';
        WFimage = 'KM44_180302_45deg_el[15 0 -15]_az[-20 0 20]_FOV.tif';
    case 'i880'
        WFroot = 'N:\home\kevin\Data\Widefield\i880_180621';
        WFimage = 'i880_180621_45deg_el[-15 0 15]_az[-20 0 20]_FOV.tif';
    case 'i881'
        WFroot = 'N:\home\kevin\Data\Widefield\i881_180621';
        WFimage = 'i881_180621_45deg_el[-15 0 15]_az[-20 0 20]_FOV.tif';
    case 'i883'
        WFroot = 'N:\home\kevin\Data\Widefield\i883_180706';
        WFimage = 'i883_180706_45deg_el[15 0 -15]_az[-20 0 20]_FOV.tif';
    case 'i884'
        WFroot = 'N:\home\kevin\Data\Widefield\i884_180705';
        WFimage = 'i884_180705_45deg_el[15 0 -15]_az[-20 0 20]_FOV.tif';
end
fixed = double(imread(fullfile(WFroot,WFimage)));
for i=1:3
    azim(:,:,i) = imread(fullfile(WFroot,'Composite_az.tif'),i);
    temp = imread(fullfile(WFroot,'Composite_az.tif'),i);
    %temp = temp-min(temp(:)); temp = temp/max(temp(:))*255;
    %azim(:,:,i) = single(temp);
    elim(:,:,i) = imread(fullfile(WFroot,'Composite_el.tif'),i);
    temp = imread(fullfile(WFroot,'Composite_el.tif'),i);
    %temp = temp-min(temp(:)); temp = temp/max(temp(:))*255;
    %elim(:,:,i) = single(temp);
end
figure(2);clf;
subplot(1,3,1)
imagesc(fixed)
subplot(1,3,2)
imagesc(azim)
subplot(1,3,3)
imagesc(elim)

%% show WF and 2P in one plot
figure(1);clf;
subplot(4,2,[1:4])
imagesc(fixed)
title(['WF FOV - Mouse: ' mouse])
axis image;colormap gray
subplot(4,2,5)
imagesc(mean_imgs{1})
%title({['Exp ' num2str(1) ' - Mouse: ' mouse ', date: ' cell2mat(expdata.date(1))],'Area=V1'})
title(['Exp ' num2str(1) ' Area=V1'])
axis image
subplot(4,2,6)
imagesc(mean_imgs{3})
%title({['Exp ' num2str(3) ' - Mouse: ' mouse ', date: ' cell2mat(expdata.date(3))],'Area=LM'})
title(['Exp ' num2str(3) ' Area=LM'])
axis image
subplot(4,2,7)
imagesc(mean_imgs{2})
%title({['Exp ' num2str(2) ' - Mouse: ' mouse ', date: ' cell2mat(expdata.date(2))],'Area=AL'})
title(['Exp ' num2str(2) ' Area=AL'])
axis image
subplot(4,2,8)
gamma=0.3;
%imagesc(mean_imgs{4})
imagesc((mean_imgs{4}/255).^gamma*255) % gamma corrected for brighter
%title({['Exp ' num2str(4) ' - Mouse: ' mouse ', date: ' cell2mat(expdata.date(4))],'Area=PM'})
title(['Exp ' num2str(4) ' Area=PM'])
axis image
%% align by manual values

fixedim = fixed-min(fixed(:));
fixedim = fixedim/max(fixedim(:))*255;
fixedim_resize=imresize(fixedim,17); % can just do this because scale same for all 4 exps
moving_overlay = 0*fixedim_resize;

figure(2);clf;
imagesc(fixed)
for iImg = 1:4
    %moving = mean_imgs{iImg};
    %moving = rot90(mean_imgs{iImg});
    moving = mean_imgs{iImg}';
    
    % define transformation parameters for each exp (manual)
    switch iImg
        case 1
            p = [139 220 58 17];
            sp = 1;
        case 2
            p = [84 173 58 17];
            sp = 3;
        case 3
            p = [64 217 58 17];
            sp = 2;
        case 4
            p = [198 167 58 17];
            sp = 4;
    end
    Xpos=p(1); Ypos=p(2); rotation=p(3); scale=p(4);
    
    % map transform and add aligned image (2P) into fixedim_resize (WF)
    xOff = Ypos*scale-size(moving,2)/2;
    yOff = Xpos*scale-size(moving,1)/2;
    [xx, yy] = ndgrid(1:size(moving_overlay,1),1:size(moving_overlay,2));
    mapX = round((xx-xOff)*cosd(rotation) + (yy-yOff)*sind(rotation) +size(moving,2)/2);
    mapY = round((yy-yOff)*cosd(rotation) - (xx-xOff)*sind(rotation) +size(moving,1)/2);
    zeromask = mapX<1 | mapX>size(moving,2) | mapY<1 | mapY>size(moving,1);
    mapX(zeromask) = 0; mapY(zeromask)=0;
    useind = find(~zeromask);
    mapInd = sub2ind(size(moving),mapY(useind),mapX(useind));
    moving_overlay(useind) = moving(mapInd);
    moving_overlay = 255*moving_overlay/max(moving_overlay(:));
    
    overlay_temp(:,:,iImg) = moving_overlay;
end

%% plot overlays for supp fig
subplot(4,4,1)
imagesc(mean_imgs{1})
%title({['Exp ' num2str(1) ' - Mouse: ' mouse ', date: ' cell2mat(expdata.date(1))],'Area=V1'})
title(['Exp ' num2str(1) ' Area=V1'])
axis image;set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);set(gca,'XTick',[]);set(gca,'YTick',[])
hold on
plot([100 230],[470 470],'w-')
text(165,420,'100 um','Color','w','HorizontalAlignment','center')
hold off
subplot(4,4,2)
imagesc(mean_imgs{3})
%title({['Exp ' num2str(3) ' - Mouse: ' mouse ', date: ' cell2mat(expdata.date(3))],'Area=LM'})
title(['Exp ' num2str(3) ' Area=LM'])
axis image;set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);set(gca,'XTick',[]);set(gca,'YTick',[])
subplot(4,4,5)
imagesc(mean_imgs{2})
%title({['Exp ' num2str(2) ' - Mouse: ' mouse ', date: ' cell2mat(expdata.date(2))],'Area=AL'})
title(['Exp ' num2str(2) ' Area=AL'])
axis image;set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);set(gca,'XTick',[]);set(gca,'YTick',[])
subplot(4,4,6)
gamma=0.3;
%imagesc(mean_imgs{4})
imagesc((mean_imgs{4}/255).^gamma*255) % gamma corrected for brighter
%title({['Exp ' num2str(4) ' - Mouse: ' mouse ', date: ' cell2mat(expdata.date(4))],'Area=PM'})
title(['Exp ' num2str(4) ' Area=PM'])
axis image;set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);set(gca,'XTick',[]);set(gca,'YTick',[])

moving_overlay = sum(overlay_temp,3);
border_overlay = (conv2(bwmorph(moving_overlay>0,'remove'),ones(10),'same')>0);
meanIm = (fixedim_resize + moving_overlay)/2;
subplot(4,4,[3 4 7 8])
imagesc(meanIm(1001:4250,201:3450));
title('WF FOV')
axis image;set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);set(gca,'XTick',[]);set(gca,'YTick',[])
hold on
plot([500 500+1420]-200,[4100 4100]-1000,'w-')
text(1210-200,3900-1000,'1 mm','Color','w','HorizontalAlignment','center')
hold off
subplot(4,4,[9 10 13 14])
azim_resize=imresize(azim,17);
meanIm = azim_resize + single(border_overlay)/2;
imagesc(meanIm(1001:4250,201:3450,:));
title('Az Map')
axis image;set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);set(gca,'XTick',[]);set(gca,'YTick',[])
subplot(4,4,[11 12 15 16])
elim_resize=imresize(elim,17);
gamma = [1.1 1.05 1];
for i=1:3
    elim_resize(:,:,i) = elim_resize(:,:,i)*gamma(i);
end
meanIm = elim_resize + single(border_overlay)/2;
imagesc(meanIm(1001:4250,201:3450,:));
title('El Map')
axis image;set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);set(gca,'XTick',[]);set(gca,'YTick',[])

% set(gcf,'PaperOrientation','landscape');
% set(gcf,'PaperUnits','normalized');
% set(gcf,'PaperPosition', [0 0 1 1]);
% filename = fullfile('N:\home\kevin\ppts\_paper figs','i838_1-4_2PtoWF_aligned.pdf');
% print(filename,'-dpdf')
%% click a guess
w=0
hold on
tempTarg=plot(size(fixed)/2,'gx');
while ~w
    w = waitforbuttonpress
    tempTarg.Visible='off';
    click = get(gca,'CurrentPoint');
    tempTarg = plot(click(1),click(3),'rx');
end
clicky = click(1); clickx = click(3);
%% automatic alignment doesnt work
fprintf('\nAligning...')
[optimizer,metric] = imregconfig('multimodal');
optimizer.InitialRadius = 0.009;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;
fixed_shad = imgaussfilt(fixed,1);
R_moving = imref2d(size(moving),0.7,0.7);
R_fixed = imref2d(size(fixed),16,16);
tx = clicky*16-size(moving,2)*0.7/2;
ty = clickx*16-size(moving,1)*0.7/2;
iTform = affine2d([1 0 0; 0 1 0; tx ty 1])
[moving_reg, R_reg] = imregister(moving, R_moving, fixed_shad, R_fixed, 'similarity', optimizer, metric, 'InitialTransformation', iTform);
fprintf('done\n')

figure(2);clf;
imshowpair(moving_reg,fixed_shad,'blend')

%% try with normxcorr2
% find max but also incrementally rotate fixed
rot = 20:5:30;
fixed_shad_resize = imresize(fixed_shad(150:250,50:150),20);
for iTh = 1:length(rot)
    temp = imrotate(fixed_shad_resize,rot(iTh),'bilinear','crop');
    c = normxcorr2(moving,temp);
    maxcorr = max(c(:));
    [ypeak, xpeak] = find(c==max(c(:)));
    figure(2);clf
    %imagesc(temp)
    temp2 = 0*temp;
    [r, c] = ind2sub(size(moving),find(moving));
    ind = sub2ind(size(temp),ypeak+r,xpeak+c);
    temp2(ind) = moving(:);
    subplot(2,3,1)
    imagesc(moving); colormap gray
    subplot(2,3,4)
    imagesc(temp)
    subplot(2,3,[2 3 5 6])
    imshowpair(temp,temp2,'blend');
    hold on
    plot(xpeak,ypeak,'rx')
    hold off
    title(['Rotation: ' num2str(rot(iTh))])
    pause
end

%% gabor filters for fig schematics
figure(3);clf
x=-500:500;
y=x;
[xx,yy]=ndgrid(x,y);
rr = sqrt(xx.^2+yy.^2);
gg = sin(0.08*yy).*(rr<139.5); %62 93 139.5 209.25
imagesc(gg*0.4)
colormap gray
axis square; axis off
clim([-1 1])
