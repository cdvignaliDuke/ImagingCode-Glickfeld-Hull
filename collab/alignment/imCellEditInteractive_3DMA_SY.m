function [bwCell_vol_numbered cell_pos] = imCellEditInteractive_3DMA_SY(fovImg_vol, bwCell_vol_numbered, zoomToAxis, cellRadius, Im_contrast,resolution)
% [bwCell_vol_numbered cell_pos]=imCellEditInteractive_3DMA_SY('E:\users\\Lindsey\091215\invitrostacks\091215_invitro_RGB.tif',[],[],-10,[],[0.27 0.27 2]);
%IMCELLEDITINTERACTIVE Add cells to a cell mask interactively
%  BWOUT = IMCELLEDITINTERACTIVE(FOVIMG, BWCELL, ZOOMTOAXIS, CELLRADIUS)
%
%   Interactively add/del cells to a mask by clicking on them.
%
%   Saves mask at each step to the UserData field of the figure: if this
%   crashes you can pull it out manually.
%   ud = get(gcf, 'UserData');
%   recoveredMask = ud.bwMask;
% previous comment is replaced by making bwCell_vol global
%
%   params:
%
%       fovImg_vol: 3D stack input to imFindCells - MA
%       bwCell_vol_numbered: bw cell mask, numbered 1, .., N for each 3D
%       cluster - MA
%       zoomToAxis: 4-element vector to be passed to AXIS, sets image zoom
%           default: []
%       cellRadius: how large of a spot to fill cells with
%           default: 1; usually you should not need to change this
%             (if you find that clicks are too sensitive to exact pixel
%             position it can help to increase this to 2)
%             set to empty to auto-compute from existing cells; this is
%             rarely helpful
%       Im_contrast: for viewing only (doesn't affect analysis), set the
%       contrast between .51 and 1, 1 scale each image to MIN/MAX of
%       volume, lower values increase contrast - MA
%
%   Hints for using:
%      Contrast threshold (fract of mean) should be between 0.9 and 1.05,
%          usually between 0.95 and 1.0.  Higher values mean more
%          restrictive/smaller masks, and avoid 'spillover'.
%          I usually use 0.95 for dim cells and 1.0 for bright cells.
%      Dilation disk size in pix: 1-3, depending on cell sizes.  For cell
%          radii of ~10 pix (512pix,360um FOV), I use 3
%          for 3-5 pix radii, I use 1
%
%      * If two cells are nearby, you can often avoid spillover by clicking
%        the bright cell before the dim cell; likewise you can sometimes use
%        bright stuff between two cells to segment them, then delete it
%        afterwards.
%      * If you're not getting a big enough mask for a cell whose center is
%        bright, try clicking a bit off center on the dimmer edge.

%   MH 5/1/08: add delete option; ui for editing theshold and disk size
%

% SY - 01/25/2010 - fixed some bugs, use sparse matrix for undo command
% SY - 02/02/2010 - removed uno, added "cell center" output

%$Id: imCellEditInteractive.m 394 2008-11-21 07:02:34Z histed $


%make an fovImg from the volume: MA
%plan: click on xy, then click on either xz or yz to decide [X Y Z]
global bwCell_vol;


X1 = 20;
Y1 = 20;
Z1 = 3;

if ischar(fovImg_vol) % if filename
    % read tif file info
    info = imfinfo(fovImg_vol,'tif');
    W=info(1).Width;
    H=info(1).Height;
    nFrames=length(info);
    tmp_im= imread(fovImg_vol,1);
    data_class=class(tmp_im(1));

    %preallocate color stacks
    Red_im=zeros(H,W,nFrames,data_class);
    Green_im=zeros(H,W,nFrames,data_class);
    Blue_im=zeros(H,W,nFrames,data_class);

    %read frames, separete color chennels
    for frame=1:nFrames
        tmp_im= imread(fovImg_vol,frame);
        Red_im(:,:,frame)=tmp_im(:,:,1);
        Green_im(:,:,frame)=tmp_im(:,:,2);
        Blue_im(:,:,frame)=tmp_im(:,:,3);
    end
%fovImg_vol=(Red_im+Green_im+Blue_im)/3;
fovImg_vol=(Red_im);

clear Red_im
clear Green_im
clear Blue_im    
end


Zproj = squeeze(fovImg_vol(:,:,Z1));
Xproj = squeeze(fovImg_vol(Y1,:,:));
Yproj = squeeze(fovImg_vol(:,X1,:));

%IMG_CONTRAST = .8; %.8 of max across entire volume, [0 1]
if nargin < 5 || isempty(Im_contrast)
    %default is 1, scale min->max of volume
    IMG_CONTRAST_INV=[0 0];
elseif length(Im_contrast)==1;
    IMG_CONTRAST_INV = [1 - Im_contrast 1 - Im_contrast];
else
    IMG_CONTRAST_INV = [Im_contrast(1) 1 - Im_contrast(2)];
end


MINMAX0 = [min(fovImg_vol(:)) max(fovImg_vol(:))]; %SY
MINMAX = [(MINMAX0(1) + IMG_CONTRAST_INV(1)*diff(MINMAX0)) ((MINMAX0(2) - IMG_CONTRAST_INV(2)*diff(MINMAX0)))];

if nargin<6
    resolution=[1 1 1];
end




%cludgy way to keep scaling on images the same across planes..
fovImg_view = [flipud(Xproj') zeros(size(Yproj,2),size(Xproj,2))'; Zproj Yproj];
fovImg_view(1,1) = MINMAX(1);
fovImg_view(1,2) = MINMAX(2);
fovImg_view(fovImg_view<MINMAX(1)) = MINMAX(1); %SY
fovImg_view(fovImg_view>MINMAX(2)) = MINMAX(2); %SY

fovImg = Zproj;

[nRows0 nCols0 nZ] = size(fovImg_vol);

[nRows nCols] = size(fovImg_view);
X_skip = nRows - nRows0;

if ~isempty(bwCell_vol_numbered);
    bwCell_vol_numbered = bwlabeln(bwCell_vol_numbered>0,6);
    bwCell_vol = bwCell_vol_numbered> 0;
else
    bwCell_vol = zeros([nRows0,nCols0, nZ]);
    bwCell_vol_numbered = zeros([nRows0,nCols0, nZ]);
end

bwCell = squeeze(bwCell_vol(:,:,Z1));
Zproj = squeeze(max(bwCell_vol,[],3));
Xproj = squeeze(bwCell_vol(Y1,:,:));
Yproj = squeeze(bwCell_vol(:,X1,:));
bwCell_view = [flipud(Xproj') zeros(size(Yproj,2),size(Xproj,2))'; Zproj Yproj];

labCell = bwlabel(bwCell,8);

%% measure existing cells
stats = regionprops(labCell, 'Area');
cellAreas = cat(1, stats.Area);

bwvol = bwlabeln(bwCell_vol,6);
nCells = max(bwvol(:));

%% process arguments
if nargin < 3 || isempty(zoomToAxis),
    % add 0.5 margin at each edge to emulate imshow / imagesc
    zoomToAxis = [1 nCols 1 nRows] + 0.5+[-1 1 -1 1];
end
if nargin < 4, cellRadius = 1; end

%% auto-computer the cell area dilation disk radius, only if requested
if isempty(cellRadius)
    if nCells < 5
        error('Too few cells to measure cell radius (%d): ', nCells);
    else
        cellRadius = sqrt(mean(cellAreas)/pi);
    end
end
cellRadius = round(cellRadius);

% draw initial figure
figH = figure;

imH = subDrawShaded(fovImg_view, bwCell_view, zoomToAxis, X1, Y1 + X_skip, Z1, nRows0, nCols0, nZ, MINMAX);
axH = get(imH, 'Parent');
set(axH, 'Position', [0.06 0.06 0.88 0.88])

% set up UI
figPos = get(figH, 'Position');
figColor = get(figH, 'Color');
set(axH, 'Units', 'normalized');
axPos = get(axH, 'Position');
utlH = uicontrol('Style', 'text', ...
    'String', {'Contrast threshold:', ' fract of mean'}, ...
    'Position', [5 figPos(4)-40, 80, 40], ...
    'BackgroundColor', figColor, ...
    'HorizontalAlignment', 'left');
utH = uicontrol('Style', 'edit', ...
    'String', '0.95', ...   % default cThresh
    'Units', 'pixels', ...
    'Position', [5 figPos(4)-70, 60, 30]);

udlH = uicontrol('Style', 'text', ...
    'String', {'Dilation disk:', ' size in pix'}, ...
    'Position', [5 figPos(4)-130, 80, 30], ...
    'BackgroundColor', figColor, ...
    'HorizontalAlignment', 'left');
udH = uicontrol('Style', 'edit', ...
    'String', '1', ...    % default diskR
    'Units', 'pixels', ...
    'Position', [5 figPos(4)-160, 60, 30]);
cdH = uicontrol('Style', 'text', ...
    'String', { 'Cell radius:', ...
    ' (initial disk)', ...
    sprintf(' %5.2fpix', cellRadius) }, ...
    'Units', 'pixels', ...
    'Position', [5 figPos(4)-260, 60, 60], ...
    'BackgroundColor', figColor, ...
    'HorizontalAlignment', 'left');


msgH = uicontrol('Style', 'text', ...
    'String', {'Msg:', ' N/A'}, ...
    'Position', [5 figPos(4)-360, 80, 30], ...
    'BackgroundColor', figColor, ...
    'HorizontalAlignment', 'left');

% title str
tStr = { [mfilename ': Click on a cell region to add'], ...
    '       rigth-click to del, q / ret to finish, u to undo' };
tHeight = (1-axPos(4))-axPos(2);
tlH = uicontrol('Style', 'text', ...
    'String', tStr, ...
    'Units', 'normalized', ...
    'Position', [0.25 1-tHeight, 0.5, tHeight*0.8], ...
    'BackgroundColor', figColor, ...
    'HorizontalAlignment', 'left');


%% iterate: get a point, add it, display it
bwCurr = bwCell;

%nActions = 0;
nTotal = nCells;

while 1  % till broken out of
    % interactively get clicks
    [X Y0 selectionType] = getAPoint(gca);
    set(msgH, 'String','');
    Y = Y0 - X_skip;
    if isnan(X) %keyboard input
        key = lower(Y0);
        if isempty(key)
            key = 'NaN';
            continue
        end

        switch key
            case char(13) % return
                break;  % out of loop, done
            case 'q' %quit
                %prune the values
                %unique(nonzeros(bwCell_vol_numbered))
                tmp2=unique(nonzeros(bwCell_vol_numbered));
                bwCell_vol_numbered2 = zeros(size(bwCell_vol_numbered));
                for count5 = 1:length(tmp2)
                    bwCell_vol_numbered2(bwCell_vol_numbered == tmp2(count5)) = count5;
                end
                bwCell_vol_numbered = bwCell_vol_numbered2;
                break;
%             case 'u'
%                 % undo
%                 if nActions == 0
%                     printf(1, '** Cannot undo: no cells added yet\n');
%                     set(msgH, 'String','Cannot undo: no cells added yet');
%                     continue;
%                 end
%                 bwCell_vol_numbered(:)=full(bwSave{nActions-1});
%                 bwCell_vol = bwCell_vol_numbered>0;
%                 bwSave{nActions}=[];
%                 cell_pos(nActions,:)=[];
%                 nActions = nActions-1;
%                 
%                 %prune the values
%                 tmp2=unique(nonzeros(bwCell_vol_numbered));
%                 nTotal=length(tmp2);
%                 bwCell_vol_numbered2 = zeros(size(bwCell_vol_numbered));
%                 for count5 = 1:nTotal
%                     bwCell_vol_numbered2(bwCell_vol_numbered == tmp2(count5)) = count5;
%                 end
%                 bwCell_vol_numbered = bwCell_vol_numbered2;
%                 clear bwCell_vol_numbered2
%                 fprintf(1, 'Undo!  %d cells total now\n', nTotal);
%                 set(msgH, 'String',['Undo! Now cells total : ' num2str(nTotal)]);
%                 continue
            case 'z'
                Z1a = Z1 + 1;
                if Z1a <= nZ && Z1a > 0
                    Z1 = Z1a;
                end
                fprintf(1, ['New Z position set: ',num2str(Z1),'\n']);
                set(msgH, 'String',['New Z position set: ' num2str(Z1)]);
            case 'x'
                %move up
                Z1a = Z1 -1;
                if Z1a <= nZ && Z1a > 0
                    Z1 = Z1a;
                end
                fprintf(1, ['New Z position set: ',num2str(Z1),'\n']);
                set(msgH, 'String',['New Z position set: ' num2str(Z1)]);
        end
    else %mouse cilck
        %% validate XY point
        X = round(X);
        Y = round(Y);
        if X <= 0 || X >= nCols || Y <= (-1*X_skip) || Y >= nRows0
            % clicked outside axes, repeat
            fprintf(1, 'Click was outside image axes, try again\n');
            set(msgH, 'String','Click was outside image axes, try again');
            continue;
        end

        %shift Z axis if click is on XZ or YZ section
        if (((Y > (-1*X_skip)) && (Y <= 0)) || ((X > nCols0) && (X <= nCols)))
            if ((Y > (-1*X_skip)) && (Y <= 0))
                Z1 = nZ + 1 - floor(Y0);
            elseif ((X > nCols0) && (X <= nCols))
                Z1 =  (X - nCols0);
            end

            if X>0 && X<=nCols0
                X1=X;
            end
        else %click is on XY section
            %%% get ui data
            diskR = str2double(get(udH, 'String'));
            if isnan(diskR)
                fprintf(1, 'Error reading disk radius: %s, try again\n', ...
                    get(udH, 'String'));
                set(msgH, 'String','Error reading disk radius: try again');
                continue
            end
            if (diskR < 1) || diskR > 50
                fprintf(1, 'Disk too small or too big, try again\n');
                continue
            elseif ~iswithintol(round(diskR), diskR, 10^2*eps)
                fprintf(1, 'Disk radius must be an integer, try again\n');
                set(msgH, 'String','Disk radius must be an integer, try again');
                continue
            end

            cThresh = str2double(get(utH, 'String'));
            if isnan(cThresh)
                fprintf(1, 'Error reading threshold: %s, try again\n', ...
                    get(utH, 'String'));
                set(msgH, 'String','Error reading threshold: try again');
                continue
            end
            if cThresh <= 0 || cThresh >= 100
                fprintf(1, 'Threshold too small or too big, try again\n');
                set(msgH, 'String','Threshold too small or too big, try again');
                continue
            end

            %%% what kind of mouse click?
            switch lower(selectionType)
                case 'alt'    % rigth-click: delete
                    % make sure we're in a cell
                    if bwCurr(Y,X) ~= 1
                        fprintf(1, '** Trying to delete, not on a cell, try again\n');
                        set(msgH, 'String','Trying to delete, not on a cell, try again');
                        continue;
                    else
                        %% do delete
                        cell_pos(bwCell_vol_numbered(Y,X,Z1),:)=[];
                        bwCell_vol_numbered(bwCell_vol_numbered == bwCell_vol_numbered(Y,X,Z1)) = 0;                                               
                        %prune the values
                        tmp2=unique(nonzeros(bwCell_vol_numbered));
                        nTotal=length(tmp2);
                        bwCell_vol_numbered2 = zeros(size(bwCell_vol_numbered));
                        for count5 = 1:nTotal
                            bwCell_vol_numbered2(bwCell_vol_numbered == tmp2(count5)) = count5;
                        end
                        bwCell_vol_numbered = bwCell_vol_numbered2;
                        clear bwCell_vol_numbered2

                        bwCell_vol = bwCell_vol_numbered>0;
                        bwCurr = squeeze(bwCell_vol(:,:,Z1));
                    %    nTotal = nTotal-1;
                        fprintf(1, 'Deleted object, %d total remain\n', nTotal);
                        set(msgH, 'String',['Deleted object, total remain: ' num2str(nTotal)]);
                    end

                case 'normal'    % left-click: add
                    % make sure not in a cell
                    if bwCurr(Y,X) == 1
                        fprintf(1, '** Trying to add in existing cell region, try again\n');
                        set(msgH, 'String','Trying to add in existing cell region, try again');
                        continue;
                    else

                        bwNew_vol0 = zeros(size(bwCell_vol));
                        bwNew_vol = bwCell_vol;
                        %fill in whole volume, use same mean as mean from current plane to get whole cell:
                        if cellRadius>0
                        [bwNew bwCell cellMean] = subAddCell(bwCurr, fovImg, X, Y, ...
                            cellRadius, cThresh, diskR);
                        end
                        
                    end
                    if cellRadius>0
                    bwCell_tmp=bwCell;
                    %step up and down in Z
                    for step_Z = [1 -1]
                        bwCell=bwCell_tmp;
                        STOP = 0;
                        if step_Z == 1
                            Z1_use = Z1;
                        elseif step_Z == -1
                            Z1_use = Z1 - 1;
                        end
                        %Z1_use = Z1_use + step_Z;
                        while STOP ~= 1 && (Z1_use > 0) && (Z1_use <= nZ)
                            fovImg_USE = squeeze(fovImg_vol(:,:,Z1_use));
                            bwCurr_USE = squeeze(bwNew_vol(:,:,Z1_use));
                            %                            [bwNew bwCell] = subAddCell(bwCurr_USE, fovImg_USE, X, Y, ...
                            %                                cellRadius, cThresh, diskR, cellMean);
                            cell_proj=double(fovImg_USE).*double(bwCell); %SY
                            max_cell_proj=max(cell_proj(:)); %SY
                            [Y2 X2]=find(cell_proj==max_cell_proj); %SY
                            [bwNew bwCell] = subAddCell(bwCurr_USE, fovImg_USE, X2(1), Y2(1), ...
                                cellRadius, cThresh, diskR, cellMean);
                            %check the plane above/below to make sure no touching..
                            Z1_use_tmp = Z1_use + step_Z;
                            %first, check if this plane exists
                            if (Z1_use_tmp > 0) && (Z1_use_tmp <= nZ)
                                bwCell_abovebelow =  squeeze(bwNew_vol(:,:,Z1_use_tmp));
                                %check for direct contact between two cells (diagonal
                                %through Z doesn't get checked, could dilate to do
                                %this)
                                tmp = sum(sum(bwCell_abovebelow.*bwCell));
                                if tmp > 0
                                    STOP = 1;
                                end
                            end

                            if sum(sum(bwCell)) > 0 && STOP == 0
                                %msg='added'
                                bwNew_vol(:,:,Z1_use) = bwNew;
                                bwNew_vol0(:,:,Z1_use) = bwCell;
                                Z1_use = Z1_use + step_Z;
                            else
                                STOP = 1;
                            end

                        end
                    end


                    bwCell_vol_numbered = bwCell_vol_numbered + (nTotal+1)*bwNew_vol0;
                    bwCell_vol = bwCell_vol_numbered>0;
                    bwCurr = squeeze(bwCell_vol(:,:,Z1));
                    else
                        Yrange=max(Y+round(cellRadius/resolution(2)),1):min(Y-round(cellRadius/resolution(2)),nRows0);
                        Xrange=max(X+round(cellRadius/resolution(1)),1):min(X-round(cellRadius/resolution(1)),nCols0);
                        Zrange=max(Z1+round(cellRadius/resolution(3)),1):min(Z1-round(cellRadius/resolution(3)),nZ);
                        
                    bwNew_vol0(Yrange,Xrange,Zrange)=1;
                    bwCell_vol_numbered(Yrange,Xrange,Zrange)=nTotal+1;
                    bwCell_vol(Yrange,Xrange,Zrange)=1;
                    bwCurr = squeeze(bwCell_vol(:,:,Z1));
                    end                    
                    
                    nTotal = nTotal+1;
                    fprintf(1, 'Added object #%d: %d pix\n',nTotal, sum(bwNew_vol0(:)));
                    set(msgH, 'String',['Added object #' num2str(nTotal) '  pix ' num2str(sum(bwNew_vol0(:)))]);
                    cell_pos(nTotal,1:3)=[X*resolution(1) Y*resolution(2) Z1*resolution(3)];
                otherwise
                    % other type of click, just continue without comment
                    %keyboard
                    fprintf(1, 'Unrecognized click occurred (matlab bug?) %s, %g %g\n', ...
                        selectionType, X, Y);
                    set(msgH, 'String',['Unrecognized click occurred (matlab bug?) ' selectionType ' X,Y ' num2str(X) ',' num2str(Y)]);
                    continue
            end

            %% save old mask and update
            %bwSave{nActions+1} = bwCell_vol_numbered;
%            bwSave{nActions+1} = sparse(bwCell_vol_numbered(:));
%            nActions = nActions+1;
            X1 = X;
            Y1 = Y;
        end
    end
    Zproj = squeeze(bwCell_vol(:,:,Z1));
    Xproj = squeeze(bwCell_vol(Y1,:,:));
    Yproj = squeeze(bwCell_vol(:,X1,:));
    bwCurr_view = [flipud(Xproj') zeros(size(Yproj,2),size(Xproj,2))'; Zproj Yproj];

    Zproj = squeeze(fovImg_vol(:,:,Z1));
    Xproj = squeeze(fovImg_vol(Y1,:,:));
    Yproj = squeeze(fovImg_vol(:,X1,:));
    fovImg_view = [flipud(Xproj') zeros(size(Yproj,2),size(Xproj,2))'; Zproj Yproj];
    fovImg = Zproj;
    fovImg_view(1,1) = MINMAX(1);
    fovImg_view(1,2) = MINMAX(2);
    fovImg_view(fovImg_view<MINMAX(1)) = MINMAX(1);
    fovImg_view(fovImg_view>MINMAX(2)) = MINMAX(2);

    % draw the new version
    subDrawShaded(fovImg_view, bwCurr_view, zoomToAxis, X1, Y1 + X_skip, Z1, nRows0, nCols0, nZ, MINMAX);
end

%bwOut = bwCurr;
fprintf(1, '%s: done, %d objects total).\n', ...
    mfilename, nTotal);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function imH = subDrawShaded(fovImg, bw, zoomToAxis, Xuse, Yuse, Zuse, nRows0, nCols0, nZ, MINMAX)

shadeimg = imShade(fovImg, bw);
cla;
imH = imagesc(shadeimg,MINMAX);
set(gca, 'XGrid', 'on', ...
    'YGrid', 'on', ...
    'Visible', 'on', ...
    'YDir', 'reverse', ...
    'DataAspectRatio', [1 1 1]);
axis(zoomToAxis);  %[397 468 212 283]);

line(Xuse*ones(1,2),zoomToAxis(3:4),'Color','w','LineWidth',1);
line(zoomToAxis(3:4),Yuse*ones(1,2),'Color','w','LineWidth',1);
Zuse2 = (nZ + 1 - Zuse);
Zuse3 = nCols0 + (nZ + 1 - Zuse2);
line(Zuse3*ones(1,2),zoomToAxis(3:4),'Color','m');
line(zoomToAxis(1:2),Zuse2*ones(1,2),'Color','m');


drawnow;

function [bwNew bwCell cellMean] = subAddCell(bwCurr, fovImg, X, Y, ...
    cellRadius, cThresh, diskR, cellMean)

%%% Real cell-adding logic is here
[nRows nCols] = size(bwCurr);

% set up dilation disks:  use 4-connected regions for all
%se = strel('disk',1,8);  % for cells-to-avoid
se = strel('square',3);  % for cells-to-avoid: all 9 pix in a square;
% avoid diagonally-connected cells.
se2 = strel('disk',round(cellRadius),4);  % for region to find mean over
seJunk = strel('disk', max(round(cellRadius/4), 1), 4);  % remove thin junk
seExpand = strel('disk', diskR, 4);  % expand thresholded region

% add a disk around start point, non-overlapping with adj cells
tempmask = false(nRows, nCols);
dilateorg = imdilate(bwCurr,se);
tempmask(Y, X) = 1;
tempmask = imdilate(tempmask,se2);
tempmask = tempmask & ~dilateorg;

% fill region around disk of similar intensity, combine with disk
if nargin < 8
    cellMean = mean(fovImg(tempmask == 1),1);
end

if fovImg(Y, X) > cellMean.*cThresh
    allMeanBw = fovImg >= cellMean.*cThresh;  % threshold by intensity
    %Npts = sum(sum(allMeanBw))
    connMeanBw = bwselect(allMeanBw &~dilateorg, X, Y, 4);
    connMeanBw = connMeanBw |tempmask & ~dilateorg;

    % erode then dilate filled to remove sharp things
    erMean = imerode(connMeanBw, seJunk);
    dilateMean = imdilate(erMean, seJunk);
    dilateMean = imdilate(dilateMean, seExpand); % (thresh is conservative)
    bwCell = dilateMean & ~dilateorg;

    %    bwNew = bwCurr | bwCell;
else
    bwCell = zeros(size(bwCurr));
    %    bwNew = bwCurr;
end
bwNew = bwCurr | bwCell;

