runFolder = expt(iexp).redChannelRun;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate);

%% load red data
fName = [runFolder '_000_000'];
redData = loadsbx_choosepmt(2,mouse,expDate,runFolder,fName);
if size(redData,3) > 1000
    redData = redData(:,:,1:1000);
end
%% register
load(fullfile(fn,'regOuts&Img.mat'))
redDataDownSamp = stackGroupProject(redData,10);
[~,redReg] = stackRegister(redDataDownSamp(:,:,1),data_corr_img);

redImage = mean(redReg,3);
figure;colormap hot; imagesc(redImage)

for i = 1:size(redReg,3)
    imagesc(redReg(:,:,i))
    drawnow
end
%% save image
save(fullfile(fn,'redChannelImage'),'redImage')

%% select cells
bwout = imCellEditInteractive(redImage);
red_mask = bwlabel(bwout);
figure;imagesc(red_mask);title(sprintf('%s cells',num2str(length(unique(red_mask(:)))-1)))

%% match gcamp mask with red label mask
load(fullfile(fn,'final_mask.mat'))
activeMask_linear = mask_cell(:);
nc = length(unique(red_mask(:)))-1;
redCellIdentity = zeros(1,nc);
for icell = 1:nc
    redpix = red_mask(:) == icell;
    matchedpix = activeMask_linear(redpix);
    isMatch = any(matchedpix > 0);
    if isMatch
        activeIdentity = unique(matchedpix(matchedpix > 0));
        if length(activeIdentity) > 1
            warning('two cells identified, red cell skipped')
            continue
        else
            redCellIdentity(icell) = activeIdentity;
        end
    end
end
save(fullfile(fn,'redCellIdentity'),'redCellIdentity')