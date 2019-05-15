fn = 'Z:\home\ashley\Manuscripts\Attention V1\Images';
antiImg = readtiff(fullfile(fn,'maxDFF_2Pimage_AW25_160309_bx_anti.tif'));
tarImg = readtiff(fullfile(fn,'maxDFF_2Pimage_AW25_160309_bx_tar.tif'));
load(fullfile('Z:\home\ashley\.TemporaryItems\Analysis\AW25\two-photon imaging\160309\data processing','final_mask_cells'))
antiShade = imShade(antiImg,mask_cell);
% antiShade(:,:,2) = antiShade(:,:,2) .* (1-.75*double(mask_cell));
figure;imagesc(antiShade) 
tarShade = imShade(tarImg,mask_cell);
figure;imagesc(tarShade)
% writetiff(antiShade,fullfile(fn,'maxDFF_2Pimage_AW25_160309_bx_anti_shade'))
% writetiff(tarShade,fullfile(fn,'maxDFF_2Pimage_AW25_160309_bx_tar_shade'))
imwrite(antiShade,fullfile(fn,'maxDFF_2Pimage_AW25_160309_bx_anti_shade.png'))
imwrite(tarShade,fullfile(fn,'maxDFF_2Pimage_AW25_160309_bx_tar_shade.png'))
