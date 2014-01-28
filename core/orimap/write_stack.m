function write_stack (stack, fname)
% write 3 dimensional matrix to many 2-dimensional tif files.
% pixel value will be between 0 and 255.
% stack: 3 dimensional matrix. Pixel value between 0 and 1 will be converted to 0-1.
% fname: e.g. 'C:\share\OG031210\OG031210_6\OG031210_6_vstim12_t'
% this program requires numextCR to run. 

dim=size(stack);
for fr = 1:dim(3)
    filename=[fname,numextCR(fr-1,3),'.tif'];
    imwrite(stack(:,:,fr),filename,'Compression','none');
end;