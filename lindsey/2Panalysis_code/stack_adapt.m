%% contrast adapt stack
dir = 'G:\users\lindsey\analysisLG\active mice\LG16\110219\';
file = '110219_LG16_850nm_zpt25_C_red_avg_reg.tif';
stack = readtiff([dir file]);
stack2 = zeros(size(stack));
stack2(:,:,10) = stack(:,:,10);
initial_mean = mean(stack(:,:,10));
for iPlane = 2:size(stack,3);
    new_mean = mean(stack(:,:,iPlane));
    scale = initial_mean/new_mean;
    stack2(:,:,iPlane) = scale*stack(:,:,iPlane);
end
fn_out = fullfile(dir,sprintf('110219_LG16_850nm_zpt25_C_red_avg_reg_adapt.tif'));
writetiff(stack2, fn_out);