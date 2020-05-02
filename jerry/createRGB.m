function [rgb_reg2ref, rgb_reg2ref_dfof, yellowTotal] = createRGB(ref, reg2ref, reg2ref_dfof, data_dfof_max_ref, sz_target)
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
    hsv = rgb2hsv(rgb_reg2ref_dfof); %create hue, saturation, and value of rgb_reg2ref_dfof image
    hue = hsv(:,:,1); %extract hue from hsv
    hue = hue > 49/360 & hue < 65/360; %thresholding for yellow hue
    yellowTotal = sum(hue(:)) %summing total pixels with yellow hue
%     yellowpixels = rgb_reg2ref_dfof(:,:,1) > 0.78 & rgb_reg2ref_dfof(:,:,2) > 0.78 & rgb_reg2ref_dfof(:,:,3) < 0.4; 
%     yellowTotal = sum(yellowpixels(:))
    end