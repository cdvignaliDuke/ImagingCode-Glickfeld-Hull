function [fov_siz um_per_pix] = frGetFovSizeInvitro(zoom_v,res);

fov_siz = 691*(zoom_v^-1);

um_per_pix = fov_siz/res

return;

% calibration with 40X obj

zoom_mat = [1 2 3 4 5 10];
pix_per_25um_mat = [19 36 52 71 89 177];
fov_mat = (25./pix_per_25um_mat).*512;
b = regress(log10(zoom_mat'), log10(fov_mat'))

