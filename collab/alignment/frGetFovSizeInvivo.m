function [fov_siz um_per_pix] = frGetFovSizeInvitro(zoom_v,res);

fov_siz = 103*(zoom_v) + 15.25;

um_per_pix = fov_siz/res

return;

% calibration with 40X obj

zoom_mat = [4 3 2 1 .5 .25];
pix_per_25um_mat = [15 20 28 52 98 178];
fov_mat = (25./pix_per_25um_mat).*256;


