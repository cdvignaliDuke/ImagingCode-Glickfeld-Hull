function [] = write_intensity_maps (params, fname, max);

% write selectivity parameter images as tif files
% auto scaled between zero and max value (for mag, ave_change, max_change)
% tune is between 0-1, without normalization.

imwrite(params.mag./max, [fname, 'mag.tif']);
imwrite(params.ave_change./max, [fname, 'ave_change.tif']);
imwrite(params.max_change./max, [fname, 'max_change.tif']);
imwrite(params.tune, [fname, 'tune.tif']);
