function rng = rangeMinMax(data)

rng1 = floor(min(data(:)));
rng2 = ceil(max(data(:)));
rng = [rng1, rng2];


end