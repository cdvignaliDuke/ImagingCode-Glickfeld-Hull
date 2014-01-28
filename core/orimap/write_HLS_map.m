function HLS = write_HLS_map(params, fname, max, tune_max)

% hue: preferred angle
% intenisity: max_change
% saturation: tune (0-tune_max)
% HLS_hc: high contrast HLS map

dim=size(params.th);

m1=params.max_change./max;
m2=m1*2;
m1(find(m1>1))=1;
m2(find(m2>1))=1;
s=params.tune./tune_max;
s(find(s>1))=1;
temp=ones(dim(1),dim(2));

HLS=hsv2rgbKO(params.th,s,m1);
HLS_hc=hsv2rgbKO(params.th,s,m2);


imwrite(HLS, [fname, 'HLS.tif']);
imwrite(HLS_hc, [fname, 'HLS_hc.tif']);
