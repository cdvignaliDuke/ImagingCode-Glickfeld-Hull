function polar = write_polar_map(params, fname, max)

% hue: preferred angle
% intenisity: vector magnitude
% polar_hc: high contrast polar map

dim=size(params.th);

m1=params.mag./max;
m2=m1*2;
m1(find(m1>1))=1;
m2(find(m2>1))=1;
temp=ones(dim(1),dim(2));

polar=hsv2rgbKO(params.th,temp,m1);
polar_hc=hsv2rgbKO(params.th,temp,m2);

imwrite(polar, [fname, 'polar.tif']);
imwrite(polar_hc, [fname, 'polar_hc.tif']);
