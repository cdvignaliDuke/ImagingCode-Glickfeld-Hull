function img_adapt = padnadapt(img)
%PADNADAPT
%IMG_ADAPT = PADNADAPT(IMG)

pad = 30;

%img_med = medfilt2(img,[2,2]);
img_med = img;         
img_pad = padarray(img_med, pad*[1 1], 'symmetric','both');
img_full = imScale(stackLocalContrastAdaptation(img_pad, pad, 1));
img_adapt = imscale(img_full(pad+1:end-pad,pad+1:end-pad));
return;
