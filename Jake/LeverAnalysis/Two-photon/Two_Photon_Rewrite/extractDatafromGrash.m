open('img92_maskraw_Day2.fig');
h = gcf;
axesObjs = get(h, 'Children');
dataObjs = get(axesObjs, 'Children');
objTypes = get(dataObjs, 'Type');
mask = get(dataObjs, 'cData');
mask_final = reshape(mask,1, 264*796);
save('mask_final.mat', 'mask_final');