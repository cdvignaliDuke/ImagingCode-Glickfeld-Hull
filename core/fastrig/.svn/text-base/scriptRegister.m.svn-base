function scriptRegister(greenSourcePn,redSourcePn,greenTargetPn,redTargetPn,shiftsTargetPn)
% scriptRegister(greenSourcePn,redSourcePn,greenTargetPn,redTargetPn,shiftsTargetPn)

dirs = frGetDirs;

%% registration / multicore
green = readtiff(greenSourcePn);
red = readtiff(redSourcePn);

[ny,nx,nframes]=size(green);

siz = 128;
sel = floor(nframes/2):floor(nframes/2)+siz-1;
target = mean(green(:,:,sel)+red(:,:,sel),3);

delete(fullfile(dirs.multicore,'*.*'));

out = RegMulticore_shifts(green+red, target, shiftsTargetPn,dirs.multicore, 200);

green_reg = stackShifts(double(green), out);
red_reg = stackShifts(double(red), out);

green_reg(find(green_reg(:)==nan)) = 0;
red_reg(find(red_reg(:)==nan)) = 0;

writetiff(green_reg,greenTargetPn,'uint16');
writetiff(red_reg,redTargetPn,'uint16');

return;