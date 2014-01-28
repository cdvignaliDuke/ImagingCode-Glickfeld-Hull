frames=500;
down=1 ;
nPlanes = 8;
base = 'G:\users\lindsey\analysisLG\active mice\LG27\110226';
expt = '110226_LG27_850nm_RG_Z6';
channel1 = '_green';
channel2 = '_red';
fn_1= fullfile(base, [expt channel1 '.tif']);
fn_2= fullfile(base, [expt channel2 '.tif']);
outDir = 'G:\users\lindsey\analysisLG\active mice\LG27\110226';
zstack_g = readtiff(fn_1);
zstack_green = stackGroupProject(zstack_g,down);
zstack_r = readtiff(fn_2);
zstack_red = stackGroupProject(zstack_r,down);

% register 1st green plane
iPlane = 5;
zstack_sub_green = zeros(240,256,frames/down);
zstack_sub_green = zstack_green(:,:,iPlane:nPlanes:end);
av = mean(zstack_sub_green(:,:,1:10),3);
[out,reg] = stackRegister(zstack_sub_green,av);
zstack2_green = zeros(240,256,nPlanes);
zstack2_green(:,:,iPlane) = mean(reg,3);

%register 1st red plane
zstack_sub_red = zeros(240,256,(frames/down));
zstack_sub_red = zstack_red(:,:,iPlane:nPlanes:end);
zstack2_red = stackShifts(double(zstack_sub_red), -out(:,4:-1:3));
zstack3_red = zeros(240,256,nPlanes);
zstack3_red(:,:,iPlane) = mean(zstack2_red,3);


% register other planes to preceding plane 
for iPlane = [6 7 8 1]
    av = mean(reg,3);
    zstack_sub_green = zstack_green(:,:,iPlane:nPlanes:end);
    [out,reg] = stackRegister(zstack_sub_green,av);
    zstack2_green(:,:,iPlane)= mean(reg,3);
    zstack_sub_red = zstack_red(:,:,iPlane:nPlanes:end);
    zstack2_red = stackShifts(double(zstack_sub_red), -out(:,4:-1:3));
    zstack3_red(:,:,iPlane) = mean(zstack2_red,3);
end
fn_out_green = fullfile(outDir,sprintf('%s_green_sub_avg_reg.tif',expt));
writetiff(zstack2_green, fn_out_green);
fn_out_red = fullfile(outDir,sprintf('%s_red_sub_avg_reg.tif',expt));
writetiff(zstack3_red, fn_out_red);