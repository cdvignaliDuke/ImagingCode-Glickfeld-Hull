frames=100;
down=1 ;
base = 'G:\users\lindsey\dataLG';
date = '110512';
mouse = 'X32';
expt = 'run5';
channel1 = '_green';
channel2 = '_red';

fn_g= fullfile(base, [date '_' mouse], expt, [expt channel1]);
fn_r= fullfile(base, [date '_' mouse], expt, [expt channel2]);
outDir = ['G:\users\lindsey\analysisLG\active mice\' mouse '\' date '\'];




zstack_green = readtiff(fn_g);

zstack_red = readtiff(fn_r);


% register 1st green plane
[a,b,c] = size(zstack_green);
nPlanes = c/(frames/down);
iPlane = 1;
zstack_sub_green = zeros(240,256,frames/down);
zstack_sub_green = zstack_green(:,:,((iPlane-1)*(frames/down))+1:iPlane*(frames/down));
av = mean(zstack_sub_green(:,:,1:10),3);
[out,reg] = stackRegister(zstack_sub_green,av,1);
zstack2_green = zeros(240,256,nPlanes);
zstack2_green(:,:,iPlane) = mean(reg,3);

%register 1st red plane
zstack_sub_red = zeros(240,256,(frames/down));
zstack_sub_red = zstack_red(:,:,((iPlane-1)*(frames/down))+1:iPlane*(frames/down));
zstack2_red = stackShifts(double(zstack_sub_red), -out(:,4:-1:3));
zstack3_red = zeros(240,256,nPlanes);
zstack3_red(:,:,iPlane) = mean(zstack2_red,3);


% register other planes to preceding plane 
for iPlane = 2:nPlanes
    av = mean(reg,3);
    zstack_sub_green = zstack_green(:,:,((iPlane-1)*(frames/down))+1:iPlane*(frames/down));
    [out,reg] = stackRegister(zstack_sub_green,av);
    zstack2_green(:,:,iPlane)= mean(reg,3);
    zstack_sub_red = zstack_red(:,:,((iPlane-1)*(frames/down))+1:iPlane*(frames/down));
    zstack2_red = stackshifts(double(zstack_sub_red), -out(:,4:-1:3));
    zstack3_red(:,:,iPlane) = mean(zstack2_red,3);
end
fn_out_green = fullfile(outDir,sprintf('%s_green_avg_reg.tif',expt));
writetiff(zstack2_green, fn_out_green);
fn_out_red = fullfile(outDir,sprintf('%s_red_avg_reg.tif',expt));
writetiff(zstack3_red, fn_out_red);