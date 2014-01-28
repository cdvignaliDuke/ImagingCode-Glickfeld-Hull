frames=100;
base = '\\zoloto\bigstorlab\Lindsey\fastrig\2012\August';
out = 'G:\users\lindsey\analysisLG\active mice\';
mouse = 'CM124';
date = '120808';
expt = 'stack_z4';
channel = '_green';
fn_1= fullfile(base, [date '_' mouse], expt, [expt channel]);
outDir = fullfile(out, mouse,date);
zstack_g = readtiff(fn_1);

% register 1st green plane
[a,b,c] = size(zstack_g);
nPlanes = c/frames;
iPlane = 1;
zstack_sub_green = zeros(240,256,frames);
zstack_sub_green = zstack_g(:,:,((iPlane-1)*(frames))+1:iPlane*(frames));
av = mean(zstack_sub_green(:,:,1:10),3);
[out,reg] = stackRegister(zstack_sub_green,av,1);
zstack2_green = zeros(240,256,nPlanes);
zstack2_green(:,:,iPlane) = mean(reg,3);
av = mean(reg,3);

% register other planes to preceding plane 
for iPlane = 2:nPlanes
    av = mean(reg,3);
    zstack_sub_green = zstack_g(:,:,((iPlane-1)*(frames))+1:iPlane*(frames));
    [out,reg] = stackRegister(zstack_sub_green,av,1);
    zstack2_green(:,:,iPlane)= mean(reg,3);
    zstack2_long = reshape(zstack2_green(:,:,iPlane-1:iPlane), [a*b 2]);
    r = corr(zstack2_long(:,1),zstack2_long(:,2));
    if r<.5
        av = mean(zstack_sub_green(:,:,1:10),3);
        [out,reg] = stackRegister(zstack_sub_green,av,1);
        zstack_new = mean(reg,3);
        [out,reg] = stackRegister(zstack_new, zstack2_green(:,:,iPlane-1));
        zstack2_green(:,:,iPlane)= reg;
    end
end

fn_out_green = fullfile(outDir,sprintf('%s_green_avg_reg.tif',expt));
writetiff(zstack2_green, fn_out_green);

%% register after averaging
zstack_down = stackGroupProject(zstack_g,20);
zstack_reg = zeros(a,b,c/20);
zstack_reg(:,:,1) = zstack_down(:,:,1);
for iframe = 2:c/20
    [out reg] = stackRegister(zstack_down(:,:,iframe),zstack_down(:,:,iframe-1),1);
    zstack_reg(:,:,iframe) = reg;
end
zstack3_green = stackGroupProject(zstack_reg,frames/20);

fn_out_green = fullfile(outDir,sprintf('%s_green_down_reg.tif',expt));
writetiff(zstack3_green, fn_out_green);