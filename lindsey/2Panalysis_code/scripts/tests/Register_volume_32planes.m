clear all
%%
date = '120913';
mouse = 'VC2';
userun = [1:3];
count_protocol = 1;
blanks = 1;
dirs = 1;
channel = '_green';
down = 32;
base = 'G:\users\lindsey\analysisLG\active mice\';
outDir = fullfile(base, mouse, date);

%% register (single plane) and decimate
for iRun =3:length(userun);
    newdir = fullfile('\\zoloto\bigstorlab\Lindsey\fastrig\2012\September', [date '_' mouse], ['run' num2str(userun(iRun))],['run' num2str(userun(iRun)) channel]);
    stack = uint16(readtiff(newdir));
    if iRun == 1;
        stack_down = stackGroupProject(stack, down);
        av = mean(stack_down(:,:,1:10),3);
        [out, reg] = stackRegister(stack_down,av,10);                
        av = mean(reg(:,:,1:500),3);
    else
        stack_down = stackGroupProject(stack, down);
        [out, reg] = stackRegister(stack_down,av,10);                        
    end
    fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '_dec_reg.tif']);
    writetiff(reg, fn_out);
end

%%
nPlanes = 32;
siz = size(stack);
vol_stack = zeros(siz(1),siz(2),nPlanes);
for iplane = 1:nPlanes 
    plane_stack = stack(:,:,iplane:nPlanes:nPlanes*50);
    av = mean(plane_stack(:,:,1:50),3);
    [out, reg] = stackRegister(plane_stack, av, 10);
    plane_av = mean(reg,3);
    vol_stack(:,:,iplane) = plane_av;
end

fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun) '_volume.tif']);
    writetiff(vol_stack, fn_out);
    
%%
nON = 5;
nOFF = 5;
pre_win = [3 5];
post_win = [6 10];
P = 2;
nPlanes = 1;
nCond = 25;
begin = 1;
PARAMS_120913_VC2
resort_seq_only

stack_sort