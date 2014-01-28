%%
clear all

%%
date = '121022';
mouse = 'CM118';
userun = [2:5];
count_protocol = 1;
blanks = 1;
dirs = 1;
channel = '_green';
down = 12;
base = 'G:\users\lindsey\analysisLG\active mice\';
outDir = fullfile(base, mouse, date);

%% register (single plane) and decimate
for iRun = 1:length(userun);
    newdir = fullfile('\\zoloto\bigstorlab\Lindsey\fastrig\2012\October', [date '_' mouse], ['run' num2str(userun(iRun))],['run' num2str(userun(iRun)) channel]);
    stack = uint16(readtiff(newdir));
    if iRun == 1;
        av = mean(stack(:,:,1:50),3);
        [out, reg] = stackRegister(stack,av,10);                
        stack_down = stackGroupProject(reg, down);
        av = mean(reg(:,:,1:1000),3);
    else
        [out, reg] = stackRegister(stack,av,10);                
        stack_down = stackGroupProject(reg, down);
    end
    fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '_dec_reg.tif']);
    writetiff(stack_down, fn_out);
    clear stack_down
    clear stack
    clear reg
end

%%
PARAMS_121022_CM118
resort_seq_only
nON = 12;
nOFF = 12;
pre_win = [7 12];
post_win = [13 24];
P = 2;
nPlanes = 1;
nCond = 25;
begin = 1;

stack_sort

%% register individual runs
for iRun = 1:length(userun);
    newdir = fullfile('K:\images', [date '_' mouse], ['run' num2str(userun(iRun))],['run' num2str(userun(iRun)) channel]);
    stack = uint16(readtiff(newdir));
    av = mean(stack(:,:,1:50),3);
    [out, reg] = stackRegister(stack,av,10);
    fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '_reg.tif']);
    writetiff(reg, fn_out);
end


%% decimate before registering

for iRun = 1:length(userun);
    newdir = fullfile('K:\images', [date '_' mouse], ['run' num2str(userun(iRun))],['run' num2str(userun(iRun)) channel]);
    stack = uint16(readtiff(newdir));
    stack_down = stackGroupProject(stack, down);
    if iRun == 1;
        av = mean(stack_down(:,:,1:10),3);
        [out, reg] = stackRegister(stack_down,av,10);                
        av = mean(stack(:,:,1:500),3);
    else
        [out, reg] = stackRegister(stack_down,av,10);                
    end
    fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '_dec_reg.tif']);
    writetiff(reg, fn_out);
end
