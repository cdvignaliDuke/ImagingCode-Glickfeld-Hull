%% expt params
date = '110509';
mouse = 'Y13';

userun = [1:2];
count_protocol = 1;

blanks = 1;
run = 0;

nCond = 36;
%%
P = 2;
nON = 12;
nOFF = 12;
nPlanes = 1;
begin = 7;
TFSFetc = [1:2];
pre_win = [1 6];
post_win = [7 18];

base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse,date);
zoloto = '\\zoloto\bigstorlab\Lindsey\fastrig\2011\May';

%% load and decimate original stack
fn_in = fullfile(zoloto, [date '_' mouse], ['run' num2str(userun(2))], ['run' num2str(userun(iRun)) '_green']);
stack = readtiff(fn_in);
stack_dec = stackGroupProject(stack,12);
clear stack

%fix lag on decimated stack
[outStack,lagP] = stackFixBidirPhase(stack_dec);
clear outStack
clear stack_dec

lagP_interp = round(interp(lagP,12));

%fix lag on original stack
list = dir(fullfile(fn_in,'*.tif')); 
nfiles = length(list);
start = 1;
for ifile = 1:nfiles;
    stack = readtiff(fn_in, ifile);
    z = size(stack,3);
    stack_shiftcorr = imBidirShift(stack, lagP_interp(:,start:start+z-1));
    fn_out = fullfile(outDir, 'shiftcorr', [date '_' mouse '_run' num2str(userun(iRun))],[date '_' mouse '_run' num2str(userun(iRun)) sprintf('_%06d.tif',ifile)]);
    writetiff(stack_shiftcorr, fn_out);
    start = start+z;
end

%decimate lag-shifted stack
fn_corr = fullfile(outDir, 'shiftcorr', [date '_' mouse '_run' num2str(userun(iRun))]);
stack_shiftcorr = readtiff(fn_corr);
stack_shiftcorr_dec = stackGroupProject(stack_shiftcorr,12);
fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun(2)) '_stack_shiftcorr_dec.tif']);
writetiff(stack_shiftcorr_dec, fn_out);
clear stack_shiftcorr

%register decimated stack
av = mean(stack_shiftcorr_dec(:,:,115:215),3);
[out reg] = stackRegister(stack,av,10);
fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '_dec_reg.tif']);
writetiff(reg, fn_out);
fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '_reg_out.mat']);
save(fn_out, 'out');
clear reg
clear stack_shiftcorr_dec

%interpolate shifts from registration
out_interp_row = round(interp1(9:18:(size(out,1)*18), out(:,3),1:1:(size(out,1)*18))*10)/10;
out_interp_col = round(interp1(9:18:(size(out,1)*18), out(:,4),1:1:(size(out,1)*18))*10)/10;
out_interp = cat(1,out_interp_col, out_interp_row)';
for iframe = 1:8
    out_interp(iframe,:) = out_interp(9,:);
end
for iframe = size(out_interp,1)-9:size(out_interp,1)
    out_interp(iframe,:) = out_interp(size(out_interp,1)-10,:);
end

figure;
plot(out_interp(1:1800,1), 'k-*');
hold on; 
plot(9:18:1800, out(1:100,4),'r-^');

%apply to original stack
fn_corr = fullfile(outDir, 'shiftcorr', [date '_' mouse '_run' num2str(userun(iRun))]);
list = dir(fullfile(fn_corr,'*.tif')); 
nfiles = length(list);
start = 0;
for ifile = 1:nfiles;
    stack = readtiff(fn_corr,ifile);
    sz = size(stack);
    X0 = 1:sz(1);
    Y0 = 1:sz(2);
    [X00,Y00] = meshgrid(X0,Y0);
    X = X00';
    Y = Y00';
    stack_reg = zeros(size(stack));
    for iframe = 1:sz(3);
        X2 = X-(out_interp(iframe+start,2));
        Y2 = Y-(out_interp(iframe+start,1));
        tmp = squeeze(double(stack(:,:,iframe)));
        stack_reg(:,:,iframe) = interp2(Y,X,tmp,Y2,X2,'cubic',0);
    end
    start = start+sz(3);
    fn_out = fullfile(outDir, 'reg', [date '_' mouse '_run' num2str(userun(2))],[date '_' mouse '_run' num2str(userun(2)) sprintf('_%06d_cubic.tif',ifile)]);
    writetiff(stack_reg, fn_out);
end