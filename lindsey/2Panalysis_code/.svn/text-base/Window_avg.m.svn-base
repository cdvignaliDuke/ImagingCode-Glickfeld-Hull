pn = 'G:\users\lindsey\dataLG\110327_L5\run4.tif';
pn_out = 'G:\users\lindsey\dataLG\110327_L5\';
Nstimtypes = 6;
Non = 10;
Noff = 10;
pre_win =[16:20]; 
post_win = [6:10];

stack = readtiff(pn);
[x y nframes] = size(stack);
Ntot = Non + Noff;
cycle = rate*(Non+Noff)*Nstimtypes;
repeat = nframes/cycle;
stack_pre = zeros(x,y,Nstimtypes);
[a,b] = size(pre_win);
ind_pre = zeros(1,(repeat*b));
for iN = 1:Nstimtypes;
   for iR = 1:repeat
        ind_pre(1,1+((iR-1)*b):b+((iR-1)*b))=(pre_win(1):pre_win(end))+((iR-1)*cycle+(iN-1)*Ntot);
   end
    stack_pre(:,:,iN)= mean(stack(:,:,ind_pre),3);
end

stack_post = zeros(x,y,Nstimtypes);
[a,b] = size(post_win);
ind_post = zeros(1,(repeat*b));
for iN = 1:Nstimtypes;
   for iR = 1:repeat
        ind_post(1,1+((iR-1)*b):b+((iR-1)*b))=(post_win(1):post_win(end))+((iR-1)*cycle+(iN-1)*Ntot);
   end
    stack_post(:,:,iN)= mean(stack(:,:,ind_post),3);
end

stack_out = stack_post./stack_pre;
avg = (mean(stack,3));
stack_avg = zeros(x,y,Nstimtypes);
for iN = 1:Nstimtypes;
    stack_avg(:,:,iN) = avg;
end
stack_post_dF = stack_post./stack_avg;


fn_out = fullfile(pn_out,sprintf('run4_post_dF'));
writetiff(stack_post_dF, fn_out);
fn_out = fullfile(pn_out,sprintf('run4_dF'));
writetiff(stack_out, fn_out);

