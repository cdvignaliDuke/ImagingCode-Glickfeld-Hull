function [green_reg, red_reg, shifts] = stackRegisterMulticore(greenStack,redStack, multicoredir, nChunks)
%STACKREGISTERMULTICORE Starts multicore registration
% [green_reg] = stackRegisterMulticore(greenStack,[], multicoredir, nchunks) 
%   Registers stack based on data from one channel. Fast.
% [green_reg, red_reg] = stackRegisterMulticore(greenStack, redStack, multicoredir, nchunks)
%   Registers stack based on data from sum of two channels. Slow.
%
% see stackRegisterJava for single core version of this code

nchan = 1;
if ~isempty(redStack);nchan = 2;end

[ny,nx,nframes]=size(greenStack);

list = dir(fullfile(multicoredir,'*.mat'));

if length(list) > 2
    disp(multicoredir);
    prompt = sprintf('Directory not empty. Do you want delete all files?\n');
    str = input(prompt,'s');
    if strcmp(str,'y')
        delete(fullfile(multicoredir,'*.*'));
        fprintf(1,'Deleted all mat files from directory %s\n',multicoredir);
    else
        fprintf(1,'Exited\n',multicoredir);
        return;
    end
end

% target is average 128 frames from middle of the stack
siz = 128;
sel = floor(nframes/2):floor(nframes/2)+siz-1;

switch nchan 
    case 1 % register from green channel only
        target = mean(greenStack(:,:,sel),3);               
        green_reg = RegMulticore(greenStack, target, [], multicoredir, nChunks);
        
    case 2 % register using green and red channels
        target = mean(greenStack(:,:,sel)+redStack(:,:,sel),3);       
        shifts = RegMulticore_shifts(greenStack+redStack, target, [], multicoredir, nChunks);

        green_reg = stackShifts(double(greenStack), shifts);
        red_reg = stackShifts(double(redStack), shifts);

        green_reg(find(green_reg(:)==nan)) = 0;
        red_reg(find(red_reg(:)==nan)) = 0;
        
    otherwise
        error('something went wrong');
end    

return;

