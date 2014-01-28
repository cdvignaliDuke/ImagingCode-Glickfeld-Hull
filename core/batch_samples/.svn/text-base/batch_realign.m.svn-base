
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batch file to realign images
%   data_dir
%       |-tif   xxxx.tif (original stack)
%       |-rtif  rxxxx.tif (realigned stack)
%
% start a few slaves (startRegSlave.m) before running this batch
%
%       Kenichi Ohki  04.28.08

UseTurboReg = 0;    % 0: Use DFT, 1: Use TurboReg

if UseTurboReg
    mcoreDIR='E:\users\kenichi\multcore';
    % set this to your multicore directory (avoid to use my folder!)
end

overwrite=0; % do not overwrite

data_dir='K:\images\mouse090727';
analysis_dir='K:\images\mouse090727';


dirs=GetDirs(data_dir, analysis_dir);

[filelist, nFiles] = GetFileList(dirs, 'data');
for i=1:nFiles
    fnames=GetFnames(filelist{i},data_dir,analysis_dir);
    rfnames=GetFnames(filelist{i},data_dir,analysis_dir,'+r');
    if exist(rfnames.rdata)==2 && overwrite~=1
         continue;
    end
    if UseTurboReg
        stackRealignTifStack(fnames.data, rfnames.rdata, mcoreDIR);
    else
        stack=readtiff(fnames.data);
        target=double(sum(stack,3));
        [outs,rstack]=stackRegister(stack,target,1);
        writetiff(rstack, rfnames.rdata);
        writetiff(mean(rstack,3),rfnames.avg);
        save (rfnames.rparams, 'outs', '-v6');
    end
end
clear all
