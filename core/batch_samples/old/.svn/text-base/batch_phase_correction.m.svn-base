
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batch file to realign images
%   data_dir
%       |-tif   xxxx.tif (original stack)
%       |-rtif  rxxxx.tif (realigned stack)
%
% start a few slaves (startRegSlave.m) before running this batch
%
%       Kenichi Ohki  04.28.08

% set this to your multicore directory (avoid to use my folder!)

overwrite=0; % do not overwrite

data_dir='E:\users\kenichi\leicatemp\';

tif_dir=[data_dir,'tif\'];
tifpc_dir=[data_dir,'tifpc\'];
mkdir(data_dir, 'tifpc');
[filelist, nFiles] = GetTifFileNames(tif_dir);
for i=1:nFiles
    [pathstr, name, ext, versn] = fileparts(filelist{i});
    fname=fullfile(tif_dir, filelist{i});
    pcfname=[filelist{i},'pc'];
    outfname=fullfile(tifpc_dir, pcfname);
    outmatname=fullfile(tifpc_dir, [name,'pc_shift.mat']);
    avgname=fullfile(tifpc_dir, [name,'pc_avg.tif']);
    stack=readtiff(fname);
    if exist(outfname)==2 && overwrite~=1
         continue;
    end
    [stack2, shift]=stackBidirShiftCorrection(stack);
    writetiff(stack2, outfname);
    writetiff(mean(stack2,3),avgname);
    save (outmatname, 'shift', '-v6');
end
clear all