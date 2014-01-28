
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batch file to realign images
%   data_dir
%       |-tif   xxxx.tif (original stack)
%       |-rtif  rxxxx.tif (realigned stack)
%
% start a few slaves (startRegSlave.m) before running this batch
%
%       Kenichi Ohki  04.28.08

mcoreDIR='E:\users\kenichi\multcore';
% set this to your multicore directory (avoid to use my folder!)

overwrite=0; % do not overwrite

data_dir='\\zmey\storlab\data\Kenichi\mouse090205\';

tif_dir=[data_dir,'tif\'];
rtif_dir=[data_dir,'rtif\'];
mkdir(data_dir, 'rtif');
[filelist, nFiles] = GetTifFileNames(tif_dir);
for i=1:nFiles
    fname=fullfile(tif_dir, filelist{i});
    rfname=['r',filelist{i}];
    outfname=fullfile(rtif_dir, rfname);
    if exist(outfname)==2 && overwrite~=1
         continue;
    end
    stackRealignTifStack(fname, outfname, mcoreDIR)
end
clear all