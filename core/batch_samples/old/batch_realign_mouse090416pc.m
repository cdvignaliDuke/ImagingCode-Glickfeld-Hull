
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

tif_dir=[data_dir,'tifpc\'];
rtif_dir=[data_dir,'rtifpc\'];
mkdir(data_dir, 'rtifpc');
[filelist, nFiles] = GetTifFileNames(tif_dir);
for i=1:nFiles
    [pathstr, name, ext, versn] = fileparts(filelist{i});
    fname=fullfile(tif_dir, filelist{i});
    rfname=['r',filelist{i}];
    outfname=fullfile(rtif_dir, rfname);
    outmatname=fullfile(rtif_dir, ['r',name,'.mat']);
    avgname=fullfile(rtif_dir, ['r',name,'_avg.tif']);
    stack=readtiff(fname);
    if exist(outfname)==2 && overwrite~=1
         continue;
    end
    target=double(sum(stack,3));
    [outs,rstack]=stackRegister(stack,target,1);
    writetiff(rstack, outfname);
    writetiff(mean(rstack,3),avgname);
    save (outmatname, 'outs', '-v6');
end
clear all