
function stackRealignTifStack (fname, outfname, mcoreDIR, target)

% fname: name of original stack (including pathname)
% outfname: name of realigned stack (including pathname)
% mcoreDIR: temporary folder 
% target: target tif image. If it is not provided, average image is
% obtained from the stack.
% start a few slaves (startRegSlave.m) before running this function

nAvgForTarget = 32;

%%%%%% file names

[pathstr, name, ext] = fileparts(outfname) ;
target_outfname=fullfile(pathstr, [name,'_target', ext]);
avg_outfname=fullfile(pathstr, [name,'_avg', ext]);

%%%%%%% read data

info=imfinfo(fname);
Nframes=length(info)

if Nframes==1
    array = readtiff(fname);   
    writetiff(array,target_outfname);
    writetiff(array,outfname);
    writetiff(array,avg_outfname);
     return;
end

if Nframes==0
    return;
end

array = readtiff(fname);      

%%%%%%%% registration
nChunks=ceil(Nframes./100);

if nargin<4
    Navg=min(nAvgForTarget,Nframes);
    target = mean(array(:,:,1:Navg),3);
end

registered = RegMulticore(array, target,[], mcoreDIR, nChunks);
registered_avg=mean(registered,3);

%%%%%%%% write files

writetiff(target,target_outfname);
writetiff(uint16(registered),outfname);
writetiff(registered_avg,avg_outfname);
