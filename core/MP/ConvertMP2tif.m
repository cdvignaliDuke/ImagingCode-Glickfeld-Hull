function params=ConvertMP2tif (fname, MPpath, tifpath, overwrite);

if nargin < 4
    overwrite=0;
end

mpfile=OpenMPfile(fullfile(MPpath, fname));
params=ReadMPparams(mpfile)

[dummypath, name, ext] = fileparts(fname) ;

channel=1;
if strcmpi(params.Enabled1,'True')
    name1=[name,'_ch',int2str(channel)];
    ext='.tif';
    outtifname=fullfile(tifpath, [name1, ext]);
    ext='.mat';
    outmatname=fullfile(tifpath, [name1, ext]);
    if overwrite==1 & exist(outtifname)==2
        delete(outtifname);
    end
    if exist(outtifname)~=2
        for i=1:params.FrameCount
            temp=uint16(ReadMPframe(mpfile, i, channel, params));
            imwrite(temp,outtifname,'tif','Compression','none','WriteMode','append');
        end
    end
    save (outmatname, 'params', '-v6');
end

channel=2;
if strcmpi(params.Enabled2,'True')
    name2=[name,'_ch',int2str(channel)];
    ext='.tif';
    outtifname=fullfile(tifpath, [name2, ext]);
    ext='.mat';
    outmatname=fullfile(tifpath, [name2, ext]);
    if overwrite==1 & exist(outtifname)==2
        delete(outtifname);
    end
    if exist(outtifname)~=2
        for i=1:params.FrameCount
            temp=uint16(ReadMPframe(mpfile, i, channel, params));
            imwrite(temp,outtifname,'tif','Compression','none','WriteMode','append');
        end
    end
    save (outmatname, 'params', '-v6');
end

close all;
