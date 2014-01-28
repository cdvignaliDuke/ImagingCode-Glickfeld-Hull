
filelist='exp_list_scmr6_2.xls';
s=readXLSfilelist(filelist);

outDirectory = 'E:\users\kenichi\eachcolormaps\scmr6_gray\';
InDirectory = 'E:\users\kenichi\scmr6\pixelmaps\'

for i=1:length(s)
    if (isempty(strfind(s(i).stim_code,'Img_flash')) & isempty(strfind(s(i).stim_code,'sweepbar')))
        continue;
    end
    load(fullfile(InDirectory,s(i).filename,'dir_dF_sm.mat'));
    [Ron,Roff,Gon,Goff,Bon,Boff,RonGoff,GonRoff]=getEachColormap3(dir_dF_sm);
    imwrite(Ron, fullfile(outDirectory, [s(i).filename, '_Lon.bmp']));
    imwrite(Roff, fullfile(outDirectory, [s(i).filename, '_Loff.bmp']));
    imwrite(Gon, fullfile(outDirectory, [s(i).filename, '_Mon.bmp']));
    imwrite(Goff, fullfile(outDirectory, [s(i).filename, '_Moff.bmp']));
    imwrite(Bon, fullfile(outDirectory, [s(i).filename, '_Son.bmp']));
    imwrite(Boff, fullfile(outDirectory, [s(i).filename, '_Soff.bmp']));
    imwrite(RonGoff, fullfile(outDirectory, [s(i).filename, '_LonMoff.bmp']));
    imwrite(GonRoff, fullfile(outDirectory, [s(i).filename, '_MonLoff.bmp']));

end
