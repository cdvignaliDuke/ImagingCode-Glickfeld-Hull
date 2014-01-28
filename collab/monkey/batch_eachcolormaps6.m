
filelist='exp_list_scmr6_2.xls';
s=readXLSfilelist(filelist);

outDirectory = 'E:\users\kenichi\eachcolormaps\scmr6\';
InDirectory = 'E:\users\kenichi\scmr6\pixelmaps\'

for i=1:length(s)
    if (isempty(strfind(s(i).stim_code,'Img_flash')) & isempty(strfind(s(i).stim_code,'sweepbar')))
        continue;
    end
    load(fullfile(InDirectory,s(i).filename,'dir_dF_sm.mat'));
    [R,G,B,RG]=getEachColormap2(dir_dF_sm);
    imwrite(R, fullfile(outDirectory, [s(i).filename, '_L.bmp']));
    imwrite(G, fullfile(outDirectory, [s(i).filename, '_M.bmp']));
    imwrite(B, fullfile(outDirectory, [s(i).filename, '_S.bmp']));
    imwrite(RG, fullfile(outDirectory, [s(i).filename, '_LM.bmp']));
end
