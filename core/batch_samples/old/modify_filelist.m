filelist='exp_list_scmr5.xls';
s=readXLSfilelist(filelist);
for i=1:length(s)
    s(i).filename=['r',s(i).filename,'_ch1'];
end
writeXLSfilelist(s,'E:\users\kenichi\filelist\exp_list_scmr5_2.xls');
