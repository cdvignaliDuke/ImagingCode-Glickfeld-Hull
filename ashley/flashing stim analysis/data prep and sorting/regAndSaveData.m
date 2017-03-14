function data_reg = regAndSaveData(data,avg_frames,fnout)
data_avg = mean(data(:,:,avg_frames),3);

[out data_reg] = stackRegister(data, data_avg);

clear data
save(fullfile(fnout,'reg_img'),'data_avg','out');
end