% dataset directory = FSAV_V1_100ms_naive_temp.m
% 
% Protocol:
% 1) run compareRegImg to choose an average image for registering the whole
% dataset
% 2) run reRegDatasets.m in sections - first run up to get getTaskMaxDFF 
% (126) and dirTuningMaxDFF (128). Then run the next sections one at a time. 
% You will first crop out the bright edges, then select cells from all of the
% max dF/F images, then you will extract and save timecourses and neuropil.
