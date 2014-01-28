function startRegSlave 
% startRegSlave

dirs = frGetDirs;

fprintf('ImageJ Version %s\n',char(ij.IJ.getVersion));
t = timer('TimerFcn',@ijGetFreeMemory, 'Period', 60,'ExecutionMode','FixedRate');
start(t);
fprintf('Starting multicore slave from %s\n',dirs.multicore);
startmulticoreslave(dirs.multicore);

return;
