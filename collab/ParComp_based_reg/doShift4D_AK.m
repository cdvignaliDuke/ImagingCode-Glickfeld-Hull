function RegInfo=doShift4D(RegInfo,trial_per_thread)

datestr(clock)

path = RegInfo.All.path;
fbase = RegInfo.All.fbase;
ch = RegInfo.All.ch;
fnums = RegInfo.All.fnums;
ftarget = RegInfo.All.ftarget;
is_rect = RegInfo.All.is_rect;
reverse_ch = RegInfo.All.reverse_ch;
nPlanes = RegInfo.All.nPlanes;
Nvolbin = 1;
Nupsamp = RegInfo.All.Nupsamp;
Nfromedge = RegInfo.All.Nfromedge; 


%ch = '__ch1.tif';
%Nvolbin = 16;
%Nupsamp = 4; %no subpixel = 1, subpixel -> integer>1
%Nfromedge = [6 6 6]; %distance in pixels for potential correg


% 
% fstart=strfind(ftarget,'\');
% ftarget2=ftarget(fstart(end):end);

target = readtiff(path,[],ftarget);

%Align and save each run
for v=1:length(fbase)
%return
display(['loading original stack ',num2str(v),' of ', num2str(length(fbase))]);
fid=fopen([path,fbase{v},'_fix.bin'],'r');
sz=RegInfo.stack_fix_sz{v};
stack=fread(fid,[sz(1)*sz(2),sz(3)],'uint16=>uint16');
fclose(fid);
display('done loading original stack');

stack=reshape(stack,sz(1),sz(2),sz(3));

nSlices = size(stack,3);
%nVol=floor(nSlices/nPlanes);
nTrial=floor(nSlices/(Nvolbin*nPlanes));
nChunks = ceil(nTrial/trial_per_thread);

parameterCell = cell(1,nChunks);

for i=1:nChunks
        start = (i-1)*(trial_per_thread*nPlanes*Nvolbin)+1;
        if i == nChunks
            stop = nSlices;
        else
            stop = start + (trial_per_thread*nPlanes*Nvolbin) - 1;
        end
        substack=stack(:,:,start:stop);
        parameterCell{1,i} = {substack,target,nPlanes,is_rect,reverse_ch,Nupsamp,Nfromedge,Nvolbin};
end

clear stack;
display('starting multicore shiftreg')

parfor i=1:size(parameterCell,2)
     parameters=parameterCell{1,i};
     out{i} = Shift4D(parameters{1}, parameters{2}, parameters{3},parameters{4},...
         parameters{5},parameters{6},parameters{7},parameters{8})
end  


ShiftMat=cell2mat(out');

RegInfo.Runs(v).ShiftMat=ShiftMat;
RegInfo.Runs(v).RunID=fnums(v);

disp(['Run',num2str(fnums(v)),' alignment complete.']);

%save([PWD,fbase,num2str(fnums(v)),'_RegInfo.mat'],'ShiftMat');

disp(['Run',num2str(fnums(v)),' output complete.']);
end

counter=0;
for v=1:length(fnums)
    si=size(RegInfo.Runs(v).ShiftMat);
    ShiftMat((counter+1):(counter+si(1)),:)=RegInfo.Runs(v).ShiftMat;
    counter=counter+si(1);
end

AllRunsStr=[];
for i=1:length(fnums)
    AllRunsStr=[AllRunsStr,num2str(fnums(i)),'_'];
end

alloutfname=[path,fbase{1},AllRunsStr];

RegInfo.All.ShiftMat=ShiftMat;
RegInfo.All.target=target;
RegInfo.All.nPlanes=nPlanes;
RegInfo.All.Nvolbin=Nvolbin;
RegInfo.All.is_rect=is_rect;
RegInfo.All.reverse_ch=reverse_ch;
RegInfo.All.path=path;
RegInfo.All.fbase=fbase;
RegInfo.All.ch=ch;
RegInfo.All.AllRunsStr=AllRunsStr;
RegInfo.All.alloutfname=alloutfname;

save([alloutfname,'RegInfo.mat'],'RegInfo');

datestr(clock)