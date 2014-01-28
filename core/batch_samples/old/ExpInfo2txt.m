function [expinfo, stiminfo, scaninfo]=ExpInfo2txt(filelist,MPparams) 
expinfo=[];
stiminfo=[];
scaninfo=[];

%%%%%%%%%%%%%%
% expinfo
names{1}=    'run';
names{2}=    'expdir';
names{3}=    'filename';
names{4}=    'depth';
names{5}=    'stim_description';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stiminfo
names{6}=    'stim_code';
names{7}=    'stim_type';
names{8}=    'eye';
names{9}=    'Noff';
names{10}=    'Non';
names{11}=    'Nstim_per_run';
names{12}=    'Nrep';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% scaninfo

names{13}=    'zoom';
names{14}=   'MATsize';
names{15}=    'Nframes';
names{16}=    'framerate';
names{17}=    'truncated';

MPnames{1}='Magnification';
MPnames{2}='FrameHeight';
MPnames{3}='FrameWidth';
MPnames{4}='FrameCount';
MPnames{5}='Rotation';
MPnames{6}='Enabled1';
MPnames{7}='Enabled2';
MPnames{8}='XPosition';
MPnames{9}='YPosition';

for i=1:5
    temp=getfield(filelist,names{i});
    if isnumeric(temp)
        temp=num2str(temp);
    end
    expinfo=[expinfo,' ',names{i},': ',temp,', '];
end

for i=6:12
    temp=getfield(filelist,names{i});
    if isnumeric(temp)
        temp=num2str(temp);
    end
    stiminfo=[stiminfo,' ',names{i},': ',temp,', '];
end

for i=13:17
    temp=getfield(filelist,names{i});
    if isnumeric(temp)
        temp=num2str(temp);
    end
    scaninfo=[scaninfo,' ',names{i},': ',temp,', '];
end

for i=1:9
    temp=getfield(MPparams,MPnames{i});
    if isnumeric(temp)
        temp=num2str(temp);
    end
    scaninfo=[scaninfo,' ',MPnames{i},': ',temp,', '];
end

expinfo(strfind(expinfo,'_'))='-';
stiminfo(strfind(stiminfo,'_'))='-';
scaninfo(strfind(scaninfo,'_'))='-';

