% look at traces of individual cells across trials, color coded by
% condition
width=[min(abs(cStimOn-cDecision)):30:max(abs(cStimOn-cDecision))];
figure; histogram(abs(cStimOn-cDecision),width)
xlabel('frames (1 every 33.3 msec)')
%%
cll=48;

sr1=1;
sucright100=zeros(50,length(idx{16}));
sl1=1;
sucleft100=zeros(50,length(idx{13}));
ur1=1;
unright100=zeros(50,length(idx{17})+length(idx{18}));
ul1=1;
unleft100=zeros(50,length(idx{14})+length(idx{15}));
sr6=1;
sucright6=zeros(50,length(idx{4}));
sl6=1;
sucleft6=zeros(50,length(idx{1}));
sr10=1;
sucright10=zeros(50,length(idx{10}));
sl10=1;
sucleft10=zeros(50,length(idx{7}));
ur10=1;
unright10=zeros(50,length(idx{11})+length(idx{12}));
ul10=1;
unleft10=zeros(50,length(idx{8})+length(idx{9}));
ur6=1;
unright6=zeros(50,length(idx{5})+length(idx{6}));
ul6=1;
unleft6=zeros(50,length(idx{2})+length(idx{3}));


for ii=1:nTrials
    if and(and(SIx(ii),~tLeftTrial(ii)),tGratingContrast(ii)==1)
        sucright100(:,sr1)=data_stim_dfof(:,cll,ii);
        sr1=sr1+1;
    end
end

for ii=1:nTrials
    if and(and(SIx(ii),tLeftTrial(ii)),tGratingContrast(ii)==1)
        sucleft100(:,sl1)=data_stim_dfof(:,cll,ii);
        sl1=sl1+1;
    end
end

for ii=1:nTrials
    if and(and(~SIx(ii),~tLeftTrial(ii)),tGratingContrast(ii)==1)
        unright100(:,ur1)=data_stim_dfof(:,cll,ii);
        ur1=ur1+1;
    end
end

for ii=1:nTrials
    if and(and(~SIx(ii),tLeftTrial(ii)),tGratingContrast(ii)==1)
        unleft100(:,ul1)=data_stim_dfof(:,cll,ii);
        ul1=ul1+1;
    end
end

for ii=1:nTrials
    if and(and(SIx(ii),~tLeftTrial(ii)),tGratingContrast(ii)==cons(1))
        sucright6(:,sr6)=data_stim_dfof(:,cll,ii);
        sr6=sr6+1;
    end
end

for ii=1:nTrials
    if and(and(SIx(ii),tLeftTrial(ii)),tGratingContrast(ii)==cons(1))
        sucleft6(:,sl6)=data_stim_dfof(:,cll,ii);
        sl6=sl6+1;
    end
end

for ii=1:nTrials
    if and(and(SIx(ii),~tLeftTrial(ii)),tGratingContrast(ii)==cons(2))
        sucright10(:,sr10)=data_stim_dfof(:,cll,ii);
        sr10=sr10+1;
    end
end

for ii=1:nTrials
    if and(and(SIx(ii),tLeftTrial(ii)),tGratingContrast(ii)==cons(2))
        sucleft10(:,sl10)=data_stim_dfof(:,cll,ii);
        sl10=sl10+1;
    end
end

%

for ii=1:nTrials
    if and(and(~SIx(ii),~tLeftTrial(ii)),tGratingContrast(ii)==cons(2))
        unright10(:,ur10)=data_stim_dfof(:,cll,ii);
        ur10=ur10+1;
    end
end

for ii=1:nTrials
    if and(and(~SIx(ii),tLeftTrial(ii)),tGratingContrast(ii)==cons(2))
        unleft10(:,ul10)=data_stim_dfof(:,cll,ii);
        ul10=ul10+1;
    end
end

for ii=1:nTrials
    if and(and(~SIx(ii),~tLeftTrial(ii)),tGratingContrast(ii)==cons(1))
        unright6(:,ur6)=data_stim_dfof(:,cll,ii);
        ur6=ur6+1;
    end
end

for ii=1:nTrials
    if and(and(~SIx(ii),tLeftTrial(ii)),tGratingContrast(ii)==cons(1))
        unleft6(:,ul6)=data_stim_dfof(:,cll,ii);
        ul6=ul6+1;
    end
end


%% pick and choose
figure
hold on
for ii=1:size(sucright100,2)
    plot(tt,sucright100(:,ii),'Color',[.8,.8,1])
end

for ii=1:size(sucleft100,2)
    plot(tt,sucleft100(:,ii),'Color',[1,.8,.8])
end

for ii=1:size(unright100,2)
    plot(tt,unright100(:,ii),'Color',[.7,.9,.7])
end

for ii=1:size(unleft100,2)
    plot(tt,unleft100(:,ii),'Color',[.8,.8,.8])
end

for ii=1:size(sucright6,2)
    plot(tt,sucright6(:,ii),'Color',[.7 .9 .9])
end

for ii=1:size(sucleft6,2)
    plot(tt,sucleft6(:,ii),'Color',[.9 .8 .9])
end

for ii=1:size(sucright10,2)
    plot(tt,sucright10(:,ii),'Color',[.9 .9 .7])
end

for ii=1:size(sucleft10,2)
    plot(tt,sucleft10(:,ii),'Color',[.7 .8 .9])
end
%
for ii=1:size(unright6,2)
    plot(tt,unright6(:,ii),'Color',[.9 .7 .9])
end

for ii=1:size(unleft6,2)
    plot(tt,unleft6(:,ii),'Color',[.8 .9 .9])
end

for ii=1:size(unright10,2)
    plot(tt,unright10(:,ii),'Color',[.9 .9 .7])
end

for ii=1:size(unleft10,2)
    plot(tt,unleft10(:,ii),'Color',[.7 .8 .9])
end

plot(tt,mean(sucleft100,2),'Color',[1,0,0])
plot(tt,mean(sucright100,2),'Color',[0,0,1])
plot(tt,mean(unright100,2),'Color',[0,.5,0])
plot(tt,mean(unleft100,2),'Color',[0,0,0])
plot(tt,mean(sucright6,2),'Color',[0 .5 .5])
plot(tt,mean(sucleft6,2),'Color',[1 0 1])
plot(tt,mean(sucright10,2),'Color',[.5 .5 0])
plot(tt,mean(sucleft10,2),'Color',[.5 .5 1])
%
plot(tt,mean(unright6,2),'Color',[ .5 0 .5])
plot(tt,mean(unleft6,2),'Color',[0 1 1])
plot(tt,mean(unright10,2),'Color',[.5 .5 0])
plot(tt,mean(unleft10,2),'Color',[.5 .5 1])
%% early/late
sr1=1;
for ii=1:nTrials/2
    if and(and(SIx(ii),~tLeftTrial(ii)),tGratingContrast(ii)==1)
        tt1(:,sr1)=data_stim_dfof(:,cll,ii);
        sr1=sr1+1;
    end
end
sr1=1;
for ii=nTrials/2:nTrials
    if and(and(SIx(ii),~tLeftTrial(ii)),tGratingContrast(ii)==1)
        tt2(:,sr1)=data_stim_dfof(:,cll,ii);
        sr1=sr1+1;
    end
end

figure
hold on
for ii=1:size(tt1,2)
    plot(tt,tt1(:,ii),'Color',[.8,.8,1])
end
for ii=1:size(tt2,2)
    plot(tt,tt2(:,ii),'Color',[.8,.8,.8])
end  
plot(tt,mean(tt1,2),'Color',[0,0,1])
plot(tt,mean(tt2,2),'Color',[0,0,0])
    
    
    
%% labels
xlabel('ms')
ylabel('df/f')
title('success v fail') %variable