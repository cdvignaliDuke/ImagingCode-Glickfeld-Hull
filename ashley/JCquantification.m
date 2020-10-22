ID = {'Miao','Grace','Ashley','Nev','Jenny'};

voting =    [1,3,5,2,4;...
    1,3,2,4,5;...
    2,1,3,4,5;...
    5,2,3,1,4;...
   3,4,5,1,2;...
   2,3,4,1,5];
        
        
[nParticipants,nArticles] = size(voting);
weighting = nArticles:-1:1;

score = zeros(1,nArticles);
for i = 1:nArticles
    rank = voting(:,i);
    c = histcounts(rank,1:(nArticles+1));
    score = score + c.*weighting(i);
end

wkday = weekday(datenum(date));
tuesdayDate = datenum(date);
while wkday~=3
    wkday = weekday(tuesdayDate+1);
    tuesdayDate = tuesdayDate+1;
end

figure
bar(1:nArticles,score)
% figXAxis([],'Article Number',[0 length(ID)+1],1:length(ID),ID)
figXAxis([],'Article Number',[0 max(voting(1,:))+1])
figYAxis([],'Score',[])
figAxForm
title({['Tues. JC ' datestr(tuesdayDate)];'score = sum of weighted votes'})
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\ashley\Analysis\_temp figs\JC voting\' datestr(tuesdayDate)],'-dpdf','-fillpage')