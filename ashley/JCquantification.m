voting =    [4,1,3,5,2;...
    4,2,1,3,5;...
    1,5,4,2,3;...
    1,2,4,3,5;...
    3,4,5,2,1;...
    1,2,4,5,3;...
    2,4,1,3,5];
        
        
[nParticipants,nArticles] = size(voting);
weighting = nArticles:-1:1;

score = zeros(1,nArticles);
for i = 1:nArticles
    rank = voting(:,i);
    c = histcounts(rank,1:(nArticles+1));
    score = score + c.*weighting(i);
end

figure
bar(1:nArticles,score)
figXAxis([],'Article Number',[])
figYAxis([],'Score',[])
figAxForm
title({date;'score = sum of weighted votes'})
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\ashley\Analysis\_temp figs\JC voting\' date],'-dpdf','-fillpage')