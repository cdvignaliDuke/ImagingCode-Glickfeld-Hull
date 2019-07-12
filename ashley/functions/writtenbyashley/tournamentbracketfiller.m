rng(457)
nTeams = 64;
nSubBrackets = 4;
%%
nTeamsPerSubBracket = nTeams/nSubBrackets;
subBracketWeights = 0.9:-(0.8/(nTeamsPerSubBracket-1)):0.1;

nSubBracketRounds = log(nTeamsPerSubBracket)./log(2);

nRound1Matchups = (nTeamsPerSubBracket/2);

roundMatchups = cell(nSubBrackets,nSubBracketRounds);
subBracketWinners = nan(1,nSubBrackets);
for ibrack = 1:nSubBrackets
    for iround = 1:nSubBracketRounds
        if iround == 1
            nMatchups = nRound1Matchups;
            roundMatchups{ibrack,iround} = cat(1,1:nRound1Matchups,...
                nTeamsPerSubBracket:-1:(nRound1Matchups+1));
        else
            nMatchups = nMatchups/2;
            roundMatchups{ibrack,iround} = cat(1,winners(1:nMatchups),...
                winners((nMatchups*2):-1:(nMatchups+1)));
        end
        winners = nan(1,nMatchups);
        for igame = 1:nMatchups
            thisGame = roundMatchups{ibrack,iround}(:,igame);
            w = subBracketWeights(thisGame);
            bucket = cat(2,ones(1,round(w(1)*1000)),ones(1,round(w(2)*1000)).*2);
            draw = bucket(randsample(length(bucket),1));
            winners(igame) = thisGame(draw);
        end
        if iround == nSubBracketRounds
            subBracketWinners(ibrack) = winners;
        end
    end
end

nSubBracketRounds = log(nSubBrackets)./log(2);
lateRoundMatchups = cell(1,nSubBracketRounds);
lateRoundDraws = cell(1,nSubBracketRounds);
for iround = 1:nSubBracketRounds
    if iround == 1
        nMatchups = length(subBracketWinners)./2;
        lateRoundMatchups{1,iround} = cat(1,...
            subBracketWinners(1:nMatchups),...
            subBracketWinners(length(subBracketWinners):-1:(nMatchups+1)));
    else
        nMatchups = nMatchups/2;
        lateRoundMatchups{1,iround} = cat(1,winners(1:nMatchups),...
            winners((nMatchups*2):-1:(nMatchups+1)));
    end
    winners = nan(1,nMatchups);
    d = nan(1,nMatchups);
    for igame = 1:nMatchups
        thisGame = lateRoundMatchups{1,iround}(:,igame);
        w = subBracketWeights(thisGame);
        bucket = cat(2,ones(1,round(w(1)*1000)),ones(1,round(w(2)*1000)).*2);
        draw = bucket(randsample(length(bucket),1));
        winners(igame) = thisGame(draw);
        d(igame) = draw;
    end
    lateRoundDraws{iround} = d;
    if iround == nSubBracketRounds
        champion = winners;
    end
end

%% count points per round
nRounds = 6;
nPointsPerRound = 320;
nCorrectPerRound_all = [25,10,3,0,0,0];

currentRound = 5;
nCorrectPerRound_current = [25,10,3,0,0,0];

bestPossibleScore = 0;
currentScore = 0;
for iround = 1:nRounds
    if iround == 1
        nResults = nTeams/2;
    else
        nResults = nResults/2;
    end
    pointsThisRound = nPointsPerRound./nResults;
    bestPossibleScore = bestPossibleScore + (pointsThisRound*nCorrectPerRound_all(iround));
    if iround <= currentRound
        currentScore = currentScore + (pointsThisRound*nCorrectPerRound_current(iround));
    end
end