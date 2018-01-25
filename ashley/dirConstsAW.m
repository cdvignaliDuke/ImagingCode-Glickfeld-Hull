function rc = dirConstsAW


tHostname = lower(hostname);
[s,tUsername] = dos('ECHO %USERNAME%');

switch tHostname
    case {'nuke'}      
        if tUsername(1:5) == 'linds'
            rc.name = 'linds';
            rootDir = '\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis';
        elseif tUsername(1:5) == 'ashle'
            rc.name = 'ashle';
            rc.behavData = 'Y:\home\andrew\Behavior\Data';
            rootDir = '\\CRASH.dhe.duke.edu\data\home\ashley';
            rc.dataPath = fullfile(rootDir,'data');
            rc.FSAVSummaryAnalysis = fullfile(rootDir, 'FSAV Summaries');
            rc.ashleyAnalysis = fullfile(rootDir,'Analysis');
        end
end

end

