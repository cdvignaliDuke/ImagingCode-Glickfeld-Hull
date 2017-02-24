function rc = behavConstsNAM


tHostname = lower(hostname);
[s,tUsername] = dos('ECHO %USERNAME%');

switch tHostname
    case {'nuke'}
        rc.pathStr = '\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data';
        rc.dataPat = 'data-i%03d-%s.mat';
        
        if tUsername(1:5) == 'linds'
            rc.name = 'linds';
            rootDir = '\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis';
            rc.eyeOutputDir = fullfile(rootDir,'Behavior\EyeTracking');
            rc.eyeInputDir = '\\CRASH.dhe.duke.edu\data\home\ashley\Analysis\';
            rc.caOutputDir = fullfile(rootDir,'2P');
            rc.ashleyAnalysis = '\\CRASH.dhe.duke.edu\data\home\ashley\Analysis';
        elseif tUsername(1:5) == 'ashle'
            rc.name = 'ashle';
            rootDir = '\\CRASH.dhe.duke.edu\data\home\ashley\Analysis';
            rc.eyeOutputDir = fullfile(rootDir, 'nAM Summaries','pupil');
            rc.eyeInputDir = '\\CRASH.dhe.duke.edu\data\home\ashley\Analysis\';
            rc.caOutputDir = fullfile(rootDir, 'nAM Summaries');
            rc.ashleyAnalysis = fullfile(rootDir);
            rc.ashleyData = '\\CRASH.dhe.duke.edu\data\home\ashley\data';
        end
end

rc.fitOutputTextCols = { 'DateStr', 'DataBlock', 'DateTimeStarted', 'PdfFigFilename', 'MatFilename' };
rc.indexTextCols = { 'DateStr', 'DataBlock', 'TrialRangeToUse', 'Notes' };

rc.fhWeibull = @(p,xs) (p(4) + (p(3)-p(4)) .* (1 - exp( -(xs./p(1)) .^ p(2))));


%%%%%%%%%%% simple simple functions

rc.computeFName = @subComputeFName;

function outFName = subComputeFName(subjNum, dateStr)
    outFName = fullfile(rc.pathStr, sprintf(rc.dataPat, subjNum, deblank(dateStr)));
end



end

