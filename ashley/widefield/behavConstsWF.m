function rc = behavConstsWF(mouse)


tHostname = lower(hostname);
[s,tUsername] = dos('ECHO %USERNAME%');

switch tHostname      
    case {'nuke'}
        if tUsername(1:7) == 'lindsey'
            rc.pathStr = 'Z:\home\andrew\Behavior\Data';
            rc.dataPat = 'data-i%03d-%s.mat';
            rootDir = 'Z:\home\lindsey\Analysis\Widefield_imaging';
            rc.indexFilename = fullfile(rootDir, mouse, [mouse '-days.xls']);
            rc.outputFilename = fullfile(rootDir, mouse, [mouse '-fits.xls']);
            rc.structOutput = fullfile(rootDir, mouse);
            rc.fitOutputPdfDir = fullfile(rootDir, mouse, 'output\pdfFits');
            rc.fitOutputMatDir = fullfile(rootDir, mouse, 'output\fitMatStats');
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
