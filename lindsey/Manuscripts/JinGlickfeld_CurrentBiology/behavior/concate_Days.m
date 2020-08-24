%% loop all the concatenate ones
% clear all
%
% clc

rootDir = 'Z:\miaomiao\behavior\Behavior';
rootDir2 = 'Z:\miaomiao\Analysis\Behavior\select and analysis';
rootDir3 = 'Z:\andrew\Behavior\Data';
%rootDir3 =  'Y:\Data\behavior\ISI-human-mat';
rc.indexFilename = fullfile(rootDir, 'experimentIndexes\subj-days-mj.xls');
rc.goodfitFilename = fullfile(rootDir2,'goodfit.xlsx');
rc.fitOutputMatDir = fullfile(rootDir, 'output\fitMatStats');
homeDir = 'Y:\Analysis\Behavior\select and analysis';
%xlsfile = fullfile(homeDir, 'LED-full.xlsx');
%xlsfile = fullfile(homeDir, 'speed.xlsx');
%xlsfile = fullfile(homeDir, 'contrast.xlsx');
%xlsfile = fullfile(homeDir, 'randomISI.xlsx');
xlsfile = fullfile(homeDir, 'randombase.xlsx');
%xlsfile = fullfile(homeDir, 'randomContrast.xlsx');
%xlsfile = fullfile(homeDir, 'ISI-LED.xlsx');
%xlsfile = fullfile(homeDir, 'V1 inhibition\summary.xlsx')
%%
%dataPat = 'data-i%03d-%s.mat';
% Need a tag for seperate sessions...
[Gdata, Gtext, Graw] = xlsread(xlsfile,7);
IDs = unique(Gdata(:,1));

% make new concatenate stuff allow for trial history analysis

for id_x = 1:length(IDs)
    
    inputcombined.trialOutcomeCell = [];
    inputcombined.reactTimesMs = [];
    inputcombined.holdTimesMs = [];
    inputcombined.tGratingDirectionDeg = [];
    inputcombined.tBaseGratingDirectionDeg = [];
    
    inputcombined.tGratingContrast = [];
    inputcombined.nCyclesOn = [];
    inputcombined.tCyclesOn = [];
    inputcombined.tCycleLaser = [];
    inputcombined.tseq = [];
    
    
    inputcombined.tStimOffTimes = [];
    inputcombined.stimOnTimeMs = [];
    inputcombined.tooFastTimeMs = [];
    inputcombined.targetStimOnMs = [];
    inputcombined.tLeverReleaseTimeMs = [];
    inputcombined.tLeverPressTimeMs = [];
    
    inputcombined.tBaseDotSpeedDPS = [];
    inputcombined.tDotSpeedDPS = [];
    inputcombined.tBaseDotCoherence = [];
    inputcombined.tDotCoherence = [];
    inputcombined.tDotDirectionDeg = [];
    inputcombined.tDotSizeDeg = [];
    inputcombined.tDotFieldSizeDeg = [];
    
    inputcombined.tBaseGratingContrast = [];
    
  
    
    inputcombined2.trialOutcomeCell = [];
    inputcombined2.reactTimesMs = [];
    inputcombined2.holdTimesMs = [];
    inputcombined2.tGratingDirectionDeg = [];
    inputcombined2.tBaseGratingDirectionDeg = [];
    inputcombined2.tGratingContrast = [];
    inputcombined2.nCyclesOn = [];
    inputcombined2.tCyclesOn = [];
    inputcombined2.tCycleLaser = [];
    inputcombined2.tseq = [];
    
    
    inputcombined2.tStimOffTimes = [];
    inputcombined2.stimOnTimeMs = [];
    inputcombined2.tooFastTimeMs = [];
    inputcombined2.targetStimOnMs = [];
    inputcombined2.tLeverReleaseTimeMs = [];
    inputcombined2.tLeverPressTimeMs = [];
    
    subjNum = IDs(id_x);
    matfiles = [];
    datestring = [];
    date = [];
    trial_range = [];
    trial = [];
    
    matfiles = Gdata(:,1)==subjNum ;
    date = Gtext(2:end, 2);
    trial = Gtext(2:end, 3);
    datestring = date(matfiles);
    trial_range = trial(matfiles);
    %     sessions = zeros(1,length(matfiles));
    
    for i = 1: sum(matfiles)
        idate = num2str(cell2mat(datestring(i,1)));
        i_range =cell2mat( trial_range(i,1));
        dataPat = 'data-i%03d-%s.mat';
        outname = [];
        filename = dir(fullfile(rootDir3, ['data-i'  num2str(subjNum) '-' idate '*']));
        outname = fullfile(rootDir3,filename.name)
        load(outname);
        % quality control of early releases, get rid of very early fidget
        % fidget is defined as trials have hold time less than 0.4s.
        hold_tmp = double(cell2mat(input.holdTimesMs(eval(i_range))))./1000; % holdtime_S of the animal
        %fidget = sum(hold_tmp<0.4); % number of fidget trials, which should be excluded for analysis
        % redefine the fidget as < 0.9 s on 2017.10.03
        fidget = sum(hold_tmp<0.9);
        earlys = (sum(strcmp(input.trialOutcomeCell(eval(i_range)),'failure'))-fidget)./(numel(eval(i_range))-fidget)
        successes = sum(strcmp(input.trialOutcomeCell(eval(i_range)),'success'))./(numel(eval(i_range))-fidget)
        
        % if it is speed detection
        if isfield(input,'doSpeedDetect')
        if earlys <= 0.35 && successes>=0.5 && input.doSpeedDetect==1 && input.doDirDetect ==0
            %             sessions(i) = 1;
            %if earlys > 0.35  && length(unique(cell2mat(input.tBlock2TrialNumber)))==1
            inputcombined.trialOutcomeCell = [inputcombined.trialOutcomeCell  input.trialOutcomeCell(eval(i_range))];
            inputcombined.reactTimesMs = [inputcombined.reactTimesMs input.reactTimesMs(eval(i_range))];
            inputcombined.holdTimesMs = [inputcombined.holdTimesMs input.holdTimesMs(eval(i_range))];
            inputcombined.tBaseDotSpeedDPS = [inputcombined.tBaseDotSpeedDPS input.tBaseDotSpeedDPS(eval(i_range))];
            inputcombined.tDotSpeedDPS = [inputcombined.tDotSpeedDPS input.tDotSpeedDPS(eval(i_range))];
            inputcombined.tBaseDotCoherence = [inputcombined.tBaseDotCoherence input.tBaseDotCoherence(eval(i_range))];
            inputcombined.tDotCoherence = [inputcombined.tDotCoherence input.tDotCoherence(eval(i_range))];
            
            inputcombined.tDotDirectionDeg = [inputcombined.tDotDirectionDeg input.tDotDirectionDeg(eval(i_range))];
            inputcombined.tDotSizeDeg = [inputcombined.tDotSizeDeg input.tDotSizeDeg(eval(i_range))];
            inputcombined.tDotFieldSizeDeg = [inputcombined.tDotFieldSizeDeg input.tDotFieldSizeDeg(i_range)];
            
            
            
            unique(double(cell2mat(inputcombined.tDotSpeedDPS(eval(i_range)))))
            
            range = [];
            range = eval(i_range);
            range = range - range(1)+1;
            inputcombined.tseq = [inputcombined.tseq,range];
            
            
            % fix the reaction time differences
            inputcombined.RTwindowMs =  input.reactTimeMs;
            inputcombined.ID = subjNum;
            input.doOriDetect = 0;
            input.doContrastDetect = 0;
            
        end
        end
        
        if isfield(input,'doRandCon')
            if earlys <= 0.35 && successes>=0.5 && length(unique(cell2mat(input.tBlock2TrialNumber)))==1 &&  input.doLaserStim~=1 && input.doOriDetect ==1 && input.doRandCon ==1
                %             sessions(i) = 1;
                %if earlys > 0.35  && length(unique(cell2mat(input.tBlock2TrialNumber)))==1
                inputcombined.trialOutcomeCell = [inputcombined.trialOutcomeCell  input.trialOutcomeCell(eval(i_range))];
                inputcombined.reactTimesMs = [inputcombined.reactTimesMs input.reactTimesMs(eval(i_range))];
                inputcombined.holdTimesMs = [inputcombined.holdTimesMs input.holdTimesMs(eval(i_range))];
                inputcombined.tGratingDirectionDeg = [inputcombined.tGratingDirectionDeg input.tGratingDirectionDeg(eval(i_range))];
                inputcombined.tBaseGratingDirectionDeg = [inputcombined.tBaseGratingDirectionDeg input.tBaseGratingDirectionDeg(eval(i_range))];
                inputcombined.tBaseGratingContrast = [inputcombined.tBaseGratingContrast input.tBaseGratingContrast(eval(i_range))];
                unique(double(cell2mat(inputcombined.tBaseGratingContrast)))
                unique(double(cell2mat(input.tGratingDirectionDeg(eval(i_range)))))
                
                inputcombined.nCyclesOn = [inputcombined.nCyclesOn input.nCyclesOn(eval(i_range))];
                inputcombined.tCyclesOn = [inputcombined.tCyclesOn input.tCyclesOn(eval(i_range))];
                range = [];
                range = eval(i_range);
                range = range - range(1)+1;
                inputcombined.tseq = [inputcombined.tseq,range];
                if isfield(input,'tStimOffTimes')
                    inputcombined.tStimOffTimes = [inputcombined.tStimOffTimes input.tStimOffTimes(eval(i_range))];
                end
                
                inputcombined.targetStimOnMs =[inputcombined.targetStimOnMs input.targetStimOnMs(eval(i_range))];
                inputcombined.tLeverReleaseTimeMs = [inputcombined.tLeverReleaseTimeMs input.tLeverReleaseTimeMs(eval(i_range))];
                inputcombined.tLeverPressTimeMs = [inputcombined.tLeverPressTimeMs input.tLeverPressTimeMs(eval(i_range))];
                inputcombined.stimOnTimeMs =  input.stimOnTimeMs;
                inputcombined.stimOffTimeMs =  input.stimOffTimeMs;
                inputcombined.tooFastTimeMs =  input.tooFastTimeMs;
                % fix the reaction time differences
                inputcombined.RTwindowMs =  input.reactTimeMs;
                inputcombined.ID = subjNum;
                
            end
            
            if earlys <= 0.35 && successes>=0.5 && length(unique(cell2mat(input.tBlock2TrialNumber)))==1 &&  input.doLaserStim~=1 && input.doOriDetect ==1 && input.doRandCon ==0
                %             sessions(i) = 1;
                %if earlys > 0.35  && length(unique(cell2mat(input.tBlock2TrialNumber)))==1
                inputcombined.trialOutcomeCell = [inputcombined.trialOutcomeCell  input.trialOutcomeCell(eval(i_range))];
                inputcombined.reactTimesMs = [inputcombined.reactTimesMs input.reactTimesMs(eval(i_range))];
                inputcombined.holdTimesMs = [inputcombined.holdTimesMs input.holdTimesMs(eval(i_range))];
                inputcombined.tGratingDirectionDeg = [inputcombined.tGratingDirectionDeg input.tGratingDirectionDeg(eval(i_range))];
                inputcombined.tBaseGratingDirectionDeg = [inputcombined.tBaseGratingDirectionDeg input.tBaseGratingDirectionDeg(eval(i_range))];
                
                unique(double(cell2mat(input.tGratingDirectionDeg(eval(i_range)))))
                
                inputcombined.nCyclesOn = [inputcombined.nCyclesOn input.nCyclesOn(eval(i_range))];
                inputcombined.tCyclesOn = [inputcombined.tCyclesOn input.tCyclesOn(eval(i_range))];
                range = [];
                range = eval(i_range);
                range = range - range(1)+1;
                inputcombined.tseq = [inputcombined.tseq,range];
                if isfield(input,'tStimOffTimes')
                    inputcombined.tStimOffTimes = [inputcombined.tStimOffTimes input.tStimOffTimes(eval(i_range))];
                end
                
                inputcombined.targetStimOnMs =[inputcombined.targetStimOnMs input.targetStimOnMs(eval(i_range))];
                inputcombined.tLeverReleaseTimeMs = [inputcombined.tLeverReleaseTimeMs input.tLeverReleaseTimeMs(eval(i_range))];
                inputcombined.tLeverPressTimeMs = [inputcombined.tLeverPressTimeMs input.tLeverPressTimeMs(eval(i_range))];
                inputcombined.stimOnTimeMs =  input.stimOnTimeMs;
                inputcombined.stimOffTimeMs =  input.stimOffTimeMs;
                inputcombined.tooFastTimeMs =  input.tooFastTimeMs;
                % fix the reaction time differences
                inputcombined.RTwindowMs =  input.reactTimeMs;
                inputcombined.ID = subjNum;
                
            end
            
        else
            
            
            if earlys <= 0.35 && successes>=0.5 && length(unique(cell2mat(input.tBlock2TrialNumber)))==1 &&  input.doLaserStim~=1 && input.doOriDetect ==1
                %             sessions(i) = 1;
                %if earlys > 0.35  && length(unique(cell2mat(input.tBlock2TrialNumber)))==1
                inputcombined.trialOutcomeCell = [inputcombined.trialOutcomeCell  input.trialOutcomeCell(eval(i_range))];
                inputcombined.reactTimesMs = [inputcombined.reactTimesMs input.reactTimesMs(eval(i_range))];
                inputcombined.holdTimesMs = [inputcombined.holdTimesMs input.holdTimesMs(eval(i_range))];
                inputcombined.tGratingDirectionDeg = [inputcombined.tGratingDirectionDeg input.tGratingDirectionDeg(eval(i_range))];
                inputcombined.tBaseGratingDirectionDeg = [inputcombined.tBaseGratingDirectionDeg input.tBaseGratingDirectionDeg(eval(i_range))];
                
                unique(double(cell2mat(input.tGratingDirectionDeg(eval(i_range)))))
                
                inputcombined.nCyclesOn = [inputcombined.nCyclesOn input.nCyclesOn(eval(i_range))];
                inputcombined.tCyclesOn = [inputcombined.tCyclesOn input.tCyclesOn(eval(i_range))];
                range = [];
                range = eval(i_range);
                range = range - range(1)+1;
                inputcombined.tseq = [inputcombined.tseq,range];
                if isfield(input,'tStimOffTimes')
                    inputcombined.tStimOffTimes = [inputcombined.tStimOffTimes input.tStimOffTimes(eval(i_range))];
                end
                
                inputcombined.targetStimOnMs =[inputcombined.targetStimOnMs input.targetStimOnMs(eval(i_range))];
                inputcombined.tLeverReleaseTimeMs = [inputcombined.tLeverReleaseTimeMs input.tLeverReleaseTimeMs(eval(i_range))];
                inputcombined.tLeverPressTimeMs = [inputcombined.tLeverPressTimeMs input.tLeverPressTimeMs(eval(i_range))];
                inputcombined.stimOnTimeMs =  input.stimOnTimeMs;
                inputcombined.stimOffTimeMs =  input.stimOffTimeMs;
                inputcombined.tooFastTimeMs =  input.tooFastTimeMs;
                % fix the reaction time differences
                inputcombined.RTwindowMs =  input.reactTimeMs;
                inputcombined.ID = subjNum;
                
            end
            
        end
        
        % seperate out the condition is single pulse inhibition or boost
        if earlys <= 0.5 && successes>=0.4 && input.doLaserStim==1 && input.doOriDetect ==1
            %if earlys <= 0.55 && successes>=0.35 && input.doLaserStim==1 && input.doOriDetect ==1
            % store the trials without pulse inhibition in inputcombined
            % if it is early trials, condition should <0; if not early <=0
            temp = input.trialOutcomeCell(eval(i_range));
            Early = strcmp(temp,'failure');
            block2 = [];
            block_Nearly =[];
            block_early = [];
            block_Nearly = (cell2mat(input.tCycleLaser(eval(i_range))) - (cell2mat(input.tCyclesOn(eval(i_range)))+1))<=0;
            block_early = (cell2mat(input.tCycleLaser(eval(i_range))) - (cell2mat(input.tCyclesOn(eval(i_range)))+1))<0;
            block2(Early) = block_early(Early);
            block2(~Early) = block_Nearly(~Early );
            
            % for non led condition store in inputcombined
            
            inputcombined.trialOutcomeCell = [inputcombined.trialOutcomeCell   temp(block2==0)];
            temp = input.reactTimesMs(eval(i_range));
            inputcombined.reactTimesMs = [inputcombined.reactTimesMs temp(block2==0)];
            temp = input.holdTimesMs(eval(i_range));
            inputcombined.holdTimesMs = [inputcombined.holdTimesMs temp(block2==0)];
            temp = input.tGratingDirectionDeg(eval(i_range));
            inputcombined.tGratingDirectionDeg = [inputcombined.tGratingDirectionDeg temp(block2==0)];
            unique(double(cell2mat(inputcombined.tGratingDirectionDeg)))
            temp = input.nCyclesOn(eval(i_range));
            inputcombined.nCyclesOn = [inputcombined.nCyclesOn temp(block2==0)];
            temp = input.tCyclesOn(eval(i_range));
            inputcombined.tCyclesOn = [inputcombined.tCyclesOn temp(block2==0)];
            % store when the laser cycle is on. no mening here
            temp = input.tCycleLaser(eval(i_range));
            inputcombined.tCycleLaser = [inputcombined.tCycleLaser temp(block2==0)];
            
            
            
            if isfield(input,'tStimOffTimes')
                temp = input.tStimOffTimes(eval(i_range));
                inputcombined.tStimOffTimes = [inputcombined.tStimOffTimes temp(block2==0)];
            end
            temp = input.targetStimOnMs(eval(i_range));
            inputcombined.targetStimOnMs =[inputcombined.targetStimOnMs temp(block2==0)];
            temp = input.tLeverReleaseTimeMs(eval(i_range));
            inputcombined.tLeverReleaseTimeMs = [inputcombined.tLeverReleaseTimeMs temp(block2==0)];
            temp = input.tLeverPressTimeMs(eval(i_range));
            inputcombined.tLeverPressTimeMs = [inputcombined.tLeverPressTimeMs temp(block2==0)];
            inputcombined.stimOnTimeMs =  input.stimOnTimeMs;
            inputcombined.stimOffTimeMs =  input.stimOffTimeMs;
            inputcombined.tooFastTimeMs =  input.tooFastTimeMs;
            % this is the react time window
            inputcombined.RTwindowMs =  input.reactTimeMs;
            inputcombined.ID = subjNum;
            inputcombined.doLaserstim =1;
            
            % for led condition store in inputcombined2
            
            temp = input.trialOutcomeCell(eval(i_range));
            inputcombined2.trialOutcomeCell = [inputcombined2.trialOutcomeCell   temp(block2==1)];
            temp = input.reactTimesMs(eval(i_range));
            inputcombined2.reactTimesMs = [inputcombined2.reactTimesMs temp(block2==1)];
            temp = input.holdTimesMs(eval(i_range));
            inputcombined2.holdTimesMs = [inputcombined2.holdTimesMs temp(block2==1)];
            temp = input.tGratingDirectionDeg(eval(i_range));
            inputcombined2.tGratingDirectionDeg = [inputcombined2.tGratingDirectionDeg temp(block2==1)];
            unique(double(cell2mat(inputcombined2.tGratingDirectionDeg)))
            temp = input.nCyclesOn(eval(i_range));
            inputcombined2.nCyclesOn = [inputcombined2.nCyclesOn temp(block2==1)];
            temp = input.tCyclesOn(eval(i_range));
            inputcombined2.tCyclesOn = [inputcombined2.tCyclesOn temp(block2==1)];
            
            % store when the laser cycle is on. no mening here
            temp = input.tCycleLaser(eval(i_range));
            inputcombined2.tCycleLaser = [inputcombined2.tCycleLaser temp(block2==1)];
            
            
            if isfield(input,'tStimOffTimes')
                temp = input.tStimOffTimes(eval(i_range));
                inputcombined2.tStimOffTimes = [inputcombined2.tStimOffTimes temp(block2==1)];
            end
            temp = input.targetStimOnMs(eval(i_range));
            inputcombined2.targetStimOnMs =[inputcombined2.targetStimOnMs temp(block2==1)];
            temp = input.tLeverReleaseTimeMs(eval(i_range));
            inputcombined2.tLeverReleaseTimeMs = [inputcombined2.tLeverReleaseTimeMs temp(block2==1)];
            temp = input.tLeverPressTimeMs(eval(i_range));
            inputcombined2.tLeverPressTimeMs = [inputcombined2.tLeverPressTimeMs temp(block2==1)];
            inputcombined2.stimOnTimeMs =  input.stimOnTimeMs;
            inputcombined2.stimOffTimeMs =  input.stimOffTimeMs;
            inputcombined2.tooFastTimeMs =  input.tooFastTimeMs;
            % this is the react time window
            inputcombined2.RTwindowMs =  input.reactTimeMs;
            inputcombined2.ID = subjNum;
            inputcombined2.doLaserstim =1;
            
            
        end
        
        
        
        
        
        
        % losen the criterion if it is led condition
        if earlys <= 0.5 && successes>=0.4 && length(unique(cell2mat(input.tBlock2TrialNumber)))==2 &&  input.doLaserStim~=1 && input.doOriDetect ==1
            block2 =[];
            block2= cell2mat(input.tBlock2TrialNumber(eval(i_range)));
            % for non led condition store in inputcombined
            temp = input.trialOutcomeCell(eval(i_range));
            inputcombined.trialOutcomeCell = [inputcombined.trialOutcomeCell   temp(block2==0)];
            temp = input.reactTimesMs(eval(i_range));
            inputcombined.reactTimesMs = [inputcombined.reactTimesMs temp(block2==0)];
            temp = input.holdTimesMs(eval(i_range));
            inputcombined.holdTimesMs = [inputcombined.holdTimesMs temp(block2==0)];
            temp = input.tGratingDirectionDeg(eval(i_range));
            inputcombined.tGratingDirectionDeg = [inputcombined.tGratingDirectionDeg temp(block2==0)];
            unique(double(cell2mat(inputcombined.tGratingDirectionDeg)))
            temp = input.nCyclesOn(eval(i_range));
            inputcombined.nCyclesOn = [inputcombined.nCyclesOn temp(block2==0)];
            temp = input.tCyclesOn(eval(i_range));
            inputcombined.tCyclesOn = [inputcombined.tCyclesOn temp(block2==0)];
            
            range = [];
            range = eval(i_range);
            range = range - range(1)+1;
            inputcombined.tseq = [inputcombined.tseq,range(block2==0)];
            
            
            if isfield(input,'tStimOffTimes')
                temp = input.tStimOffTimes(eval(i_range));
                inputcombined.tStimOffTimes = [inputcombined.tStimOffTimes temp(block2==0)];
            end
            temp = input.targetStimOnMs(eval(i_range));
            inputcombined.targetStimOnMs =[inputcombined.targetStimOnMs temp(block2==0)];
            temp = input.tLeverReleaseTimeMs(eval(i_range));
            inputcombined.tLeverReleaseTimeMs = [inputcombined.tLeverReleaseTimeMs temp(block2==0)];
            temp = input.tLeverPressTimeMs(eval(i_range));
            inputcombined.tLeverPressTimeMs = [inputcombined.tLeverPressTimeMs temp(block2==0)];
            inputcombined.stimOnTimeMs =  input.stimOnTimeMs;
            inputcombined.stimOffTimeMs =  input.stimOffTimeMs;
            inputcombined.tooFastTimeMs =  input.tooFastTimeMs;
            % this is the react time window
            inputcombined.RTwindowMs =  input.reactTimeMs;
            inputcombined.ID = subjNum;
            
            % for led condition store in inputcombined2
            temp = input.trialOutcomeCell(eval(i_range));
            inputcombined2.trialOutcomeCell = [inputcombined2.trialOutcomeCell   temp(block2==1)];
            temp = input.reactTimesMs(eval(i_range));
            inputcombined2.reactTimesMs = [inputcombined2.reactTimesMs temp(block2==1)];
            temp = input.holdTimesMs(eval(i_range));
            inputcombined2.holdTimesMs = [inputcombined2.holdTimesMs temp(block2==1)];
            temp = input.tGratingDirectionDeg(eval(i_range));
            inputcombined2.tGratingDirectionDeg = [inputcombined2.tGratingDirectionDeg temp(block2==1)];
            unique(double(cell2mat(inputcombined2.tGratingDirectionDeg)))
            temp = input.nCyclesOn(eval(i_range));
            inputcombined2.nCyclesOn = [inputcombined2.nCyclesOn temp(block2==1)];
            temp = input.tCyclesOn(eval(i_range));
            inputcombined2.tCyclesOn = [inputcombined2.tCyclesOn temp(block2==1)];
            inputcombined2.tseq = [inputcombined2.tseq,range(block2==1)];
            
            if isfield(input,'tStimOffTimes')
                temp = input.tStimOffTimes(eval(i_range));
                inputcombined2.tStimOffTimes = [inputcombined2.tStimOffTimes temp(block2==1)];
            end
            temp = input.targetStimOnMs(eval(i_range));
            inputcombined2.targetStimOnMs =[inputcombined2.targetStimOnMs temp(block2==1)];
            temp = input.tLeverReleaseTimeMs(eval(i_range));
            inputcombined2.tLeverReleaseTimeMs = [inputcombined2.tLeverReleaseTimeMs temp(block2==1)];
            temp = input.tLeverPressTimeMs(eval(i_range));
            inputcombined2.tLeverPressTimeMs = [inputcombined2.tLeverPressTimeMs temp(block2==1)];
            inputcombined2.stimOnTimeMs =  input.stimOnTimeMs;
            inputcombined2.stimOffTimeMs =  input.stimOffTimeMs;
            inputcombined2.tooFastTimeMs =  input.tooFastTimeMs;
            % this is the react time window
            inputcombined2.RTwindowMs =  input.reactTimeMs;
            inputcombined2.ID = subjNum;
            
            
        end
        
        % if it is contrast detection
        
        if earlys <= 0.5 && successes>=0.4 && length(unique(cell2mat(input.tBlock2TrialNumber)))==2 &&  input.doLaserStim~=1 && input.doContrastDetect ==1
            
            if isfield(input,'nCyclesOn') % if use FS code
                
                block2 =[];
                block2= cell2mat(input.tBlock2TrialNumber(eval(i_range)));
                % for non led condition store in inputcombined
                temp = input.trialOutcomeCell(eval(i_range));
                inputcombined.trialOutcomeCell = [inputcombined.trialOutcomeCell   temp(block2==0)];
                temp = input.reactTimesMs(eval(i_range));
                inputcombined.reactTimesMs = [inputcombined.reactTimesMs temp(block2==0)];
                temp = input.holdTimesMs(eval(i_range));
                inputcombined.holdTimesMs = [inputcombined.holdTimesMs temp(block2==0)];
                temp = input.tGratingDirectionDeg(eval(i_range));
                inputcombined.tGratingDirectionDeg = [inputcombined.tGratingDirectionDeg temp(block2==0)];
                
                
                temp = input.tGratingContrast(eval(i_range));
                inputcombined.tGratingContrast = [inputcombined.tGratingContrast temp(block2==0)];
                unique(double(cell2mat(inputcombined.tGratingContrast)))
                
                range = [];
                range = eval(i_range);
                range = range - range(1)+1;
                inputcombined.tseq = [inputcombined.tseq,range(block2==0)];
                
                
                if isfield(input,'tStimOffTimes')
                    temp = input.tStimOffTimes(eval(i_range));
                    inputcombined.tStimOffTimes = [inputcombined.tStimOffTimes temp(block2==0)];
                end
                temp = input.targetStimOnMs(eval(i_range));
                inputcombined.targetStimOnMs =[inputcombined.targetStimOnMs temp(block2==0)];
                temp = input.tLeverReleaseTimeMs(eval(i_range));
                inputcombined.tLeverReleaseTimeMs = [inputcombined.tLeverReleaseTimeMs temp(block2==0)];
                temp = input.tLeverPressTimeMs(eval(i_range));
                inputcombined.tLeverPressTimeMs = [inputcombined.tLeverPressTimeMs temp(block2==0)];
                
                
                temp = input.nCyclesOn(eval(i_range));
                inputcombined.nCyclesOn = [inputcombined.nCyclesOn temp(block2==0)];
                temp = input.tCyclesOn(eval(i_range));
                inputcombined.tCyclesOn = [inputcombined.tCyclesOn temp(block2==0)];
                
                
                
                
                inputcombined.stimOnTimeMs =  input.stimOnTimeMs;
                inputcombined.stimOffTimeMs =  input.stimOffTimeMs;
                inputcombined.tooFastTimeMs =  input.tooFastTimeMs;
                % this is the react time window
                inputcombined.RTwindowMs =  input.reactTimeMs;
                inputcombined.ID = subjNum;
                
                % for led condition store in inputcombined2
                temp = input.trialOutcomeCell(eval(i_range));
                inputcombined2.trialOutcomeCell = [inputcombined2.trialOutcomeCell   temp(block2==1)];
                temp = input.reactTimesMs(eval(i_range));
                inputcombined2.reactTimesMs = [inputcombined2.reactTimesMs temp(block2==1)];
                temp = input.holdTimesMs(eval(i_range));
                inputcombined2.holdTimesMs = [inputcombined2.holdTimesMs temp(block2==1)];
                temp = input.tGratingDirectionDeg(eval(i_range));
                inputcombined2.tGratingDirectionDeg = [inputcombined2.tGratingDirectionDeg temp(block2==1)];
                unique(double(cell2mat(inputcombined2.tGratingDirectionDeg)))
                
                temp = input.tGratingContrast(eval(i_range));
                inputcombined2.tGratingContrast = [inputcombined2.tGratingContrast temp(block2==1)];
                unique(double(cell2mat(inputcombined2.tGratingContrast)))
                
                temp = input.nCyclesOn(eval(i_range));
                inputcombined2.nCyclesOn = [inputcombined2.nCyclesOn temp(block2==1)];
                temp = input.tCyclesOn(eval(i_range));
                inputcombined2.tCyclesOn = [inputcombined2.tCyclesOn temp(block2==1)];
                
                inputcombined2.tseq = [inputcombined2.tseq,range(block2==1)];
                
                if isfield(input,'tStimOffTimes')
                    temp = input.tStimOffTimes(eval(i_range));
                    inputcombined2.tStimOffTimes = [inputcombined2.tStimOffTimes temp(block2==1)];
                end
                temp = input.targetStimOnMs(eval(i_range));
                inputcombined2.targetStimOnMs =[inputcombined2.targetStimOnMs temp(block2==1)];
                temp = input.tLeverReleaseTimeMs(eval(i_range));
                inputcombined2.tLeverReleaseTimeMs = [inputcombined2.tLeverReleaseTimeMs temp(block2==1)];
                temp = input.tLeverPressTimeMs(eval(i_range));
                inputcombined2.tLeverPressTimeMs = [inputcombined2.tLeverPressTimeMs temp(block2==1)];
                inputcombined2.stimOnTimeMs =  input.stimOnTimeMs;
                inputcombined2.stimOffTimeMs =  input.stimOffTimeMs;
                inputcombined2.tooFastTimeMs =  input.tooFastTimeMs;
                % this is the react time window
                inputcombined2.RTwindowMs =  input.reactTimeMs;
                inputcombined2.ID = subjNum;
                
                
                
            else  % if use HDC code
                block2 =[];
                block2= cell2mat(input.tBlock2TrialNumber(eval(i_range)));
                % for non led condition store in inputcombined
                temp = input.trialOutcomeCell(eval(i_range));
                inputcombined.trialOutcomeCell = [inputcombined.trialOutcomeCell   temp(block2==0)];
                temp = input.reactTimesMs(eval(i_range));
                inputcombined.reactTimesMs = [inputcombined.reactTimesMs temp(block2==0)];
                temp = input.holdTimesMs(eval(i_range));
                inputcombined.holdTimesMs = [inputcombined.holdTimesMs temp(block2==0)];
                
                temp = input.tGratingDirectionDeg(eval(i_range));
                inputcombined.tGratingDirectionDeg = [inputcombined.tGratingDirectionDeg temp(block2==0)];
                
                
                temp = input.tGratingContrast(eval(i_range));
                inputcombined.tGratingContrast = [inputcombined.tGratingContrast temp(block2==0)];
                unique(double(cell2mat(inputcombined.tGratingContrast)))
                
                range = [];
                range = eval(i_range);
                range = range - range(1)+1;
                inputcombined.tseq = [inputcombined.tseq,range(block2==0)];
                
                
                
                
                inputcombined.tooFastTimeMs =  input.tooFastTimeMs;
                % this is the react time window
                inputcombined.RTwindowMs =  input.reactTimeMs;
                inputcombined.ID = subjNum;
                
                % for led condition store in inputcombined2
                temp = input.trialOutcomeCell(eval(i_range));
                inputcombined2.trialOutcomeCell = [inputcombined2.trialOutcomeCell   temp(block2==1)];
                temp = input.reactTimesMs(eval(i_range));
                inputcombined2.reactTimesMs = [inputcombined2.reactTimesMs temp(block2==1)];
                temp = input.holdTimesMs(eval(i_range));
                inputcombined2.holdTimesMs = [inputcombined2.holdTimesMs temp(block2==1)];
                temp = input.tGratingDirectionDeg(eval(i_range));
                inputcombined2.tGratingDirectionDeg = [inputcombined2.tGratingDirectionDeg temp(block2==1)];
                unique(double(cell2mat(inputcombined2.tGratingDirectionDeg)))
                
                temp = input.tGratingContrast(eval(i_range));
                inputcombined2.tGratingContrast = [inputcombined2.tGratingContrast temp(block2==1)];
                unique(double(cell2mat(inputcombined2.tGratingContrast)))
                
                inputcombined2.tseq = [inputcombined2.tseq,range(block2==1)];
                
                inputcombined2.tooFastTimeMs =  input.tooFastTimeMs;
                % this is the react time window
                inputcombined2.RTwindowMs =  input.reactTimeMs;
                inputcombined2.ID = subjNum;
                
                
            end
            
        end
        
    end
    
    
    
    
    if length(unique(cell2mat(input.tBlock2TrialNumber)))==1&& input.doLaserStim==0
        save([num2str(subjNum) 'concate-randomISI'], 'inputcombined')
        %save([num2str(subjNum) 'concate-randomISI-bad'], 'inputcombined')
    elseif length(unique(cell2mat(input.tBlock2TrialNumber)))==2 && input.doOriDetect ==1
        save([num2str(subjNum) 'concate-randomISI-LED'], 'inputcombined','inputcombined2')
    elseif length(unique(cell2mat(input.tBlock2TrialNumber)))==2 && input.doContrastDetect ==1
        save([num2str(subjNum) 'concate-contrast-LMLED'], 'inputcombined','inputcombined2')
    elseif length(unique(cell2mat(input.tBlock2TrialNumber)))==1&& input.doLaserStim==1
        save([num2str(subjNum) 'concate-randomISI' '-' num2str(unique(cell2mat(input.laserPowerMw))) 'mW' '.mat'], 'inputcombined','inputcombined2')
    elseif isfield(input,'doSpeedDetect') ==1
        save([num2str(subjNum) 'concate-speed'], 'inputcombined')
    end
    
    
    
end





