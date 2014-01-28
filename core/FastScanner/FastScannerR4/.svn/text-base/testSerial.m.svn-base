function testSerial(testType)
%$Id$
global fs;

if nargin < 1, testType = 'init'; end

switch testType
    case 'init'
        nReps = 1000;
        for iR = 1:nReps
            subReInitStage;
            disp(fs.position.z);
        end
        
    case 'random'
        nReps = 10000;
        subReInitStage;
        for iR = 1:nReps
            subDoRandTest
        end
    case 'stack'
        subReInitStage;
        subGetStack;
end

%%%%%%%%%%%%%%%%%%%%%%

function subGetStack
global fs
randStep = round(rand*10)+1;
randFrN = round(rand*10+2)*32;
randNSteps = round(rand*20);

subSetStrWCb(fs.handles.txtFramesPerStep, num2str(randFrN));

subSetStrWCb(fs.handles.txtStepZ, num2str(randStep));
subSetStrWCb(fs.handles.txtStartZ, num2str(fs.position.z));
subSetStrWCb(fs.handles.txtNofSteps, num2str(randNSteps));

% remove files
assert(strcmp(fs.DAQ.BaseFileName, 'E:\Data\tst'));
[succ,msg]=rmdir('E:\Data\tst', 's');  % recursively delete
%delete(fullfile(fs.DAQ.BaseFileName, '*.dad'));

% run acq
%tF = feval(str2func(fs.handles.btnStartCycles, 'ButtonDownFcn'));
FastScanner('btnStartCycles_Callback', fs.handles.btnStartCycles, ...
    [], guidata(fs.handles.btnStartCycles));

%%%%%%%%%%%%%%%%%%%%%%

function subReInitStage(doChangeParams)
global fs;

out1 = fs.stage.comport;
if ~isempty(out1)
    if isvalid(out1)
        set(fs.stage.comport, ...
            'Parity', 'none' , ...
            'Terminator', 'CR', ...
            'StopBits', 1, ...
            'DataBits',8, ...
            'FlowControl','hardware', ...
            'BaudRate',9600);  % was 19200
    
    end
end

out1=instrfind('Type', 'serial');
if ~isempty(out1)
    fclose(out1);
    delete(out1);
end
fs.stage.comport = [];

out1=instrfind;
if ~isempty(out1)
    keyboard;
end
IniStage;

%%%%%%%%%%%%%%%%%%%%%%

function subDoRandTest
% things we can do
%   re init stage
%   move to random pos (cannot use _nowait - it will crash stage and turn
%   off motor)
%
global fs;

nTests = 3;
rN = ceil(rand*(nTests));

switch rN
    case 1
        % init
        subReInitStage(true);
        fprintf(1, 'done: re-init.  Position: %d\n', fs.position.z);
    case 2
        % move some random interval
        currPos = fs.position.z;
        randInt = rand*5000-2500;  % uniform
        %randInt = randn*100; % gaussian sigma = 100um
        newPos = round(currPos + randInt);  
        SetPositionZ(newPos);
        fprintf(1, 'done: moved.  Old: %d; new %d; diff %d\n', ...
            currPos, newPos, newPos-currPos);
    case 3
        subGetStack
        fprintf(1, 'Done: acq stack\n');
    otherwise
        return
end

%%%%%%%5

function subSetStrWCb(h, newStr)
get(h, 'Type');
set(h, 'String', newStr);
gcbo = h;
cb = get(h, 'Callback');
eval(cb);
clear gcbo;


        
        
