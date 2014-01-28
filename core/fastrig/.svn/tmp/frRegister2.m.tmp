function frRegister2(newdir,varargin)
%FRREGISTER
% FRREGISTER(NEWDIR)
% FRREGISTER(...,VARARGIN) where VARARGIN is list of property pairs 
%  Overwrite - Overwrites existing image directory [ true | {false} ]
%  Oversampling - Subpixel precision [{1},2,...]
%  DoRecurse - [ {true} | false ] 
%  Engine - [ {'subpixel'} | 'turboreg' ]

%%
if iscell(newdir)
    for item = 1:length(newdir)
        frRegister2(newdir{item},varargin{:});
    end;
    return;
end

% process all subdirectories
if ~sum(newdir==filesep)
    exptnames=frGetAllExpts(newdir)
    for index = 1:length(exptnames)
        frRegister2(exptnames{index},varargin{:});
    end
    return;
end;
    
%%
if ~verLessThan('matlab','7.5')
    prevn = maxNumCompThreads(2);
end;

defaultopts = {'Overwrite',false,'Oversampling',1,'DoRecurse',true,...
               'Engine','subpixel','nPlanes',1,'TargetFileNs',[]};

options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);

for iarg = 1:2:length(varargin)
    options.(varargin{iarg}) =varargin{iarg+1}
end

% dirs = frGetDirs;
% 
% imdir = fullfile(dirs.images,newdir);
% 
% if exist(imdir) ~= 7
%     str = sprintf('directory ''%s'' does not exist',imdir);
%     error(str);
% end

%% path names
[animal,name]=fileparts(newdir);
expt = frGetExpt(animal,name);
greenTargetPn = expt.dirs.reggreenpn;

[dummy,greenTargetFn] = fileparts(expt.filenames.registered);

%redTargetPn = fullfile(dirs.registered,expt.animal,expt.name,'red',expt.filenames.registeredred);

list = dir(fullfile(greenTargetPn,'*.tif'));
if length(list) && ~options.Overwrite
    fprintf('ignored. Already exists.\n');
	return;
end

if ~exist(expt.dirs.regrootpn)
    mkdir(expt.dirs.regrootpn);
end

if length(list) && options.Overwrite
    fprintf('deleted existing files');
    delete(fullfile(greenTargetPn,'*.tif'));    
    delete(fullfile(greenTargetPn,'*.mat'));
    delete(fullfile(expt.dirs.regrootpn,'*.log'));
end

fprintf('\n');

diary(fullfile(expt.dirs.regrootpn,'diary.log'));
if ~verLessThan('matlab','7.8')
    c = onCleanup(@()diary('off'));
end;
disp(options);

list = dir(fullfile(expt.dirs.rawgreenpn,'*.tif')); 
nFiles = length(list);

%% compute target, 128 (or less) frames from middle of the stack
fprintf('Computed target\n');

if isempty(options.TargetFileNs)
    options.TargetFileNs = 1:min(8,nFiles);
end

green = readtiff(expt.dirs.rawgreenpn,options.TargetFileNs);
[ny,nx,nframes] = size(green);

nP = options.nPlanes ;
os = options.Oversampling;

for iP = 1:nP
    targets(:,:,iP) =mean(green(:,:,iP:nP:end),3);
end;

% for iP = 1:nP
%     this = mean(green(:,:,iP:nP:end),3);
%     [out,reg] = stackRegister(green(:,:,iP:nP:end),this,os);
%     targets1(:,:,iP) = this;
%     targets(:,:,iP) = mean(reg,3);
% end;
% stackAnimator(cat(1,targets1,targets2))

%% register one file at a time
fprintf('Registering\n');
nn = 0;
t = cputime;

for iFile =1:nFiles            
    stack = readtiff(expt.dirs.rawgreenpn,iFile);
    siz = size(stack);
    
    iFr = 1:siz(3);
    
    % figure;imagesc(stack(:,:,1))
    % figure;imagesc(targets(:,:,1))
    
    fn = sprintf('%s%06d.tif',greenTargetFn,iFile);
    for iS = 1:nP % loop over starting indices
        iP = mod(nn + iS - 1,nP) + 1;

%         figure;
%         subplot(2,1,1);imagesq(mean(stack(:,:,iS:nP:end),3));
%         subplot(2,1,2);imagesq(targets(:,:,iP));
%         drawnow;
        
        if strcmp(options.Engine,'subpixel')
            [outputs{iFile},stack_reg] = stackRegister(stack(:,:,iS:nP:end),targets(:,:,iP), os);
            
        elseif strcmp(options.Engine,'turboreg')
            stack_reg = RegMulticore(stack(:,:,iS:nP:end),targets(:,:,iP), [], dirs.multicore, 3);
            outputs{iFile}=[];
        end;

        if nP > 1
            thisPn = fullfile(greenTargetPn,sprintf('Plane%i',iP),fn);
        else
            thisPn = fullfile(greenTargetPn,fn);
        end;
        writetiff(stack_reg,thisPn,'uint16');
%         save(fullfile(greenTargetPn,expt.filenames.shifts),'outputs');
    end; % iS    
    nn = nn + siz(3);
    fprintf('%s%06d.tif (%2.1f fps)\n',greenTargetFn,iFile,nn/(cputime-t));    
end; % iFile
  
if ~verLessThan('matlab','7.5')    
	maxNumCompThreads(prevn);  % restore previous value
end

diary off;

return;
