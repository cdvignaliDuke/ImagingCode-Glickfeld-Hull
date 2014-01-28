function frRegister(newdir,varargin)
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
        frRegister(newdir{item},varargin{:});
    end;
    return;
end

% process all subdirectories
if ~sum(newdir==filesep)
    exptnames=frGetAllExpts(newdir)
    for index = 1:length(exptnames)
        frRegister(exptnames{index},varargin{:});
    end
    return;
end
    
%%
if ~verLessThan('matlab','7.5')
    prevn = maxNumCompThreads(2);
end

defaultopts = {'Overwrite',false,'Oversampling',1,'DoRecurse',true,...
               'Engine','subpixel'};

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
nfiles = length(list);

%% compute target, 128 (or less) frames from middle of the stack
fprintf('Computed target\n');
mid = fix(nfiles/2);
if mid == 0 && nfiles == 1
    mid = 1;
elseif nfiles == 0
    error('No files found');
end
green = readtiff(expt.dirs.rawgreenpn,mid);
[ny,nx,nframes] = size(green);
mid = fix(nframes/2);
sel = intersect(1:nframes,mid:mid+127);
target = mean(green(:,:,sel),3);

%% register one file at a time
fprintf('Registering\n');

nn = 0;
t = cputime;

for iFile =1:nfiles
    green = readtiff(expt.dirs.rawgreenpn,iFile);
    siz = size(green);
    
    if strcmp(options.Engine,'subpixel')
        [outputs{iFile},green_reg] = stackRegister(green,target, options.Oversampling);
    elseif strcmp(options.Engine,'turboreg')
        green_reg = RegMulticore(green, target, [], dirs.multicore, 3);
        outputs{iFile}=[];
    end
    
    fn = sprintf('%s%06d.tif',greenTargetFn,iFile);
    writetiff(green_reg,fullfile(greenTargetPn,fn),'uint16');
    save(fullfile(greenTargetPn,expt.filenames.shifts),'outputs');
    nn = nn + siz(3);
    fprintf('%s%06d.tif (%2.1f fps)\n',greenTargetFn,iFile,nn/(cputime-t));
end

if ~verLessThan('matlab','7.5')    
	maxNumCompThreads(prevn);  % restore previous value
end

diary off;

return;
