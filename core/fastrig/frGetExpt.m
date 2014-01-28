function expt = frGetExpt(animal,exptname)
%FRGETEXPT Creates experiment structure
% EXPT = FRGETEXPT(NEWDIR) where NEWDIR is string with 
% format $ANIMAL\$EXPTNAME
% EXPT = FRGETEXPT(ANIMAL,EXPTNAME)

% 10/01 id uses windows directory format

if nargin < 2
    [animal,exptname]=fileparts(animal);
end

expt.animal = animal;
expt.name = exptname;

expt.id = sprintf('%s\\%s',expt.animal,expt.name);

if nargin > 2
    expt.prot = protname;
else
    expt.prot = '';
end

expt.dirs = frGetDirs;

%% image root
expt.dirs.rawrootpn = fullfile(expt.dirs.images,expt.animal,expt.name);

if exist(expt.dirs.rawrootpn)~=7 
    expt.dirs.rawrootpn = fullfile(expt.dirs.backup,'images',expt.animal,expt.name);
end

%% raw data directories
scanlog = frParseScanLog(expt.dirs.rawrootpn); 

if isempty(scanlog)   
    pn = fullfile(expt.dirs.data,expt.animal,expt.name);
    if exist(pn)==7
        expt.dirs.dadsrootpn = pn;
        % read scan log if exists
        scanlog = frParseScanLog(expt.dirs.dadsrootpn);    
    end
end

expt.scanlog = scanlog;

if ~isempty(expt.scanlog)
    expt.desc = [expt.scanlog.start_date ' ' expt.id];
else
    expt.desc = expt.id;
end

%% green / red directories
list = dir(expt.dirs.rawrootpn);
if ~sum(strmatchl('_green',{list.name}))
    fprintf(1,'the directory \n')
    disp(expt.dirs.rawrootpn);
    fprintf(1,'needs a ''*_green'' and ''*_red'' subdirectory\n');
    error('with dir structure');
end;

if isempty(list)    
    expt.dirs.rawgreenpn = fullfile(expt.dirs.rawrootpn,'_green');
    expt.dirs.rawredpn = fullfile(expt.dirs.rawrootpn,'_red');
else    
    list(1:2)=[];
    expt.dirs.rawgreenpn = fullfile(expt.dirs.rawrootpn,list( strmatchl('_green',{list.name}) ).name);
    expt.dirs.rawredpn = fullfile(expt.dirs.rawrootpn,list( strmatchl('_red',{list.name}) ).name);
end

expt.dirs.regrootpn = fullfile(expt.dirs.registered,expt.animal,expt.name);

list = dir(fullfile(expt.dirs.regrootpn,'*'));

if ~isempty(list);list(1:2)=[];end;

if ~any([list.isdir])

    expt.dirs.reggreenpn = fullfile(expt.dirs.regrootpn,[expt.name '_green']);
    expt.dirs.regredpn = fullfile(expt.dirs.regrootpn,[expt.name '_red']);
    
else
    ind = find(strmatchl('green',{list.name}));
    if length(ind)
        expt.dirs.reggreenpn = fullfile(expt.dirs.regrootpn,list( ind ).name);
    else 
        expt.dirs.reggreenpn=[];
    end
    
    ind = find(strmatchl('red',{list.name}));
    if length(ind)
        expt.dirs.regredpn = fullfile(expt.dirs.regrootpn,list( ind ).name);
    else
        expt.dirs.regredpn =[];
    end
end

expt.dirs.pcarootpn = fullfile(expt.dirs.pca,expt.animal,expt.name);
list = dir(expt.dirs.pcarootpn);

if isempty(list)
    expt.dirs.pcagreenpn = fullfile(expt.dirs.pcarootpn,[expt.name '_green']);
else
    list(1:2)=[];
    expt.dirs.pcagreenpn = fullfile(expt.dirs.pcarootpn,list( strmatchl('_green',{list.name}) ).name);
end

%% analysis directories
expt.dirs.analrootpn = fullfile(expt.dirs.analysis,expt.animal,expt.name);

if exist(expt.dirs.analrootpn)~=7
    [success]=mkdir(expt.dirs.analrootpn);
    if ~success
        warning(sprintf('Could not create directory: %s',expt.dirs.analrootpn))
    end
end

% expt.dirs.figsOut = fullfile(expt.dirs.figsOut,expt.animal);
% 
% if exist(expt.dirs.figsOut)~=7
%     [success]=mkdir(expt.dirs.figsOut);
%     if ~success
%         warning(sprintf('Could not create directory: %s',expt.dirs.figs))
%     end;
% end;

%% standard filenames for derived data

expt.filenames.registered = [exptname,'_green_reg.tif'];
expt.filenames.pca = [exptname,'_green_pca.tif'];
% expt.filenames.smoothed = [exptname,'_green_reg_sm.tif'];
expt.filenames.registeredred = [exptname,'_red_reg.tif'];
expt.filenames.shifts = [exptname,'_shifts.mat'];
expt.filenames.masks = 'masks_neurons.mat';
expt.filenames.masks_auto = 'masks_neurons_auto.mat';
expt.filenames.timecourses = 'timecourses.mat';
expt.filenames.timecoursespca = 'timecourses_pca.mat';
expt.filenames.pca_usv = 'pca_usv.mat';
expt.filenames.ica = 'ica.mat';
expt.filenames.imgs = 'imgs.mat';
expt.filenames.deconv= 'deconv.mat';
expt.filenames.revcorr = 'revcorr.mat';

%% visual stimulation directories
temp = fullfile(expt.dirs.presentation,expt.animal);

% presentation log directory
if exist(temp)==7 
    expt.dirs.presentation = temp;
end

% presentation log file
temp = fullfile(expt.dirs.presentation, [expt.name,'.log']);
if exist(temp)==2
    expt.filenames.presentationlog = temp;
end

return;
