function out = RegMulticore_shifts(source, target, filename, multicoredir, nChunks)
%REGMULTICORE_SHIFTS Stack registration on multiple cores
% OUT = REGMULTICORE(SOURCE, TARGET, FILENAME, MULTICOREDIR, NCHUNKS) where
%  SOURCE is 3d array specifying the stack to be registered, 
%  TARGET is 2d array reference image, 
%  FILENAME is string specifiying output file name
%  MULTICOREDIR string specifiying where to store parameter files, and
%  NCHUNKS is the number of multicore parameters files.
%  OUT = shifts in pixels

% Aaron Kerlin, Vincent Bonin
%
% 08/01/20 vb uses readtiff, writetiff instead of FastLoadStackAK
%          vb Cropping from source image
%          vb replaced call to RegMcoreSlaveAK with call to RegMcoreSlave 
%          vb always writes tif file if filename specified, deleted writeOut flag
%
% 08/05/02 ak generates RegMulticore_shifts.m, RegMcoreSlave_shifts.m, and
%           IJAlign_shifts.class based on vb current working set
%
%Options

%Cropping = '0 0 126 126'; % left- top- right- and bottom-most pixels to consider 
Transformation = '-translation'; % translation is currently the only suppported transformation, speak to AK if you need others.

%End Options

if (ischar(source))
    disp('Reading source stack from disk...'); 
    source = readtiff(source);
    disp('done.');
end    

[h,w,l]=size(source);
Cropping = num2str([0 0 w-1 h-1]);

center=num2str(fix(size(source,1)/2)-1);
landmarks = [center,' ',center,' ',center,' ',center];
cmdstr=['-align -window s ',Cropping,' -window t ', Cropping,' ', Transformation, ' ', landmarks,' -hideOutput']

tic
disp('Starting alignment...');
nSlices = size(source,3);
frames_per_thread = fix(nSlices/nChunks);

parameterCell = cell(1,nChunks);

for i=1:nChunks
        start = (i-1)*frames_per_thread+1;
        if i == nChunks
            stop = nSlices;
        else
            stop = start + frames_per_thread - 1;
        end
        substack=source(:,:,start:stop);
        parameterCell{1,i} = {cmdstr,substack,target};
end

clear source;

maxMasterEvaluations = 0;
out = startmulticoremaster(@RegMcoreSlave_shifts, parameterCell, multicoredir, maxMasterEvaluations);
out = permute(out,[2,1]); % to allow cat along third dim in cell2mat
out = cell2mat(out);
disp('Alignment complete.')
toc

if ~isempty(filename)
   disp('Saving shifts to disk...'); 
   [pathstr]=fileparts(filename);
   if exist(pathstr)~=7
       mkdir(pathstr);
   end
   save(filename, 'out');
   disp('done.');
end

return;
