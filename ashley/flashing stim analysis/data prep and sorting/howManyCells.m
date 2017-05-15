function ncells = howManyCells(rc,expt,varargin)

if ~isempty(varargin)
    doDendrites = varargin{1};
else
    doDendrites = 0;
end
    

n = zeros(1,size(expt,2));

for iexp = 1:size(expt,2)

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
down = 10;

if doDendrites
    fnbx = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate,'dendrite_mask.mat');
else
    fnbx = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate,'final_mask.mat');
end

load(fnbx)

nc = length(unique(mask_cell(:)))-1;

n(iexp) = nc;


end
ncells = n;
end