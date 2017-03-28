function ncells = howManyCells(rc,expt)

n = zeros(1,size(expt,2));

for iexp = 1:size(expt,2)

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
down = 10;

fnbx = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate,'final_mask.mat');

load(fnbx)

nc = length(unique(mask_cell(:)))-1;

n(iexp) = nc;


end
ncells = n;
end