stack = [];
for iRun = 1:length(userun);
    substack = readtiff(fullfile(base,mouse,date, [date '_' mouse '_run' num2str(userun(iRun)) '.tif']));
    stack = cat(3, stack, substack);
end
clear('substack');
