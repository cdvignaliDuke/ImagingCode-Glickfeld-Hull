function out = tcContrast(timeCourses,ind)

[nsamples,ncells]=size(timeCourses);

if nargin < 2
    ind = 1:nsamples;
end

baseline = repmat( mean(timeCourses(ind,:),1) , nsamples,1);

out = (timeCourses - baseline)./baseline*100;

return


