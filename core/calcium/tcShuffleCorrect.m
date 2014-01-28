function [out,y] = tcShuffleCorrect (timeCourses, nframes_per_repeat)
%TCSHUFFLECORRECT Subtracts trial average from time courses
% [out,y] = tcShuffleCorrect (timeCourses, nframes_per_repeat)

[nframes,ncells]=size(timeCourses);
ntrials=nframes./nframes_per_repeat;

out = zeros(size(timeCourses));
temp = zeros(size(timeCourses));

for icell = 1:ncells
    thisTimeCourse=reshape(timeCourses(:,icell),nframes_per_repeat,ntrials);
    averageTimeCourse=sum(thisTimeCourse,2)./ntrials;
    temp = repmat(averageTimeCourse,ntrials,1);
    out(:,icell)=timeCourses(:,icell)-temp;
end

return;
