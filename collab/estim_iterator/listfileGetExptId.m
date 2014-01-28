function [exptName seriesName] = listfileGetExptid(listStruct, listN)
%
% customize this for your listfiles
%
%$Id: listfileGetExptId.m 191 2008-04-25 05:06:01Z histed $

ls = listStruct;

fNames = fieldnames(ls);
assert(all(ismember({ 'ExptName', 'SeriesEstimNum' }, fNames)), ...
       'Missing fields in listStruct');

exptName = ls.ExptName{listN};
seriesName = sprintf('estim%d', ls.SeriesEstimNum(listN));


