% buildQuadProgUNIX: Builds the quadprog mex files on Unix systems
%
% Please note that usually on Unix the MEX files do have to be rebuild.

% Information:
% ============
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org
% GNU GPL >= 2

disp('Building quadprog MEX functions.');

% build the MEX files using the mexfSB function included in the SBTOOLBOX2
cd src
mex aind.f aindg.f
mex solveQP.f solveQPg.f util.f
mexfSB solveQPcompact.f solveQPcompactg.f util.f

% clean up
eval(sprintf('!mv aind ../mex/aind.%s',mexext));
eval(sprintf('!mv solveQP ../mex/solveQP.%s',mexext));
eval(sprintf('!mv solveQPcompact  ../mex/solveQPcompact.%s',mexext));
cd ..

