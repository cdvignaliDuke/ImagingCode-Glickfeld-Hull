% This small script solves the following quadratic programming problem
% using the quadprogSB and the quadprogcompactSB solver:
%
%      min ( x1^2 + x2^2 + x3^2 - 5*x2 ) 
%    [x1,x2,x3]'
%
%  subject to:
%
%    4*x1+3*x2 < 8
%   -2*x1-1*x2 < -2
%    2*x2-1*x3 < 0      
%
% The standard form is:
%
%     min(-d'*x + 1/2*x'*D*x)  
%      x
%
%     subject to 
%                  A x <= b

clc; clear all;

%% Running the standard solver quadprogSB
D    = eye(3);
d    = [0,-5,0];
A    = [4, 3, 0;
       -2,-1, 0;
        0, 2,-1];
b    = [8,
       -2,
        0];
output = quadprogSB(D, d, A, b);
disp('optimal solution x:')
x = output.sol'
disp('optimal value:')
optval = output.optval
disp('indices of active constraints:');
output.actconst
disp('A*x at optimal solution x:');
[4*x(1)+3*x(2),  -2*x(1)-1*x(2),  2*x(2)-1*x(3)]
disp('Checking constraints manually (1 means active, 0 means inactive):');
[4*x(1)+3*x(2) == 8,  -2*x(1)-1*x(2) == -2,  2*x(2)-1*x(3) == 0]

%% Running the compact solver quadprogSB
% In case of a huge sparsely populated A matrix it is advantageous to use
% the compact solver quadprogcompactSB instead. The matrix A is then
% defined by the two matrices Amat and Aind. Their definition can be found
% in the help text of the quadprogcompactSB function.
% The function convertAquadprogSB converts an A into this format.
%
% In the following example we use the same settings as above.
[Amat,Aind] = convertAquadprogSB(A);
output = quadprogcompactSB(D, d, Amat, Aind, b);
% Identical result:
output.sol
output.optval

