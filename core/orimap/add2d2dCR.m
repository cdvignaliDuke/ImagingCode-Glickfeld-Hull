function Aout = add2d2dCR(A, Bx)

% add2d.m from add2d.m***** 
% Clay Reid version of TP and BS program.  sums over 2d boxes
% Function that adds over Bx columns in an image array. 
%
% Function that bins by averaging over Bx (columns) elements along the x-axis.  A is a 2D intensity matrix.
% Bx must be divisible into the array dimensions.
% 
% By is the Binning factor for rows (lines)
% Function Form: add2d(A,Bx)
%
% Written By: Thomas Pologruto
% Cold Spring Harbor Labs
% January 26, 2001.
%
% Edited By: Bernardo Sabatini
% January 26, 2001
% Cold Spring Harbor Labs
	
if Bx == 1
	Aout = A;
else
	Aout = add2d(A, Bx);
	Aout = add2d(Aout',Bx);
	Aout = Aout';

end
	
	
	
	

	