C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C MEX INTERFACE FOR THE aind FUNCTION
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calling from MATLAB:
C [Aindok] = aind(Aind, anrow+1, q, n)
C
C Copyright (C) 2008 Henning Schmidt, henning@sbtoolbox2.org
C
C This program is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published by
C the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C 
C This program is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C 
C You should have received a copy of the GNU General Public License
C along with this program; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
C USA.
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     The gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      integer plhs(*), prhs(*)
      integer nlhs, nrhs
      integer ind_pr, m_pr, q_pr, n_pr, ok_pr
      integer m, q, n, ok
      real*8 mr, qr, nr, okr

C     Check for proper number of arguments. 
      if(nrhs .ne. 4) then
         call mexErrMsgTxt('Incorrect number of input arguments.')
      elseif(nlhs .ne. 1) then
         call mexErrMsgTxt('Incorrect number of input arguments.')
      endif

C     Get input arguments
      ind_pr = mxGetPr(prhs(1))
      m_pr = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(m_pr, mr, 1)
      m = INT(mr)
      q_pr = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(q_pr, qr, 1)
      q = INT(qr)
      n_pr = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(n_pr, nr, 1)
      n = INT(nr)

C     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(1, 1, 0)
      ok_pr = mxGetPr(plhs(1))
      
C     Call the computational subroutine.
      call aind(%VAL(ind_pr),m,q,n,ok)
      okr = REAL(ok)
 
C     Load the data into ok_pr, which is the output to MATLAB.
      call mxCopyReal8ToPtr(okr, ok_pr, 1)

      return
      end