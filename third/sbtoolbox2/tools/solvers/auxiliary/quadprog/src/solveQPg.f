C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C MEX INTERFACE FOR THE solveQP FUNCTION
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calling from MATLAB:
C [sol,soluc,crval,iact,nact,iter,ierr] = solveQPcompact(Dmat,dvec,n,n,Amat,bvec,n,q,meq,factorized)
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

C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     The gateway routine
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      integer plhs(*), prhs(*)
      integer nlhs, nrhs

      integer Dmat_pr, dvec_pr, n1_pr, n2_pr, Amat_pr
      integer bvec_pr, n3_pr, q_pr, meq_pr, factorized_pr

      integer n1, n2, n3, q, meq, factorized
      real*8 n1r, n2r, n3r, qr, meqr, factorizedr

      integer soluc_pr, sol_pr, crval_pr, iact_pr, nact_pr
      integer iter_pr, ierr_pr
      real*8 crvalr, nactr

      integer work, work_pr

C CCCCCCCCCCCCCCCCCCC
C     Check for proper number of arguments. 
C CCCCCCCCCCCCCCCCCCC
      if(nrhs .ne. 10) then
         call mexErrMsgTxt('Incorrect number of input arguments.')
      elseif(nlhs .ne. 7) then
         call mexErrMsgTxt('Incorrect number of output arguments.')
      endif

C CCCCCCCCCCCCCCCCCCC
C     Get input arguments
C CCCCCCCCCCCCCCCCCCC
C     Get pointers
      Dmat_pr = mxGetPr(prhs(1))
      dvec_pr = mxGetPr(prhs(2))
      n1_pr = mxGetPr(prhs(3))
      n2_pr = mxGetPr(prhs(4))
      Amat_pr = mxGetPr(prhs(5))
      bvec_pr = mxGetPr(prhs(6))
      n3_pr = mxGetPr(prhs(7))
      q_pr = mxGetPr(prhs(8))
      meq_pr = mxGetPr(prhs(9))
      factorized_pr = mxGetPr(prhs(10))
C     Get integers
      call mxCopyPtrToReal8(n1_pr, n1r, 1)
      n1 = INT(n1r)
      call mxCopyPtrToReal8(n2_pr, n2r, 1)
      n2 = INT(n2r)
      call mxCopyPtrToReal8(n3_pr, n3r, 1)
      n3 = INT(n3r)
      call mxCopyPtrToReal8(q_pr, qr, 1)
      q = INT(qr)
      call mxCopyPtrToReal8(meq_pr, meqr, 1)
      meq = INT(meqr)
      call mxCopyPtrToReal8(factorized_pr, factorizedr, 1)
      factorized = INT(factorizedr)

C CCCCCCCCCCCCCCCCCCC
C     Create output elements
C CCCCCCCCCCCCCCCCCCC
      plhs(1) = mxCreateDoubleMatrix(1, n1, 0)
      sol_pr = mxGetPr(plhs(1))
      plhs(2) = mxCreateDoubleMatrix(1,n1,0);
      soluc_pr = mxGetPr(plhs(2))
      plhs(3) = mxCreateDoubleMatrix(1, 1, 0)
      crval_pr = mxGetPr(plhs(3))
      call mxCopyPtrToReal8(crval_pr, crvalr, 1)
      plhs(4) = mxCreateDoubleMatrix(1, q, 0)
      iact_pr = mxGetPr(plhs(4))
      plhs(5) = mxCreateDoubleMatrix(1, 1, 0)
      nact_pr = mxGetPr(plhs(5))
      plhs(6) = mxCreateDoubleMatrix(1, 2, 0)
      iter_pr = mxGetPr(plhs(6))
      plhs(7) = mxCreateDoubleMatrix(1, 1, 0)
      ierr_pr = mxGetPr(plhs(7))

C CCCCCCCCCCCCCCCCCCC
C     Get memory for the work vector
C CCCCCCCCCCCCCCCCCCC
      work = mxCreateDoubleMatrix(1, 
     *     2*n1+min(n1,q)*(min(n1,q)+5)/2+2*q+1, 0)
      work_pr = mxGetPr(work)

C CCCCCCCCCCCCCCCCCCC
C     Call the computational subroutine.
C CCCCCCCCCCCCCCCCCCC
      call qpgen2(%VAL(Dmat_pr), %VAL(dvec_pr), n1, n2, %VAL(sol_pr),
     *     crvalr, %VAL(Amat_pr),
     *     %VAL(bvec_pr), n3, q, meq, %VAL(iact_pr), nact, 
     *     %VAL(iter_pr), %VAL(work_pr), factorized)  

C CCCCCCCCCCCCCCCCCCC
C     Convert scalars back to doubles
C CCCCCCCCCCCCCCCCCCC
      nactr = REAL(nact)
      factorizedr = REAL(factorized)

C CCCCCCCCCCCCCCCCCCC
C     Return the values to MATLAB
C CCCCCCCCCCCCCCCCCCC
      call mxCopyReal8ToPtr(crvalr, crval_pr, 1)
      call mxCopyReal8ToPtr(nactr, nact_pr, 1)
C   dvec now contains the unconstrained solution which needs to be returned
      call mxCopyReal8ToPtr(%VAL(dvec_pr), soluc_pr, n1)
C   factorized now contains the ierr output which needs to be returned
      call mxCopyReal8ToPtr(factorizedr, ierr_pr, 1)

      return
      end