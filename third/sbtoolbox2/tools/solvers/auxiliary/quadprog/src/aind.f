C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  2008, Adapted to be used for a MATLAB MEX function by Henning Schmidt
C  henning@sbtoolbox2.org
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     aind subroutine
      subroutine aind(ind,m,q,n,ok)
      implicit none
      integer m, q, n, i, j, ok
      real*8 ind(m,*)
      ok = 0
      do i=1,q
         if( INT(ind(1,i)) .LT. 1 .OR. INT(ind(1,i)) .GT. n ) return
         do j=2,INT(ind(1,i))+1
            if( INT(ind(j,i)) .LT. 1 .OR. INT(ind(j,i)) .GT. n ) return
         enddo
      enddo
      ok = 1
      return
      end

