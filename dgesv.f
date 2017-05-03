c-----------------------------------------------------------------------
c     dgesv inverts a matrix a by LU decomposition and backsubstitution,
c     in order to solve a set of linear equations a x = b.
c     Ported from LH's version, which was probably taken directly from
c     Numerical Recipes....sjs 11/27/07
c-----------------------------------------------------------------------
c     Modified to match flucq code's (netlib's) call interface.
c     ...sjs 1/11/10
c-----------------------------------------------------------------------
c
      subroutine dgesv(n, nrhs, a, nra, ipivot, b, nrb, istatus)

      implicit none

c     n       = order of array a
c     nrhs    = number of right-hand sides, i.e. 2nd args of array b
c     a       = array from a x = b
c     nra     = number of rows (leading dimension) of array a
c     ipivot  = defines permutation of a
c     b       = constant vector b on input, result vector x on output
c     nrb     = number of rows (leading dimension) of array b
c     istatus = success (0 => success, 1 => singular matrix)

      integer n
      integer nrhs
      real*8  a(nra,n)
      integer nra
      integer ipivot(n)
      real*8  b(n)
      integer nrb
      integer istatus

c     check to see that calling routine is not using any of the Netlib call
c     interface variables that are not used here

      if (nrhs .ne. 1) then
         write(0, *) 'dgesv: ERROR! not ready for nrhs = ', nrhs,
     .      ' != 1'
      endif

      if (nrb .lt. n) then
         write(0, *) 'dgesv: ERROR! nrb = ', nrb,
     .      ' is less than n = ', n
      endif

c     LU decomposition

      call ludcmp(n, nra, a, ipivot, istatus)

      if (istatus .ne. 0) then
         return
      endif

c     back-substitution

      call lubksb(n, nra, a, ipivot, b)

      return
      end
c
c-----------------------------------------------------------------------
c     ludcmp performs an LU decomposition of a row-wise permutation of a
c     matrix supplied as input.
c     Ported from LH's version, which was probably taken from Numerical
c     Recipes.
c     Information on parity of row exchanges is no longer kept.
c     ...sjs 11/27/07
c-----------------------------------------------------------------------
c     
      subroutine ludcmp(n, nra, a, ipivot, istatus)

      implicit none

c     a       = n x n matrix (input), LU decomposition (output)
c     nra     = number of rows (leading dimension) of a
c     n       = order of matrix a
c     ipivot  = permutation of rows used in LU decomposition (output)
c     istatus = success of LU decomposition (0 => successful, 1 => singular)

      integer n
      integer nra
      real*8  a(nra,n)
      integer ipivot(n)
      integer istatus

c     local variables

      integer nmax
      parameter (nmax = 10000)

      integer i, imax, j, k
      real*8  scale(nmax),
     .     biggst, recip, score, sum, swap

      if (n .gt. nmax) then
         write(0, *) 'ludcmp: n = ', n, ' rows exceeds nmax = ', nmax
         write(0, *) '        increase nmax and recompile'
         stop
      endif

c     so far so good

      istatus = 0

c     loop over rows, finding the maximum element in each row
c     and using it to calculate a scaling factor

      do 110 i = 1, n
         biggst = 0.d0
         do 100 j = 1, n
            if (abs(a(i,j)) .gt. biggst) then
               biggst = abs(a(i,j))
            endif
 100     continue
         if (biggst .eq. 0.d0) then
c           all elements in row i are zero => matrix is singular
            istatus = 1
            go to 999
         endif
         scale(i) = 1.d0 / biggst
 110  continue

c     Crout's method

c     loop over columns

      do 260 j = 1, n

c        Eq. 2.3.12 without i=j term
         do 210 i = 1, j-1
            sum = a(i,j)
            do 200 k = 1, i-1
               sum = sum - a(i,k) * a(k,j)
 200        continue
            a(i,j) = sum
 210     continue
         
         biggst = 0.d0

c        i=j term in Eq. 2.3.12 and i=j+1 ... n in Eq. 2.3.13
         do 230 i = j, n
            sum = a(i,j)
            do 220 k = 1, j-1
               sum = sum - a(i,k) * a(k,j)
 220        continue
            a(i,j) = sum

c           figure of merit for the pivot

            score = scale(i) * abs(sum)

c           is this better than the best so far?

            if (score .ge. biggst) then
               imax = i
               biggst = score
            endif
 230     continue

         if (j .ne. imax) then

c           interchange rows, and scale factor

            do 240 k = 1, n
               swap = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = swap
 240        continue
            scale(imax) = scale(j)

         endif
         ipivot(j) = imax
         if (a(j,j) .eq. 0.d0) then
c           the matrix is singular
            istatus = 1
            go to 999
         endif

c        finally, divide by the pivot element

         if (j .ne. n) then
            recip = 1.d0 / a(j,j)
            do 250 i = j+1, n
               a(i,j) = a(i,j) * recip
 250        continue
         endif

 260  continue
               
 999  return
      end
c
c-----------------------------------------------------------------------
c     lubksb performs the forward and back substitution to obtain the
c     solution x for a system of linear equations represented by
c     a x = b.
c     Ported from LH's version, which was probably taken from Numerical
c     Recipes....sjs 11/27/07
c-----------------------------------------------------------------------
c
      subroutine lubksb(n, nra, a, ipivot, b)

      implicit none

c     n      = order of matrix a
c     nra    = number of rows (leading dimension) of array a
c     a      = LU decomposition of original matrix a
c     ipivot = permutation of rows used to obtain LU decomposition
c     b      = constant vector b on input, solution x on output

      integer n
      integer nra
      real*8 a(nra,n)
      integer ipivot(n)
      real*8 b(n)

c     local variables

      integer i, ii, j, ll
      real*8  sum

c     Do the forward substitution, Eq. 2.3.6, unscrambling the permutation
c     as we go.

c     This routine takes into account the possibility that b will begin
c     with many zero elements, so it is efficient for use in matrix
c     inversion.

      ii = 0
      do 110 i = 1, n
         ll = ipivot(i)
         sum = b(ll)
         b(ll) = b(i)
         if (ii .ne. 0) then
            do 100 j = ii, i-1
               sum = sum - a(i,j) * b(j)
 100        continue
         else if (sum .ne. 0.d0) then

c           a nonzero element was encountered, so from now on we will
c           have to do the sums in the loop above

            ii = i
         endif
         b(i) = sum
 110  continue

c     Now do the back-substitution, Eq. 2.3.7

      do 210 i = n, 1, -1
         sum = b(i)
         do 200 j = i+1, n
            sum = sum - a(i,j) * b(j)
 200     continue

c        This is one component of the solution vector x

         b(i) = sum / a(i,i)

 210  continue

      return
      end
c
