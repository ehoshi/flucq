c***********************************************************************
c     tfunc calculates the value of Student's t-function (two-sided), 
c     for a given
c     t and n.  Taken from Numerical Recipes....sjs 5/20/93
c***********************************************************************
c
      function tfunc(t, n)
c
      include 'implic'
c
      tfunc = 1 - betinc(0.5d0 * n, 0.5d0, n / (n + t ** 2))
      return
      end
c
c***********************************************************************
c     tconf calculates the t value for which tfunc, with n degrees of
c     freedom, would return a value of conf / 2.  This is the t value
c     needed to construct a (2-sided) 100*conf% confidence interval.
c     ...sjs 8/19/93
c***********************************************************************
c
      function tconf(conf, n)
c
      include 'implic'
c
      ttol = 0.01
      tconf = 2.00d0
      hconf = 0.5d0 + 0.5d0 * conf
      temp = tfunc(tconf, n)
      if (abs(temp - hconf) .lt. ttol) then
         return
      else if (temp .gt. hconf) then
         tconhi = tconf
         tconlo = tconf
  200    tconlo = 0.5d0 * tconlo
         temp = tfunc(tconlo, n)
         if (abs(temp - hconf) .lt. ttol) then
            return
         else if (temp .gt. hconf) then
            go to 200
         endif
      else
         tconlo = tconf
         tconhi = tconf
  300    tconhi = 2.0d0 * tconhi
         temp = tfunc(tconhi, n)
         if (abs(temp - hconf) .lt. ttol) then
            return
         else if (temp .lt. hconf) then
            go to 300
         endif
      endif
  400 tconmi = (tconhi - tconlo) * 0.5d0
      temp = tfunc(tconmi, n)
      if (abs(temp - hconf) .lt. ttol) then
         return
      else if (temp .gt. hconf) then
         tconhi = tconmi
      else
         tconlo = tconmi
      endif
      go to 400
      end
c      
c***********************************************************************
c     betinc calculates the incomplete beta function (beta function of 
c     a and b, but integrated only up to x, not 1).  Taken from
c     Numerical Recipes....sjs 5/20/93
c***********************************************************************
c
      function betinc(a, b, x)
c
      include 'implic'
c
      if (x .lt. 0 .or. x .gt. 1) stop 'betinc: bad x argument'
      if (x .eq. 0 .or. x .eq. 1) then
         temp = 0
      else
         temp = exp(gammln(a + b) - gammln(a) - gammln(b) +
     $        a * dlog(x) + b * dlog(1 - x))
      endif
      if (x .lt. (a + 1) / (a + b + 2)) then
         betinc = temp * betacf(a, b, x) / a
      else
         betinc = 1 - temp * betacf(b, a, 1 - x) / b
      endif
      return
      end
c
c***********************************************************************
c     betacf evaluates a particular continued fraction (by iteration)
c     that is needed in the betinc routine.  Taken from Numerical
c     Recipes....sjs 5/20/93
c***********************************************************************
c
      function betacf(a, b, x)
c
      include 'implic'
c
      parameter (itmax = 100, eps = 3.0d-7)
c
      aplusb = a + b
      aplus1 = a + 1
      asub1 = a - 1
      Aold = 1
      Bold = 1
      Acurr = 1
      Bcurr = 1 - aplusb * x / aplus1
      do 200 iit = 1, itmax
         rm = iit
         twom = rm + rm
         d = rm * (b - m) * x / ((asub1 + twom) * (a + twom))
         Anext = Acurr + d * Aold
         Bnext = Bcurr + d * Bold
         d = -(a + rm) * (aplusb + rm) * x / 
     $        ((a + twom) * (aplus1 + twom))
         Anext2 = Anext + d * Acurr
         Bnext2 = Bnext + d * Bcurr
         Asave = Acurr
         Aold = Anext / Bnext2
         Bold = Bnext / Bnext2
         Acurr = Anext2 / Bnext2
         Bcurr = 1
         if (dabs(Acurr - Asave) .lt. eps * dabs(Acurr)) go to 300
  200 continue
      stop 'betacf: iterations exceeded itmax'
  300 betacf = Acurr
      return
      end
c
c***********************************************************************
c     gammln evaluates the natural log of the gamma function of the
c     argument passed to it.  Taken from Numerical Recipes.
c     ...sjs 5/20/93
c***********************************************************************
c
      function gammln(x)
c
      include 'implic'
c
      real*8 cof(6)
c
      data cof/76.18009173d0, -86.50532033d0, 24.01409822d0, 
     $     -1.231739516d0, 0.120858003d-2, -0.536382d-5/
      data rt2pi/2.50662827465d0/
      data half, fpf/0.5d0, 5.5d0/
c
      z = x - 1.0d0
      gammln = z + 5.5d0
      gammln = (z + 0.5d0) * dlog(gammln) - gammln
      gammln = gammln + dlog((1.0d0 + 76.18009173d0 / (z + 1.0d0) -
     $     86.50532033d0 / (z + 2.0d0) + 24.01409822d0 / (z + 3.0d0) -
     $     1.231739516d0 / (z + 4.0d0) + 0.120858003d-2 / (z + 5.0d0) -
     $     0.536382d-5 / (z + 6.0d0)) * rt2pi)
      return
      end
c
c***********************************************************************
c     herfc evaluates the complementary error function of its argument,
c     using an approximate formula with advertised fractional error of
c     1.2e-7 or less everywhere.  Taken from Numerical Recipes.
c     ...sjs 8/25/93
c***********************************************************************
c
      function herfc(x)
c
      include 'implic'
c
      z = abs(x)
      t = 1.0d0 / (1.0d0 + 0.5d0 * z)
      herfc = t * exp(-z * z - 1.26551223 + t * (1.00002368 + t * 
     $     (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t *
     $     (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t *
     $     (-0.82215223 + t * 0.17087277)))))))))
      if (x .lt. 0.0d0) herfc = 2.0d0 - herfc
      return
      end
c***********************************************************************
c     fit will fit a set of ndata data points stored in x(), y() with a
c     straight line y = a + b * x by minimizing chi^2.  The std devs of
c     the data are stored in sig().  a and b are returned, along with
c     siga and sigb, their respective probably uncertainties.  chi2 is
c     also returned, as is q, the probability that a perfect fit (with
c     correct model and normal errors) would have this chi2 or better.
c     If  mwt = 0 on input, then the std devs are assumed to be
c     unavailable, q is returned as 1.0, and the normalization of chi2
c     is to unit std dev on all pts.  Taken from Numerical Recipes.  The
c     formulae aren't the standard ones - S_xx and S_xy don't show up.
c     These are apparently less susceptible to roundoff errors.
c     ...sjs 11/30/93
c***********************************************************************
c
      subroutine fit(x, y, ndata, sig, mwt, a, b, siga, sigb, chi2, q)
c
      include 'implic'
c
      real*8 x(ndata), y(ndata), sig(ndata)
c
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
      sx = 0.0d0
      sy = 0.0d0
      st2 = 0.0d0
      b = 0.0d0
      chi2 = 0.0d0
c***********************************************************************
c     Accumulate the weighted or unweighted sums of 1, x, and y data:
c***********************************************************************
      if (mwt .ne. 0) then
         ss = 0.0d0
         do 200 i = 1, ndata
            wt = 1.0d0 / (sig(i) * sig(i))
            ss = ss + wt
            sx = sx + x(i) * wt
            sy = sy + y(i) * wt
  200    continue
      else
         do 210 i = 1, ndata
            sx = sx + x(i)
            sy = sy + y(i)
  210    continue
         ss = dble(ndata)
      endif
      sxoss = sx / ss
c***********************************************************************
c     Start calculating b, and some factors to be used below:
c***********************************************************************
      if (mwt .ne. 0) then
         do 300 i = 1, ndata
            t = (x(i) - sxoss) / sig(i)
            st2 = st2 + t * t
            b = b + t * y(i) / sig(i)
  300    continue
      else
         do 310 i = 1, ndata
            t = x(i) - sxoss
            st2 = st2 + t * t
            b = b + t * y(i)
  310    continue
      endif
c***********************************************************************
c     Get a, b, siga and sigb:
c***********************************************************************
      b = b / st2
      a = (sy - sx * b) / ss
      siga = sqrt((1.0d0 + sx * sx / (ss * st2)) / ss)
      sigb = sqrt(1.0d0 / st2)
c***********************************************************************
c     Calculate chi2:
c***********************************************************************
      if (mwt .eq. 0) then
         do 500 i = 1, ndata
            dy = y(i) - a - b * x(i)
            chi2 = chi2 + dy * dy
  500    continue
         q = 1.0d0
         sigdat = sqrt(chi2 / (ndata - 2))
         siga = siga * sigdat
         sigb = sigb * sigdat
      else
         do 510 i = 1, ndata
            wdy = (y(i) - a - b * x(i)) / sig(i)
            chi2 = chi2 + wdy * wdy
  510    continue
         q = gammq(0.5d0 * (ndata - 2), 0.5d0 * chi2)
      endif
      return
      end
c
c***********************************************************************
c     gammq returns the imcomplete gamma function q(a,x) = 1 - p(a,x).
c     Depending on the relative sizes of a and x, either a series or a
c     continued fraction equation is used to calculate q.  Taken from 
c     Numerical Recipes....sjs 11/30/93
c***********************************************************************
c
      function gammq(a, x)
c
      include 'implic'
c
      if (x .lt. 0.0d0 .or. a .le. 0.0d0) stop 'gammq: invalid args'
      if (x .lt. a + 1.0d0) then
         call gser(gamser, a, x, gln)
         gammq = 1.0d0 - gamser
      else
         call gcf(gammcf, a, x, gln)
         gammq = gammcf
      endif
      return
      end
c
c***********************************************************************
c     gser returns the incomplete gamma function gamser = p(a, x) as 
c     evaluated by its series representation.  gln is the ln of the
c     gamma fxn of a....sjs 11/30/93
c***********************************************************************
c
      subroutine gser(gamser, a, x, gln)
c
      include 'implic'
c
      parameter (itmax = 100, eps = 3.0d-7)
c
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
      gln = gammln(a)
c***********************************************************************
c     Take care of errors and a trivial case:
c***********************************************************************
      if (x .le. 0.0d0) then
         if (x .lt. 0.0d0) stop 'gser: invalid arg'
         gamser = 0.0d0
         return
      endif
c***********************************************************************
c     Initialize some more stuff:
c***********************************************************************
      ap = a
      sum = 1.0d0 / a
      del = sum
c***********************************************************************
c     Iterate a series until it converges or itmax is hit:
c***********************************************************************
      do 400 n = 1, itmax
         ap = ap + 1.0d0
         del = del * x / ap
         sum = sum + del
         if (abs(del) .lt. abs(sum) * eps) go to 500
  400 continue
      stop 'gser: itmax was hit without converging'
c***********************************************************************
c     Finish up:
c***********************************************************************
  500 gamser = sum * exp(-x * a * log(x) - gln)
      return
      end
c
c***********************************************************************
c     gcf calculates the incomplete gamma fxn gammcf = q(a, x) as
c     evaluated by a continued fraction representation.  gln is the
c     ln of the gamma fxn of a.  Taken from Numerical Recipes.
c     ...sjs 11/30/93
c***********************************************************************
c
      subroutine gcf(gammcf, a, x, gln)
c
      include 'implic'
c
      parameter (itmax = 100, eps = 3.0d-7)
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
      gln = gammln(a)
      gold = 0.0d0
      a0 = 1.0d0
      a1 = x
      b0 = 0.0d0
      b1 = 1.0d0
      fac = 1.0d0
c***********************************************************************
c     Iterate the continued fraction until convergence or itmax:
c***********************************************************************
      do 200 n = 1, itmax
         an = dble(n)
         ana = an - a
         a0 = (a1 + a0 * ana) * fac
         b0 = (b1 + b0 * ana) * fac
         anf = an * fac
         a1 = x * a0 + anf * a1
         b1 = x * b0 + anf * b1
         if (a1 .ne. 0.0d0) then
            fac = 1.0d0 / a1
            g = b1 * fac
            if (abs((g - gold) / g) .lt. eps) go to 300
            gold = g
         endif
  200 continue
      stop 'gcf: hit itmax without converging'
c***********************************************************************
c     Finish up:
c***********************************************************************
  300 gammcf = exp(-x + a * log(x) - gln) * g
      return
      end
c
c***********************************************************************
c     ceil returns (as a real*8) the rounded-up value of the (real*8) 
c     argument....sjs 4/13/95
c***********************************************************************
c
      function ceil(arg)
c
      include 'implic'
c
      if (arg - int(arg) .eq. 0.d0) then
         ceil = arg
      else
         ceil = dble(int(arg+1.d0))
      endif
      return
      end
c
c***********************************************************************
c     dcsint is a klugey replacement for the ESSL routine, which was
c     used for cubic spline interpolation.  This routine merely does
c     linear interpolation....sjs 2/18/99
c***********************************************************************
c
      subroutine dcsint(xvec, ygrid, coeffs, nbins, initsp, 
     .     xval, yval, njunk)
c
      include 'implic'
      include 'qpar'
c
      real*8 xvec(numbin), ygrid(numbin), coeffs(numbin)
c
      if (initsp .eq. 3) then
         write(iuout, *) 'dcsint:  WARNING! ',
     .        'cubic spline routine will use linear interpolation'
         return
      else if (initsp .eq. 4) then
         deltx = dlimr / (nbins - 1)
         xbin = xval / deltx + 1
         jgrid = int(xbin)
         kgrid = jgrid + 1
         yval = ygrid(jgrid) + (xval - (jgrid - 1) * deltx) /
     .        deltx * (ygrid(kgrid) - ygrid(jgrid))
         return
      endif
      return
      end
c

