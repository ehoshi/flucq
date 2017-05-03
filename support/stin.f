      program stin
c***********************************************************************
c     stin calculates the statistical inefficiency for a scalar 
c     parameter stored in a single column in a file....sjs 4/21/94
c***********************************************************************
c
      include 'implic'
      include 'genpar'
      include 'anapar'
c
      character*80 datfil(mdatfl), filnam, stnfil, corfil
      integer      ibsize(mbatch)
c$$$      logical      consfl
      real*8       window(maxwin),
     $             sum(mbatch), ssqbs(mbatch), 
     $             stinef(mbatch), corrf(mbatch),
     $             propa(mdatpt)
c
c***********************************************************************
c     initialize some stuff:
c***********************************************************************
      bigsum = 0.0d0
      do 100 ibatch = 1, mbatch
         sum(ibatch) = 0.0d0
         ssqbs(ibatch) = 0.0d0
 100  continue
c***********************************************************************
c     read input file stin.in:
c***********************************************************************
      write(filnam, *) 'stin.in'
      open(iutmp, file = filnam, err = 9990, iostat = ioval)
      read(iutmp, *, err = 9980, end = 9970, iostat = ioval) ndatfl
      if (ndatfl .gt. mdatfl) then
         write(iuout, *) 'stin: maximum number of data files (mdatfl) ',
     $        'exceeded'
         stop
      else if (ndatfl .lt. 1) then
         write(iuout, *) 'stin: at least one data file must be read'
         stop
      endif
c$$$      if (ndatfl .gt. 1)
c$$$     $     read(iutmp, *, err = 9980, end = 9970, iostat = ioval) consfl
      do 200 idatfl = 1, ndatfl
         read(iutmp, *, err = 9980, end = 9970, iostat = ioval) 
     $        datfil(idatfl)
 200  continue
      read(iutmp, *, err = 9980, end = 9970, iostat = ioval) stnfil
      read(iutmp, *, err = 9980, end = 9970, iostat = ioval) corfil
      read(iutmp, *, err = 9980, end = 9970, iostat = ioval) dt
      read(iutmp, *, err = 9980, end = 9970, iostat = ioval) icol
      read(iutmp, *, err = 9980, end = 9980, iostat = ioval) nwin
      read(iutmp, *, err = 9980, end = 9970, iostat = ioval) nbset
      ibatch = 0
      nbatch = 0
      do 220 ibset = 1, nbset
         read(iutmp, *, err = 9980, end = 9970, iostat = ioval) ibtchn
         read(iutmp, *, err = 9980, end = 9970, iostat = ioval) ibtchi
         read(iutmp, *, err = 9980, end = 9970, iostat = ioval) ibtchs
         nbatch = nbatch + ibtchn
         if (nbatch .gt. mbatch) then
            write(6, *) 'stin: maximum number of batch sizes ',
     $           '(mbatch) exceeded'
            stop
         endif
         do 210 iibtch = 1, ibtchn
            ibatch = ibatch + 1
            ibsize(ibatch) = ibtchs + iibtch * ibtchi
 210     continue
 220  continue
      if (nwin .gt. maxwin) then
         write(iuout, *) 'stin: too large a moving average window'
         stop
      endif
      if (ibsize(1) .ne. 1) then
         write(6, *) 'stin: first batch size must be 1'
         stop
      endif
      close(iutmp)
      write(iuout, '(1x,a,i4,a)') 'There are ', nbatch, ' batch sizes'
      write(iuout, '(1x,a,i2,a,i2,a)') 'Reading data from column ',
     $     icol, ' of ', ndatfl, ' data files.'
c***********************************************************************
c     read numbers from the specified column of the specified files,
c     accumulating them into sum().  When a batch is full, shuffle the
c     square of that sum into ssqbs().
c***********************************************************************
      icount = 0
      do 330 idatfl = 1, ndatfl
         write(filnam, *) datfil(idatfl)
         open(iutmp, file = filnam, err = 9990, iostat = ioval)
         write(iuout, *) 'Reading from ', 
     $        datfil(idatfl)(1:index(datfil(idatfl), ' ')), '...'
 300     continue
         read(iutmp, *, err = 9980, end = 320, iostat = ioval) 
     $        (rjunk, i = 1, icol - 1), prop
         icount = icount + 1
         propa(icount) = prop
         prop2 = prop * prop
         bigsum = bigsum + prop
         do 310 ibatch = 1, nbatch
            sum(ibatch) = sum(ibatch) + prop
            if (mod(icount, ibsize(ibatch)) .eq. 0) then
               ssqbs(ibatch) = ssqbs(ibatch) + sum(ibatch) * sum(ibatch)
               sum(ibatch) = 0.0d0
            endif
 310     continue
         go to 300
 320     close(iutmp)
 330  continue
      write(iuout, *) 'found ', icount, ' data points'
c***********************************************************************
c     Scale down ssqbs so that it is a 
c     mean of the squares of the batch means, use this to get the batch
c     variances, and write the statistical inefficiency data out to the
c     output file (stnfil).  Square root and inverse of the batch sizes
c     are written out for use in plotting the data.
c***********************************************************************
      write(iuout, *) 'calculating statistical inefficiency...'
      write(iuout, '(a,g10.3,a)') ' data for t > ', icount / 10 * dt, 
     $     ' will be too noisy to use'
      write(iuout, '(a,g10.3)') 
     $     ' max. believable s value = ', icount / 100 * dt
      write(filnam, *) stnfil
      open(iutmp, file = filnam, err = 9990, iostat = ioval)
      rnsinv = 1.0d0 / icount
      stinef(1) = rnsinv * (ssqbs(1) - rnsinv * bigsum * bigsum)
c      write(iutmp, *) dt, sqrt(dt), 1.d0 / dt, dt, dt/nwin
c      window(1) = dt
c      do 380 iwin = 2, nwin
c         window(iwin) = 0.d0
c 380  continue
c      rmovav = dt / nwin
      do 400 ibatch = 2, nbatch
         if (icount .ge. 2 * ibsize(ibatch)) then
            rbsinv = 1.0d0 / ibsize(ibatch)
            rnbinv = 1.0d0 / (icount / ibsize(ibatch))
            prtsum = bigsum - sum(ibatch)
            stinef(ibatch) = rnbinv * rbsinv *
     $           (ssqbs(ibatch) - rnbinv * prtsum * prtsum) / stinef(1)
c            rmovav = rmovav + 
c     $           (stinef(ibatch) * dt - window(nwin)) / nwin
c            do 390 iwin = nwin, 2, -1
c               window(iwin) = window(iwin - 1)
c 390        continue
c            window(1) = stinef(ibatch) * dt
c            write(iutmp, *) ibsize(ibatch) * dt, 
c     $           sqrt(dt * ibsize(ibatch)),
c     $           rbsinv / dt, stinef(ibatch) * dt, rmovav
         endif
 400  continue
      stinef(1) = 1.d0
      do 500 ibatch = 1, nbatch
         winwid =2.d0 * sqrt(dble(ibsize(ibatch)))
         iwinlo = ibatch - winwid / 2
         if (iwinlo .lt. 1) then
            iwinlo = 1
         endif
 410     continue
         if (ibsize(iwinlo) + winwid / 2 .lt. ibsize(ibatch)) then
            iwinlo = iwinlo + 1
            go to 410
         endif
         iwinhi = ibatch + winwid / 2
         if (iwinhi .gt. nbatch) then
            iwinhi = nbatch
         endif
 420     continue
         if (ibsize(iwinhi) - winwid / 2 .gt. ibsize(ibatch)) then
            iwinhi = iwinhi - 1
            go to 420
         endif
         rmovav = 0.d0
         do 430 iwin = iwinlo, iwinhi
            rmovav = rmovav + stinef(iwin)
 430     continue
         rmovav = rmovav / (iwinhi - iwinlo + 1)
         write(iutmp, *) ibsize(ibatch) * dt,
     $        sqrt(dt * ibsize(ibatch)),
     $        1.d0 / (dt * ibsize(ibatch)), stinef(ibatch) * dt, 
     $        rmovav * dt
 500  continue
      close(iutmp)
c***********************************************************************
c     some correlation fxn stuff:
c***********************************************************************
      write(iuout, *) 'calculating autocorrelation function...'
      do 610 ipt = 1, icount
         prop = propa(ipt) - bigsum / icount
         iptmax = min(icount,ipt + nbatch)
         corrf0 = corrf0 + prop * prop
         do 600 ipt2 = ipt + 1, iptmax
            idiff = ipt2 - ipt
            corrf(idiff) = corrf(idiff) + prop * 
     $           (propa(ipt2) - bigsum / icount)
 600     continue
 610  continue
      write(filnam, *) corfil
      open(iutmp, file = filnam, err = 9990, iostat = ioval)
      corrf0 = corrf0 / icount
      write(iutmp, *) 0.d0, 1.d0
      do 700 idiff = 1, nbatch
         corrf(idiff) = corrf(idiff) / corrf0 / (icount - idiff)
         write(iutmp, *) idiff * dt, corrf(idiff)
 700  continue
      close(iutmp)
      stop
c***********************************************************************
c     error handling:
c***********************************************************************
 9970 write(iuout, *) 'stin: end of file reached in ', filnam
      call ioerr(ioval)
      stop
 9980 write(iuout, *) 'stin: error writing to ', filnam
      call ioerr(ioval)
      stop
 9990 write(iuout, *) 'stin: error reading from ', filnam
      call ioerr(ioval)
      stop
      end



