      program tcorr
c***********************************************************************
c     tcorr calculates a time correlation function for a scalar
c     parameter stored in a single column in a file....sjs 9/26/94
c***********************************************************************
c
      include 'implic'
c
      parameter (iuout = 6)
      parameter (iutmp = 7)
      parameter (maxrec = 10000)
      parameter (mbatch = 10000)
      parameter (mdatfl = 20)
c
      character*80 datfil(mdatfl), filnam, outfil
      integer      ibsize(mbatch)
      real*8       data(maxrec)
c
c***********************************************************************
c     initialize some stuff:
c***********************************************************************
      ssq = 0.0
c***********************************************************************
c     read input file tcorr.in:
c***********************************************************************
      write(filnam, *) 'tcorr.in'
      open(iutmp, file = filnam, err = 9990, iostat = ioval)
      read(iutmp, *, err = 9980, end = 9970, iostat = ioval) ndatfl
      if (ndatfl .gt. mdatfl) then
         write(iuout, *) 'tcorr: maximum number of data files ',
     $        '(mdatfl) exceeded'
         stop
      else if (ndatfl .lt. 1) then
         write(iuout, *) 'tcorr: at least one data file must be read'
         stop
      endif
      do 200 idatfl = 1, ndatfl
         read(iutmp, *, err = 9980, end = 9970, iostat = ioval) 
     $        datfil(idatfl)
 200  continue
      read(iutmp, *, err = 9980, end = 9970, iostat = ioval) outfil
      read(iutmp, *, err = 9980, end = 9970, iostat = ioval) dt
      read(iutmp, *, err = 9980, end = 9970, iostat = ioval) icol
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
      if (ibsize(1) .ne. 1) then
         write(6, *) 'tcorr: first batch size must be 1'
         stop
      endif
      close(iutmp)
      write(iuout, '(1x,a,i4,a)') 'There are ', nbatch, ' batch sizes'
      write(iuout, '(1x,a,i2,a,i2,a)') 'Reading data from column ',
     $     icol, ' of ', ndatfl, ' data files.'
c***********************************************************************
c     read numbers from the specified column of the specified files,
c     storing it in the data() array.
c***********************************************************************
      icount = 0
      rmean = 0
      do 330 idatfl = 1, ndatfl
         write(filnam, *) datfil(idatfl)
         open(iutmp, file = filnam, err = 9990, iostat = ioval)
         write(iuout, *) 'Reading from ',
     $        datfil(idatfl)(1:index(datfil(idatfl), ' ')), '...'
 300     continue
         icount = icount + 1
         if (icount .gt. maxrec) then
            write(iuout, *) 'tcorr: maximum number of data points ',
     $           'exceeded - either raise maxrec or write a direct ',
     $           'method'
            stop
         endif
         read(iutmp, *, err = 9980, end = 320, iostat = ioval)
     $        (rjunk, i = 1, icol - 1), data(icount)
         ssq = ssq + data(icount) * data(icount)
         rmean = rmean + data(icount)
         go to 300
 320     icount = icount - 1
         close(iutmp)
 330  continue
      ssq = ssq / icount
      rmean = rmean / icount
      write(6,*) rmean
c***********************************************************************
c     move the data set by the mean:
c***********************************************************************
      do 340 i = 1, icount
         data(i) = data(i) - rmean
  340 continue
c***********************************************************************
c     Calculate the time correlation function and write it out:
c***********************************************************************
      write(filnam, *) outfil
      open(iutmp, file = filnam, err = 9990, iostat = ioval)
      write(iutmp, *) 0.0d0, ssq - rmean * rmean
      do 410 ibatch = 1, nbatch
         corrfn = 0.0
         do 400 irec = 1, icount - ibsize(ibatch)
            corrfn = corrfn + data(irec) * data(irec + ibsize(ibatch))
 400     continue
         corrfn = corrfn / (icount - ibsize(ibatch))
c         corrfn = corrfn - rmean * rmean
         write(iutmp, *) ibsize(ibatch) * dt, corrfn
 410  continue
      close(iutmp)
      stop
c***********************************************************************
c     error handling:
c***********************************************************************
 9970 write(iuout, *) 'tcorr: end of file reached in ', filnam
      call ioerr(ioval)
      stop
 9980 write(iuout, *) 'tcorr: error writing to ', filnam
      call ioerr(ioval)
      stop
 9990 write(iuout, *) 'tcorr: error reading from ', filnam
      call ioerr(ioval)
      stop
      end


