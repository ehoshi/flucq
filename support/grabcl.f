      program grabcl
c***********************************************************************
c     grabcl will grab a cluster from a qdyn-generated geometry file.
c     it chooses a random seed atom from the eor file provided as input,
c     and chooses the N-1 closest molecules in order to create an 
c     N-molecule cluster.  The code is sloppy, and makes mistakes such
c     as assuming the run that generated the eor file was periodic.
c     It is also very fragile wrt eor file versions.
c     ...sjs 7/13/95
c***********************************************************************
c
      include 'implic'
      include 'genpar'
c
      parameter (maxmol = 256)
c
      character filnam*80, eorfl*80, datdir*80, outfil*80
      integer*4 iclmol(maxmol)
      logical   incl(maxmol), ranflg
      real*8    pos(maxmol,3),
     $          dist(maxmol),
     $          boxsiz(3)
c
c***********************************************************************
c     Read in the input file grabcl.in
c***********************************************************************
      write(filnam, '(a)') 'grabcl.in'
      open(iutmp, file = filnam, err = 9970, iostat = ioval)
      read(iutmp, *, err = 9970, end = 9960, iostat = ioval) datdir
      read(iutmp, *, err = 9970, end = 9960, iostat = ioval) eorfl
      read(iutmp, *, err = 9970, end = 9960, iostat = ioval) outfil
      read(iutmp, *, err = 9970, end = 9960, iostat = ioval) nclmol
      read(iutmp, *, err = 9970, end = 9960, iostat = ioval) ranflg
      read(iutmp, *, err = 9970, end = 9960, iostat = ioval) iseed
      close(iutmp)
c***********************************************************************
c     scan through the eor file once, keeping the positions of all of
c     the head atoms:
c***********************************************************************
      write(filnam, '(2a)') datdir(1:index(datdir, ' ')-1), 
     $     eorfl(1:index(eorfl, ' ')-1)
      open(iutmp, file = filnam, err = 9990, iostat = ioval)
      read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $     ijunk
      read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $     (boxsiz(i), i = 1, 3)
      read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $     nmol
      do 200 imol = 1, nmol
         read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $        ijunk, ijunk, nfrtmp, ntatmp
         read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $        ijunk, ijunk
         read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $        (pos(imol,i), i = 1, 3)
         read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $        rjunk, rjunk, rjunk
         read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $        rjunk, rjunk
         do 190 ifr = 2, nfrtmp
            read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $           ijunk, ijunk
            read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $           rjunk, rjunk, rjunk
            read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $           rjunk, rjunk, rjunk
            read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $           rjunk, rjunk
  190    continue
         do 191 itg = 1, ntatmp
            read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $           ijunk, ijunk
            read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $           rjunk, rjunk
  191    continue
  200 continue
c***********************************************************************
c     go through the list of molecules, calculating their distances to a 
c     randomly-generated seed molecule, and sorting them to keep only
c     the N-1 closest:
c***********************************************************************
      if (ranflg) then
         iseedm = int(nmol * gran(iseed)) + 1
      else
         iseedm = iseed
      endif
      nincl = 0
      do 300 imol = 1, nmol
         incl(imol) = .false.
         if (imol .eq. iseedm) then
            incl(imol) = .true.
            go to 300
         endif
         dx = pos(imol,1) - pos(iseedm,1)
         dy = pos(imol,2) - pos(iseedm,2)
         dz = pos(imol,3) - pos(iseedm,3)
         dxinc = -1 * dnint(dx / boxsiz(1)) * boxsiz(1)
         dyinc = -1 * dnint(dy / boxsiz(2)) * boxsiz(2)
         dzinc = -1 * dnint(dz / boxsiz(3)) * boxsiz(3)
         dx = dx + dxinc
         dy = dy + dyinc
         dz = dz + dzinc
         r2 = dx ** 2 + dy ** 2 + dz ** 2
         r = sqrt(r2)
         if (nincl .eq. 0) then
            nincl = nincl + 1
            dist(nincl) = r
            iclmol(nincl) = imol
            incl(imol) = .true.
         else if (nincl .lt. nclmol - 1) then
            do 280 iincl = nincl, 1, -1
               if (r .lt. dist(iincl)) then
                  dist(iincl + 1) = dist(iincl)
                  iclmol(iincl + 1) = iclmol(iincl)
               else
                  dist(iincl + 1) = r
                  iclmol(iincl + 1) = imol
                  incl(imol) = .true.
                  go to 281
               endif
  280       continue
            dist(1) = r
            iclmol(1) = imol
            incl(imol) = .true.
  281       continue
            nincl = nincl + 1
         else
            if (r .lt. dist(nincl)) then
               incl(iclmol(nincl)) = .false.
               do 290 iincl = nincl - 1, 1, -1
                  if (r .lt. dist(iincl)) then
                     dist(iincl + 1) = dist(iincl)
                     iclmol(iincl + 1) = iclmol(iincl)
                  else
                     dist(iincl + 1) = r
                     iclmol(iincl + 1) = imol
                     incl(imol) = .true.
                     go to 291
                  endif
  290          continue
               dist(1) = r
               iclmol(1) = imol
               incl(imol) = .true.
  291          continue
            endif
         endif
  300 continue
c***********************************************************************
c     now read the (still open) eor file once more, spitting back out 
c     only the molecules that made the cut, in sys file form:
c***********************************************************************
      icount = 0
      dxs = -pos(iseedm,1)
      dys = -pos(iseedm,2)
      dzs = -pos(iseedm,3)
      rewind(iutmp)
      write(filnam, '(2a)') datdir(1:index(datdir, ' ')-1), 
     $     outfil(1:index(outfil, ' ')-1)
      open(iutmp2, file = filnam, err = 9990, iostat = ioval)
      read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $     ijunk
      write(iutmp2, *, err = 9970, iostat = ioval)
     $     ijunk
      read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $     rjunk, rjunk, rjunk
      write(iutmp2, *, err = 9970, iostat = ioval)
     $     (boxsiz(i), i = 1, 3)
      read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $     ijunk
      write(iutmp2, *, err = 9970, iostat = ioval)
     $     nclmol
      do 400 imol = 1, nmol
         read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $        ijunk1, ijunk2, nfrtmp, ntatmp
         if (incl(imol)) then
            icount = icount + 1
            write(iutmp2, *, err = 9970, iostat = ioval)
     $           icount, ijunk2, nfrtmp, ntatmp
         endif
         read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $        ijunk1, ijunk2
         if (incl(imol)) then
            write(iutmp2, *, err = 9970, iostat = ioval)
     $           ijunk1, ijunk2
         endif
         read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $        x, y, z
         if (incl(imol)) then
            dx = x - pos(iseedm,1)
            dy = y - pos(iseedm,2)
            dz = z - pos(iseedm,3)
            dxinc = -1 * dnint(dx / boxsiz(1)) * boxsiz(1)
            dyinc = -1 * dnint(dy / boxsiz(2)) * boxsiz(2)
            dzinc = -1 * dnint(dz / boxsiz(3)) * boxsiz(3)
            write(iutmp2, *, err = 9970, iostat = ioval)
     $           x + dxinc + dxs, y + dyinc + dys, z + dzinc + dzs
         endif
         read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $        rjunk, rjunk, rjunk
         read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $        rjunk, rjunk
         do 390 ifr = 2, nfrtmp
            read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $           ijunk1, ijunk2
            if (incl(imol)) then
               write(iutmp2, *, err = 9970, iostat = ioval)
     $              ijunk1, ijunk2
            endif
            read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $           x, y, z
            if (incl(imol)) then
               write(iutmp2, *, err = 9970, iostat = ioval)
     $              x + dxinc + dxs, y + dyinc + dys, z + dzinc + dzs
            endif
            read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $           rjunk, rjunk, rjunk
            read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $           rjunk, rjunk
  390    continue
         do 391 itg = 1, ntatmp
            read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $           ijunk1, ijunk2
            if (incl(imol)) then
               write(iutmp2, *, err = 9970, iostat = ioval)
     $              ijunk1, ijunk2
            endif
            read(iutmp, *, end = 9960, err = 9970, iostat = ioval)
     $           rjunk, rjunk
  391    continue
  400 continue
      close(iutmp)
      close(iutmp2)
      write(iuout, '(i3,a,i3,a,a,a,a)') icount, 
     $     '-molecule cluster from around molecule #',
     $     iseedm, ' of ', 
     $     eorfl(1:index(eorfl, ' ')-1), 
     $     ' written to ', 
     $     outfil(1:index(outfil, ' ')-1)
      go to 9999
 9960 write(iuout, *) 'grabcl: end of file in ', filnam
      call ioerr(ioval)
      stop
 9970 write(iuout, *) 'grabcl: error reading from ', filnam
      call ioerr(ioval)
      stop
 9980 write(iuout, *) 'grabcl: error writing to ', filnam
      call ioerr(ioval)
      stop
 9990 write(iuout, *) 'grabcl: error opening ', filnam
      call ioerr(ioval)
      stop
 9999 continue
      end
