      program hdpe
c***********************************************************************
c     hdpe calculates the heterodimer PE surface for a Cl-/H2O 
c     heterodimer.  There are only 2 degrees of freedom:  the separation
c     and the angle of the H bond....sjs 5/4/95
c***********************************************************************
c     added capability to do a water dimer.  Here the two degrees
c     of freedom are the separation and the angle of the donor from the
c     plane.  The acceptor is assumed to be linear....sjs 6/9/95
c***********************************************************************
c
      include 'implic'      
      include 'genpar'
      include 'qpar'
      include 'commons'
c     
      character filenm*80
      integer*4 nfrmid(maxmid), ntamid(maxmid),
     $          icjid(mxcstj,3)
      logical   knwmid(maxmid), 
     $          getnr, getfar, getlt, gethvy, gete, anaflg
c
      real*8 J
c
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
      neede = .true.
      do 110 imid = 1, maxmid
         knwmid(imid) = .false.
 110  continue
c***********************************************************************
c     Read in a header file, even though most info isn't needed.  It's
c     just convenient to use qdyn's routines to build the cross
c     reference arrays, etc.  Note that this file has to be for a 
c     Cl-/H2O heterodimer system (specified as the syfile in hdpe.dat)
c     with the Cl- as the first molecule.
c***********************************************************************
      write(filnam, '(a)') 'hdpe.in'
      open(iutmp, file = filnam, err = 9980, iostat = ioval)
      read(iutmp, *) datdir
      read(iutmp, *) nrstep, rlow, rhigh
      read(iutmp, *) nastep, alow, ahigh
      alow = alow * pi / 180.d0
      ahigh = ahigh * pi / 180.d0
      write(filenm, *) datdir(1:index(datdir, ' ')-1), 'hdpe.dat'
      write(iuout, *) 'reading header info from hdpe.dat'
      call rddet(filenm, icjid)
      call chkdet(icjid)
      call setdet
      call rdmodb(knwmid, nfrmid, ntamid)
      write(filnam, *) datdir(1:index(datdir, ' ')-1), syfile
      anaflg = .false.
      call rdeor(filnam, knwmid, nfrmid, ntamid, anaflg)
      call setsys(anaflg)
      call chksys(icjid)
      call rdmodd
      call setmod(icjid)
      call tagup(imixed)
      call conchk
      call trmjst
      do 205 imol = 1, nmol
         imty = molty(imol)
         do 200 iqind = 1, nrqmol(imty)
            imind = irqmol(imty,iqind)
            iatom = iatmol(imol,imind)
            iaty = iatype(imty,imind)
            q(iatom) = defq(iaty)
  200    continue
  205 continue
      do 210 imty = 1, nmty
         call mpeset(imty)
  210 continue
      basepe = 0.
      do 220 imol = 1, nmol
         basepe = basepe + rmonpe(molty(imol))
  220 continue
      basepe = 2. * basepe
      imty = molty(2)
      roh = sqrt(r2cons(imty,1))
      rhh = sqrt(r2cons(imty,3))
      hohb2 = asin((rhh/2.d0)/roh)
      moving = .false.
      getnr = .true.
      getfar = .true.
      getlt = .true.
      gethvy = .true.
      gete = .true.
      pliter = .true.
c$$$      dofld = .true.
c$$$      write(filnam, *) datdir(1:index(datdir, ' ')-1), 'testJ.', suffix
c$$$      open(iutmp, file = filnam, err = 9990, iostat = ioval)
c$$$      do 230 i = 0, 97
c$$$         r = 0.d0 + i * dlimr / 97.d0
c$$$         test = J(molty(1), molty(2), 1, 2, r, Jreg)
c$$$         diff = test - quadj(iatype(molty(1),1), iatype(molty(2),2), r)
c$$$         test2 = dJdr(molty(1), molty(2), 1, 2, r, Jreg)
c$$$         diff2 = test2 - 
c$$$     $        quaddj(iatype(molty(1),1), iatype(molty(2),2), r)
c$$$         write(iutmp, *) r, test, diff, test2, diff2
c$$$  230 continue
c$$$      close(iutmp)
c$$$      stop
      write(filnam, *) datdir(1:index(datdir, ' ')-1), 'hdpe.', suffix
      open(10, file = filnam, err = 9990, iostat = ioval)
      write(filnam, *) datdir(1:index(datdir, ' ')-1), 'hdbe.', suffix
      open(11, file = filnam, err = 9990, iostat = ioval)
      if (maxmid .ne. 11) then
         write(iuout, *) 'hdpe: program not up to date'
         stop
      endif
      if (hasmid(5) .or. hasmid(6) .or. hasmid(9) .or. hasmid(10) .or.
     $     hasmid(11)) then
         write(iuout, *) 'epsilon: ', trmlj(iatype(2,1),iatype(2,1),1)
         write(iuout, *) 'sigma: ', trmlj(iatype(2,1),iatype(2,1),4)
      endif
c***********************************************************************
c     loop over all geometries, calculating energy for each:
c***********************************************************************
      alowe = bigpos
      do 340 irocl = 0, nrstep
         rocl = rlow + irocl * (rhigh - rlow) / nrstep
         rlowe = bigpos
         do 330 iangle = 0, nastep
            angle = alow + iangle * (ahigh - alow) / nastep
c***********************************************************************
c     This sets the geometry for the Cl-/H2O dimer:
c***********************************************************************
            if ((hasmid(5) .or. hasmid(6) .or. hasmid(9) .or.
     $           hasmid(10) .or. hasmid(11)) .and. 
     $           (hasmid(1) .or. 
     $           hasmid(2) .or. hasmid(3) .or. hasmid(4) .or.
     $           hasmid(7) .or. hasmid(8))) then
               iatom = iatmol(1,1)
               pos(iatom,1) = 0.d0
               pos(iatom,2) = 0.d0
               pos(iatom,3) = 0.d0
               if (dslvfl) then
                  iatom = iatmol(1,2)
                  pos(iatom,1) = 0.d0
                  pos(iatom,2) = 0.d0
                  pos(iatom,3) = 0.d0
               endif
               iatom = iatmol(2,1)
               pos(iatom,1) = rocl
               pos(iatom,2) = 0.d0
               pos(iatom,3) = 0.d0
               iatom = iatmol(2,2)
               pos(iatom,1) = rocl - roh * cos(hohb2 + angle)
               pos(iatom,2) = roh * sin(hohb2 + angle)
               pos(iatom,3) = 0.d0
               iatom = iatmol(2,3)
               pos(iatom,1) = rocl - roh * cos(hohb2 - angle)
               pos(iatom,2) = -roh * sin(hohb2 - angle)
               pos(iatom,3) = 0.d0
c***********************************************************************
c     This sets the geometry for the water dimer:
c***********************************************************************
            else if (.not. hasmid(5) .and. .not. hasmid(6) .and. 
     $              .not. hasmid(9) .and. .not. hasmid(10) .and.
     $              .not. hasmid(11)) then
               iatom = iatmol(1,1)
               pos(iatom,1) = 0.d0
               pos(iatom,2) = 0.d0
               pos(iatom,3) = 0.d0
               iatom = iatmol(1,2)
               pos(iatom,1) = -roh * cos(hohb2) * cos(angle)
               pos(iatom,2) = -roh * sin(hohb2)
               pos(iatom,3) = roh * cos(hohb2) * sin(angle)
               iatom = iatmol(1,3)
               pos(iatom,1) = -roh * cos(hohb2) * cos(angle)
               pos(iatom,2) = roh * sin(hohb2)
               pos(iatom,3) = roh * cos(hohb2) * sin(angle)
               iatom = iatmol(2,1)
               pos(iatom,1) = rocl
               pos(iatom,2) = 0.d0
               pos(iatom,3) = 0.d0
               iatom = iatmol(2,2)
               pos(iatom,1) = rocl - roh
               pos(iatom,2) = 0.d0
               pos(iatom,3) = 0.d0
               iatom = iatmol(2,3)
               pos(iatom,1) = rocl - roh * cos(hohb2 * 2.d0)
               pos(iatom,2) = 0.d0
               pos(iatom,3) = -roh * sin(hohb2 * 2.d0)
            endif
            call tagup(imixed)
            if (qslvfl) then
               qwrong = .true.
            else
               qwrong = .false.
            endif
            if (dslvfl) then
               dwrong = .true.
            else
               dwrong = .false.
            endif
  300       continue
            if (qwrong .and. dwrong) then
               if (qsiter) then
                  if (dslvfl) then
                     dsiter = .true.
                  endif
                  qsiter = .false.
               else
                  dsiter = .false.
                  qsiter = .true.
               endif
            else if (dwrong) then
               dsiter = .true.
               qsiter = .false.
            else if (qwrong) then
               dsiter = .false.
               qsiter = .true.
            else
               dsiter = .false.
               qsiter = .false.
            endif
c$$$            if (dsiter) then
c$$$               dofld = .true.
c$$$            endif
            movlt = .true.
            movhvy = .true.
            call getf(getnr, getfar, getlt, gethvy, gete)
            if (qsiter) then
               call qsolv2
               rms = 0.d0
               do 320 imol = 1, nmol
                  imty = molty(imol)
                  do 310 iqind = 1, nfqmol(imty)
                     iatom = iatmol(imol,ifqmol(imty,iqind))
                     dq = q(iatom) - qans(iatom)
                     rms = rms + dq * dq
                     q(iatom) = qans(iatom)
  310             continue
  320          continue
               rms = sqrt(rms / nfqatm)
               if (rms .gt. qrmsct .and. dslvfl) then
                  qwrong = .true.
                  dwrong = .true.
               else
                  qwrong = .false.
               endif
            endif
            if (dsiter) then
               call dsolv
            endif
            if (qsiter .or. dsiter) then
               go to 300
            endif
c$$$            apetot = qpe(ilight) + qpe(imixed) + qpe(iheavy) +
c$$$     $           ape(ilight) + ape(imixed) + ape(iheavy)
            apetot = ape(ilight) + ape(imixed) + ape(iheavy)
c***********************************************************************
c     this factor of two is to turn the energy into a per-molecule
c     energy:
c***********************************************************************
            apetot = apetot * 0.5d0
            if (dslvfl) then
               dr = sqrt((pos(1,1) - pos(2,1)) ** 2 +
     $              (pos(1,2) - pos(2,2)) ** 2 + 
     $              (pos(1,3) - pos(2,3)) ** 2)
            endif
            if (apetot .lt. rlowe) then
               rlowe = apetot
               rlowan = angle
               rlowdr = dr
            endif
            write(10, *) rocl, angle * 180.d0 / pi, apetot * fjm2kc
  330    continue
         if (rlowe .lt. alowe) then
            alowe = rlowe
            alowan = rlowan
            alowr = rocl
         endif
         write(10, '(\)')
         write(11, *) rocl, rlowan * 180.d0 / pi, rlowe * fjm2kc, rlowdr
  340 continue
      write(iuout, '(a,f9.3,a,f6.3,a,f6.1,a)') 'lowest energy: ', 
     $     alowe * fjm2kc, ' kcal/mol at r = ', 
     $     alowr, ' A and an angle of ', alowan * 180.d0 / pi, 
     $     ' degrees'
      close(10)
      close(11)
      go to 9999
 9980 write(iuout, *) 'props: error opening file ',filnam
      call ioerr(ioval)
      stop
 9990 write(iuout, *) 'rdsim: error opening ', filnam
      rdsim = ioval
      return
 9999 continue
      end
