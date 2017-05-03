c***********************************************************************
c     rdsim reads in the simulation details from the specified input
c     file.  The return value is zero, unless an error is encountered.
c     A return value of 1 means the end of file was hit, -1 means there
c     was an error.  If called with ifull = 1, only the first part of
c     the file will be read.  ifull = 2 reads more....sjs 2/24/93
c***********************************************************************
      integer*4 function rdsim()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
c
      integer*4 nfrmid(maxmid), ntamid(maxmid),
     $          icjid(mxcstj,3)
      logical   anaflg,
     $          knwmid(maxmid)
      real*8    cmvel(3)
c
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
c$$$      symass = 0.
      do 100 imid = 1, maxmid
         knwmid(imid) = .false.
 100  continue
c***********************************************************************
c     Read in the simulation details from STDIN and initialize 
c     associated stuff:
c***********************************************************************
      write(filnam, '(a)') 'STDIN'
      call rddet(iuin, icjid)
      call chkdet(icjid)
      call setdet
c***********************************************************************
c     Read in some basic details about all of the known molecule ids,
c     from models.dat:
c***********************************************************************
      call rdmodb(knwmid, nfrmid, ntamid)
c***********************************************************************
c     Read in the startup or .sys file and initialize associated stuff:
c***********************************************************************
      if (stfile .eq. 'none') then
         if (syfile .eq. 'none') then
            write(filnam, '(2a)') datdir(1:index(datdir, ' ')-1), 
     $           '/system.dat'
         else
            write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), '/',
     $           syfile
         endif
      else
         if (syfile .ne. 'none') write(iuout, *)
     $        'rdsim: WARNING! using ', stfile, ' ignoring ', syfile
         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), '/', 
     .        stfile
      endif
      anaflg = .false.
      call rdeor(filnam, knwmid, nfrmid, ntamid, anaflg)
c***********************************************************************
c     set all the variables related to the stuff just read in, and check
c     for illegal values and incompatiblities:
c***********************************************************************
      call setsys(anaflg)
      call chksys(icjid)
c***********************************************************************
c     set up the molecule id <-> molecule type cross-references:
c***********************************************************************
      call rdmodd
      call setmod(icjid)
c$$$      call tagup(trumsk)
      call tagup(imixed)
      call gtcomv(cmvel)
      if (abs(cmvel(1)) .gt. 1.d-6 .or. abs(cmvel(2)) .gt. 1.d-6 .or.
     $     abs(cmvel(3)) .gt. 1.d-6) then
         write(iuout, *) 'rdsim: WARNING!  the c.o.m. was drifting at',
     $        sqrt(cmvel(1) ** 2 + cmvel(2) ** 2 + 
     $        cmvel(3) ** 2) * fA2cm * fs2smd, ' cm/s'
         do 190 imol = 1, nmol
            imty = molty(imol)
            do 180 imind = 1, nframe(imty)
               iatom = iatmol(imol,imind)
               vel(iatom,1) = vel(iatom,1) - cmvel(1)
               vel(iatom,2) = vel(iatom,2) - cmvel(2)
               vel(iatom,3) = vel(iatom,3) - cmvel(3)
  180       continue
  190    continue
      endif
c$$$      if (pbflag) call foldr(trumsk)
      if (pbflag) then
         call foldr(imixed)
      endif
      call conchk
      call trmjst
      do 410 imty = 1, nmty
         call mpeset(imty)
 410  continue
      basepe = 0.
      do 411 imol = 1, nmol
         basepe = basepe + rmonpe(molty(imol))
 411  continue
      basepe = 2. * basepe
c***********************************************************************
c     If we're doing charge dynamics, take charges from qifile if
c     possible, or stfile next, or choose them randomly if necessary.
c     Charge velocities and accelerations are zero unless taken from
c     stfile.:
c***********************************************************************
      if (qifile .ne. 'none') then
         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), '/', 
     .        qifile
         open(iutmp, file = filnam, err = 9990, iostat = ioval)
         do 600 iatom = 1, natoms
            read(iutmp, *) jatom, ityp, imol, temp
            if (jatom .ne. iatom) then
               write(isterr, *) 'atom ', iatom, ' != ', jatom
               stop
            endif
            q(iatom) = temp
            if (molec(iatom) .ne. imol .or. 
     $           ident(iatom) .ne. ityp) then
               write(iuout, *) 'rdsim: mismatch in ', qifile
               stop
            endif
            qvel(iatom) = 0
c$$$            qacc(iatom) = 0
 600     continue
         close(iutmp)
         if (stfile .ne. 'none') write(iuout, *) 'rdsim: WARNING! ',
     $        'using charges from ', qifile, ', not from ', stfile
      else if (stfile .eq. 'none') then
         do 615 imol = 1, nmol
            imty = molty(imol)
            do 610 iqind = 1, nfqmol(imty)
               iatom = iatmol(imol,ifqmol(imty,iqind))
               q(iatom) = qlowg(ident(iatom)) + gran(iseed) * 
     $              (qhighg(ident(iatom)) - qlowg(ident(iatom)))
               qvel(iatom) = 0
c$$$               qacc(iatom) = 0
 610        continue
            do 612 iqind = 1, nrqmol(imty)
               imind = irqmol(imty,iqind)
               iatom = iatmol(imol,imind)
               iaty = iatype(imty,imind)
               q(iatom) = defq(iaty)
 612        continue
 615     continue
         if (ctflag) then
            temp = 0.0d0
            do 625 imol = 1, nmol
               imty = molty(imol)
               do 620 iqind = 1, nqmol(imty)
                  iatom = iatmol(imol,iqmol(imty,iqind))
                  temp = temp + q(iatom)
 620           continue
 625        continue
            temp = temp - qsys
            temp = temp / nfqatm
            do 635 imol = 1, nmol
               imty = molty(imol)
               do 630 iqind = 1, nfqmol(imty)
                  iatom = iatmol(imol,ifqmol(imty,iqind))
                  q(iatom) = q(iatom) - temp
 630           continue
 635        continue
         else
            do 660 imol = 1, nmol
               temp = 0.0d0
               imty = molty(imol)
               do 640 iqind = 1, nqmol(imty)
                  imind = iqmol(imty,iqind)
                  temp = temp + q(iatmol(imol,imind))
 640           continue
               temp = temp - qmol(imty)
               if (hasfqt(imty)) then
                  temp = temp / nfqmol(imty)
                  do 650 iqind = 1, nfqmol(imty)
                     imind = ifqmol(imty,iqind)
                     q(iatmol(imol,imind)) = 
     $                    q(iatmol(imol,imind)) - temp
 650              continue
               else
                  if (temp .ne. 0.0) then
                     write(iuout, '(a,i4,a)') 
     $                    'rdsim:  fixed-q molecule ',
     $                    imol, ' has the wrong charge'
                     stop
                  endif
               endif
 660        continue
         endif
      else
         do 670 imol = 1, nmol
            imty = molty(imol)
            do 665 iqind = 1, nrqmol(imty)
               imind = irqmol(imty,iqind)
               iatom = iatmol(imol,imind)
               iaty = iatype(imty,imind)
               if (q(iatom) .ne. defq(iaty)) then
                  write(6, '(a,i4,a)') 
     $                 'rdsim:  WARNING! atom ', iatom, 
     $                 ' charge incorrect; charge has been fixed'
               endif
               q(iatom) = defq(iaty)
  665       continue
  670    continue
      endif
      rdsim = 0
      return
 9970 write(iuout, *) 'rdsim: error reading from ', filnam
      rdsim = ioval
      return
 9980 write(iuout, *) 'rdsim: end of file reached in ', filnam
      rdsim = ioval
      return
 9990 write(iuout, *) 'rdsim: error opening ', filnam
      rdsim = ioval
      return
      end
c
c***********************************************************************
c     wrtsim writes out the simulation details to the file filnam.
c     ...sjs 2/24/93
c***********************************************************************
      subroutine wrtsim()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
     $     '/simdat.', suffix
      write(iuout, *) 'writing simulation details to simdat.', suffix
      open(iutmp, file = filnam, err = 9990, iostat = ioval)
      write(iutmp, *, err = 9980, iostat = ioval) iosimv
      write(iutmp, '(5a)', err = 9980, iostat = ioval) '''', 
     $     datdir(1:index(datdir, ' ')-1), ''' ''',
     $     scrdir(1:index(scrdir, ' ')-1), ''''
      write(iutmp, '(7a)', err = 9980, iostat = ioval) '''',
     .     stfile(1:index(stfile, ' ')-1), ''' ''',
     .     syfile(1:index(syfile, ' ')-1), ''' ''',
     .     qifile(1:index(qifile, ' ')-1), ''''
      write(iutmp, '(3a)', err = 9980, iostat = ioval) '''', rinfo, ''''
      write(iutmp, *, err = 9980, iostat = ioval)
     $     ctflag, qslvfl, nqsolv, qmass
      write(iutmp, *, err = 9980, iostat = ioval) 
     $     ndstep, pbflag, dtbig, nl, nlf, nhf
      write(iutmp, *, err = 9980, iostat = ioval)
     $     dslvfl, ndsolv
      write(iutmp, '(5f11.5,i8,f11.5)', err = 9980, iostat = ioval) 
     $     fcutr, fcutd, fcutb, fnearr, fneard, iplint, hvywt
      write(iutmp, *, err = 9980, iostat = ioval) 
     $     interJ, intraJ, ncint
      nccint = 0
      if (hscstj) then
         do 110 imty1 = 1, nmty
            do 100 imty2 = 1, imty1
               if (icustj(imty1,imty2) .ne. interJ) then
                  write(iutmp, *, err = 9980, iostat = ioval)
     $                 midmty(imty1), midmty(imty2), icustj(imty1,imty2)
                  nccint = nccint + 1
               endif
  100       continue
  110    continue
         if (nccint .ne. ncint) then
            write(iuout, *) 'wrtsim:  found ', nccint, ' custom ',
     $           'J(r) interactions, not ', ncint
            stop
         endif
      endif
      write(iutmp, *, err = 9980, iostat = ioval)
     $     ewlflg, ewlkfl, ewsrfl, ewlkap, ewlkmx
      write(iutmp, 
     .    '(f10.5,1x,f10.5,1x,i5,1x,f11.5,1x,l,i8,i8,f10.5,f11.5,1x,l)',
     .     err = 9980, iostat = ioval) T, Ttol, nscint,stint, vscflg,
     .     nsceql, mxsceq, qT, qTrang, qvscfl
      write(iutmp, *, err = 9980, iostat = ioval) cbnflg, 
     $     dbflag, iseed, nioint, hbonde * fJm2kc, 
cmov     $     gstinf,
     $     ocflag, (Efield(i), i = 1, 3)
      write(iutmp, '(l,1x,l,1x,i8,1x,3a)', err = 9980, iostat = ioval) 
     .     ioflg, binflg, nbnint,
cmov     $     iobwfl, 
     $     '''', suffix, ''''
      close(iutmp)
      return
 9980 write(iuout, *) 'wrtsim: error writing to file ', filnam
      call ioerr(ioval)
      stop
 9990 write(iuout, *) 'wrtsim: error opening file ', filnam
      call ioerr(ioval)
      stop
      end
c
cmovc***********************************************************************
cmovc     sumbin scales and flushes all the binning stats....sjs 3/18/93
cmovc     it also now accumulates some averaging stats....sjs 8/18/93
cmovc***********************************************************************
cmovc     minor correction: some of the scaling stats should only be done 
cmovc     once, not every time this routine is called.
cmovc***********************************************************************
cmov      subroutine sumbin()
cmovc
cmov      include 'implic'
cmov      include 'qpar'
cmov      include 'genpar'
cmov      include 'commons'
cmovc     
cmov      character name*12
cmov      real*4 rio(maxbin), riovar(maxbin)
cmov      real*8 gscal2(maxaty,maxaty)
cmovc
cmovc***********************************************************************
cmovc     Open the binary binning data file and write out the header info:
cmovc***********************************************************************
cmov      write(filnam, *) datdir, '/bins.', suffix
cmov      open(iutmp, file = filnam, form = 'unformatted', err = 9990, 
cmov     $     iostat = ioval)
cmov      write(iutmp) iobinv
cmov      itemp = nrlaty * (nrlaty + 1) / 2
cmov      if (iobwfl) then
cmov         itemp = itemp + 5
cmov      endif
cmov      if (.not. pbflag) then
cmov         itemp = itemp + nmol
cmov      endif
cmovc$$$      itemp = 5 + nmaaty * (nmaaty + 1) / 2
cmovc$$$      itemp = 5 + naty * (naty + 1) / 2 +
cmovc$$$     $     2 * naty * (naty - 1) / 2 + naty
cmovc$$$      if (watflg .and. lflex) itemp = itemp + 3
cmovc$$$      if (watflg .and. hasbdt(iwatty)) itemp = itemp + 3
cmov      if (qdflag) itemp = itemp + nfqaty
cmov      write(iutmp, err = 9980, iostat = ioval) itemp
cmovc***********************************************************************
cmovc     Write out the pair energy stats:
cmovc***********************************************************************
cmov      if (iobwfl) then
cmov         write(name, *) 'pairE'
cmov         write(iutmp, err = 9980, iostat = ioval) name, maxbin, 
cmov     $        prEbmn * fJm2kc, prEbmx * fJm2kc, .false., .false.
cmov         temp = 2 * prEbwi / (nbinst * nmol * (nmol - 1) * fJm2kc)
cmov         do 108 ibin = 1, maxbin
cmov            rio(ibin) = npaire(ibin) * temp
cmov 108     continue
cmov         if (irdump(iutmp, rio, maxbin) .ne. 0) go to 9980
cmov      endif
cmovc***********************************************************************
cmovc     Write out the bonding energy stats:
cmovc***********************************************************************
cmov      if (iobwfl) then
cmov         write(name, *) 'bondE'
cmov         write(iutmp, err = 9980, iostat = ioval) name, maxbin, 
cmov     $        bdEbmn * fJm2kc, bdEbmx * fJm2kc, .false., .false.
cmov         temp = bdEbwi / (nbinst * nmol * fJm2kc)
cmov         do 110 ibin = 1, maxbin
cmov            rio(ibin) = nbonde(ibin) * temp
cmov 110     continue
cmov         if (irdump(iutmp, rio, maxbin) .ne. 0) go to 9980
cmov      endif
cmovc***********************************************************************
cmovc     Write out the H-bond count stats:
cmovc***********************************************************************
cmov      if (iobwfl) then
cmov         write(name, *) 'Hbond'
cmov         write(iutmp, err = 9980, iostat = ioval) name, mxhbp1, 
cmov     $        -0.5d0, mxhbp1 - 0.5d0, .false., .false.
cmov         temp = 1.0d0 / dble(nmol * nbinst)
cmov         do 111 ihbp1 = 1, mxhbp1
cmov            rio(ihbp1) = nwhb(ihbp1) * temp
cmov 111     continue
cmov         if (irdump(iutmp, rio, mxhbp1) .ne. 0) go to 9980
cmov      endif
cmovc***********************************************************************
cmovc     Write out the G_K(r) stats:
cmovc***********************************************************************
cmov      if (iobwfl) then
cmov         write(name, *) 'gkr'
cmov         write(iutmp, err = 9980, iostat = ioval) name, maxbin, 
cmov     $        rsqbmn, rsqbmx, .true., .true.
cmov         temp = feA2D * feA2D / (nbinst * nmol)
cmov         do 112 ibin = 1, maxbin
cmov            rio(ibin) = gkr(ibin) * temp
cmov 112     continue
cmov         if (irdump(iutmp, rio, maxbin) .ne. 0) go to 9980
cmov      endif
cmovc***********************************************************************
cmovc     Write out the angle pair correlation stats:
cmovc***********************************************************************
cmov      if (iobwfl) then
cmov         write(name, *) 'angcor'
cmov         write(iutmp, err = 9980, iostat = ioval) 
cmov     $        name, maxbin, rsqbmn, rsqbmx, .false., .true.
cmov         do 113 ibin = 1, maxbin
cmov            if (nang(ibin) .gt. 0) then
cmov               rio(ibin) = ang(ibin) / nang(ibin)
cmov            else
cmov               rio(ibin) = 0.0d0
cmov            endif
cmov 113     continue
cmov         if (irdump(iutmp, rio, maxbin) .ne. 0) go to 9980
cmov      endif
cmovc$$$      if (watflg .and. hasbdt(iwatty)) then
cmovc$$$c***********************************************************************
cmovc$$$c     Write out the water O-H bond length distribution stats:
cmovc$$$c***********************************************************************
cmovc$$$         write(name, *) 'rOH'
cmovc$$$         write(iutmp, err = 9980, iostat = ioval) name, maxbin, 
cmovc$$$     $        rOHbmn, rOHbmx, .false., .false.
cmovc$$$         temp = rOHbwi / (2 * nbinst * nmol)
cmovc$$$         do 120 ibin = 1, maxbin
cmovc$$$            rio(ibin) = nrOH(ibin) * temp
cmovc$$$  120    continue
cmovc$$$         if (rdump(iutmp, rio, maxbin) .ne. 0) go to 9980
cmovc$$$c***********************************************************************
cmovc$$$c     Write out the water H-H bond length distribution stats:
cmovc$$$c***********************************************************************
cmovc$$$         write(name, *) 'rHH'
cmovc$$$         write(iutmp, err = 9980, iostat = ioval) name, maxbin,
cmovc$$$     $        rHHbmn, rHHbmx, .false., .false.
cmovc$$$         temp = rHHbwi / (nbinst * nmol)
cmovc$$$         do 130 ibin = 1, maxbin
cmovc$$$            rio(ibin) = nrHH(ibin) * temp
cmovc$$$  130    continue
cmovc$$$         if (rdump(iutmp, rio, maxbin) .ne. 0) go to 9980
cmovc$$$c***********************************************************************
cmovc$$$c     Write out the water H-O-H bond angle distribution stats:
cmovc$$$c***********************************************************************
cmovc$$$         write(name, *) 'HOH'
cmovc$$$         write(iutmp, err = 9980, iostat = ioval) name, maxbin,
cmovc$$$     $        HOHbmn, HOHbmx, .false., .false.
cmovc$$$         temp = HOHbwi / (nbinst * nmol)
cmovc$$$         do 140 ibin = 1, maxbin
cmovc$$$            rio(ibin) = nHOH(ibin) * temp
cmovc$$$  140    continue
cmovc$$$         if (rdump(iutmp, rio, maxbin) .ne. 0) go to 9980
cmovc$$$      endif
cmov      if (qdflag) then
cmovc***********************************************************************
cmovc     Write out the charge distribution stats:
cmovc***********************************************************************
cmov         do 150 iaty = 1, naty
cmov            if (isfqat(iaty)) then
cmov               iaid = ideaty(iaty)
cmov               write(name, '(a,i2.2)') 'q', iaid
cmov               write(iutmp, err = 9980, iostat = ioval) name, maxbin,
cmov     $              qlowg(iaid), qhighg(iaid), .false., .false.
cmov               temp = qbwi(iaty) / (nbinst * nwaty(iaty))
cmov               do 145 ibin = 1, maxbin
cmov                  rio(ibin) = nq(ibin,iaty) * temp
cmov 145           continue
cmov               if (irdump(iutmp, rio, maxbin) .ne. 0) go to 9980
cmov            endif
cmov 150     continue
cmov      endif
cmovc***********************************************************************
cmovc     Write out all of the pair correlation function stats, and the
cmovc     pair charge correlation function stats:
cmovc***********************************************************************
cmovc     minor correction: these calcs can be done faster
cmovc***********************************************************************
cmovc     correction: fix for heterogeneous molecules
cmovc***********************************************************************
cmov      rideal = 4.0d0 * pi * (rndens / 3.0d0) / 3.0d0
cmov      do 165 iaty1 = 1, naty
cmov         if (isrl(iaty1)) then
cmov            imty1 = imtaty(iaty1)
cmov            iaid1 = ideaty(iaty1)
cmov            n1 = 0
cmov            do 161 imind1 = 1, natmty(imty1)
cmov               if (iatype(imty1,imind1) .eq. iaty1) n1 = n1 + 1
cmov 161        continue
cmov            gscal2(iaty1,iaty1) = 1. / (n1 * n1)
cmov            do 164 iaty2 = 1, iaty1 - 1
cmov               if (isrl(iaty2)) then
cmov                  imty2 = imtaty(iaty2)
cmov                  iaid2 = ideaty(iaty2)
cmov                  n2 = 0
cmov                  do 162 imind2 = 1, natmty(imty2)
cmov                     if (iatype(imty2,imind2) .eq. iaty2) n2 = n2 + 1
cmov 162              continue
cmov                  gscal2(iaty1,iaty2) = 0.5 / (n1 * n2)
cmov               endif
cmov 164        continue
cmov         endif
cmov 165  continue
cmov      do 190 ityp1 = 1, naty
cmov         if (isrl(ityp1)) then
cmov            do 180 ityp2 = 1, ityp1
cmov               if (isrl(ityp2)) then
cmov                  if (ideaty(ityp1) .ge. ideaty(ityp2)) then
cmov                     iaid1 = ideaty(ityp1)
cmov                     iaid2 = ideaty(ityp2)
cmov                  else
cmov                     iaid1 = ideaty(ityp2)
cmov                     iaid2 = ideaty(ityp1)
cmov                  endif
cmov                  write(name, '(''g'', 2i2.2)') iaid1, iaid2
cmov                  write(iutmp, err = 9980, iostat = ioval) 
cmov     $                 name, maxbin, rbmn, rbmx, .false., .false.
cmov                  do 170 ibin = 1, maxbin
cmov                     gr = (ibin - 0.5d0) * rbw
cmov                     if (cubflg .and. pbflag) then
cmov                        if (gr .le. hbox) then
cmov                           gscale = 1.0d0
cmov                        else if (gr .le. hbox * dsqrt(2.0d0)) then
cmov                           gscale = 1.0d0 / (3 * hbox / gr - 2)
cmov                        else
cmov                           gscale = 1.0d0
cmov                        endif
cmov                     else
cmov                        gscale = 1.0d0
cmov                     endif
cmov                     temp = gscale / (rideal * nmol *
cmov     $                    ((ibin * rbw) ** 3 - ((ibin - 1) * rbw) ** 3))
cmov                     rio(ibin) = rg(ibin,ityp1,ityp2) * temp / nbinst *
cmov     $                    gscal2(ityp1,ityp2)
cmov                     riovar(ibin) = rg2(ibin,ityp1,ityp2) * temp *
cmov     $                    temp / nbinst * gscal2(ityp1,ityp2) *
cmov     $                    gscal2(ityp1,ityp2) - rio(ibin) * rio(ibin)
cmov 170              continue
cmov                  if (irdump(iutmp, rio, maxbin) .ne. 0) go to 9980
cmov                  if (irdump(iutmp, riovar, maxbin) .ne. 0) go to 9980
cmovc$$$            if (ityp1 .eq. ityp2) then
cmovc$$$               write(name, '(''gq'', 2i2.2, ''-'', i2.2)')
cmovc$$$     $              iaid1, iaid2, iaid1
cmovc$$$               write(iutmp, err = 9980, iostat = ioval) name, maxbin,
cmovc$$$     $              rbmn, rbmx, .false., .false.
cmovc$$$               do 173 ibin = 1, maxbin
cmovc$$$c$$$                  if (ng(ibin,ityp1,ityp2) .eq. 0) then
cmovc$$$                  if (rg(ibin,ityp1,ityp2) .eq. 0.0d0) then
cmovc$$$                     rio(ibin) = 0.0d0
cmovc$$$                  else
cmovc$$$                     rio(ibin) = (gq(ibin,ityp1,ityp2,1) +
cmovc$$$     $                    gq(ibin,ityp1,ityp2,2)) /
cmovc$$$     $                    rg(ibin,ityp1,ityp2)
cmovc$$$c$$$     $                    ng(ibin,ityp1,ityp2)
cmovc$$$                  endif
cmovc$$$  173          continue
cmovc$$$               if (rdump(iutmp, rio, maxbin) .ne. 0) go to 9980
cmovc$$$            else
cmovc$$$               write(name, '(''gq'', 2i2.2, ''-'', i2.2)') 
cmovc$$$     $              iaid1, iaid2, iaid1
cmovc$$$               write(iutmp, err = 9980, iostat = ioval) 
cmovc$$$c$$$     $              name, rbmn, rbmx, .false., .false.
cmovc$$$     $              name, maxbin, rbmn, rbmx, .false., .false.
cmovc$$$               do 175 ibin = 1, maxbin
cmovc$$$c$$$                  if (ng(ibin,ityp1,ityp2) .eq. 0) then
cmovc$$$                  if (rg(ibin,ityp1,ityp2) .eq. 0.0d0) then
cmovc$$$                     rio(ibin) = 0.0d0
cmovc$$$                  else
cmovc$$$                     rio(ibin) = gq(ibin,ityp1,ityp2,1) / 
cmovc$$$     $                    (rg(ibin,ityp1,ityp2) / 2.0d0)
cmovc$$$c$$$     $                    (ng(ibin,ityp1,ityp2) / 2)
cmovc$$$                  endif
cmovc$$$  175          continue
cmovc$$$               if (rdump(iutmp, rio, maxbin) .ne. 0) go to 9980
cmovc$$$               write(name, '(''gq'', 2i2.2, ''-'', i2.2)') 
cmovc$$$     $              iaid1, iaid2, iaid2
cmovc$$$               write(iutmp, err = 9980, iostat = ioval) name, maxbin,
cmovc$$$     $              rbmn, rbmx, .false., .false.
cmovc$$$               do 177 ibin = 1, maxbin
cmovc$$$c$$$                  if (ng(ibin,ityp1,ityp2) .eq. 0) then
cmovc$$$                  if (rg(ibin,ityp1,ityp2) .eq. 0.0d0) then
cmovc$$$                     rio(ibin) = 0.0d0
cmovc$$$                  else
cmovc$$$                     rio(ibin) = gq(ibin,ityp1,ityp2,2) /
cmovc$$$     $                    (rg(ibin,ityp1,ityp2) / 2.0d0)
cmovc$$$c$$$     $                    (ng(ibin,ityp1,ityp2) / 2)
cmovc$$$                  endif
cmovc$$$  177          continue
cmovc$$$               if (rdump(iutmp, rio, maxbin) .ne. 0) go to 9980
cmovc$$$            endif
cmov               endif
cmov 180        continue
cmov         endif
cmov 190  continue
cmovc***********************************************************************
cmovc     for clusters,
cmovc     calculate the density as a function of r (from the c.o.m.), 
cmovc     normalized by water's liquid-state density:
cmovc***********************************************************************
cmov      if (.not. pbflag) then
cmov         do 220 imol = 1, nmol
cmov            write(name, '(''rho'', i2.2)') imol
cmov            write(iutmp, err = 9980, iostat = ioval)
cmov     $           name, maxbin, rbmn, rbmx, .false., .false.
cmov            do 200 ibin = 1, maxbin
cmov               temp = 1.d0 / (4.d0 * pi * 0.0334273 / 3.d0 * 
cmov     $              ((ibin * rbw) ** 3 - ((ibin - 1) * rbw) ** 3))
cmov               rio(ibin) = rcg(ibin,imol) * temp / nbinst
cmov               riovar(ibin) = rcg2(ibin,imol) * temp * temp / nbinst -
cmov     $              rio(ibin) * rio(ibin)
cmov 200        continue
cmov            if (irdump(iutmp, rio, maxbin) .ne. 0) go to 9980
cmov            if (irdump(iutmp, riovar, maxbin) .ne. 0) go to 9980
cmov 220     continue
cmov      endif
cmov      close(iutmp)
cmov      return
cmov 9980 write(iuout, *) 'sumbin: error writing to ', filnam
cmov      call ioerr(ioval)
cmov      stop
cmov 9990 write(iuout, *) 'sumbin: error opening file ', filnam
cmov      call ioerr(ioval)
cmov      stop
cmov      end
cmovc
c***********************************************************************
c     wrteor writes out an end-of-run file, either at the end of the 
c     run, or earlier, in case of a crash....sjs 10/5/93
c***********************************************************************
c
      subroutine wrteor()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
     $     '/eor.', suffix
      open(iutmp, file = filnam, err = 9990, iostat = ioval)
      write(iutmp, *, err = 9980, iostat = ioval)
     $     iosysv
      write(iutmp, *, err = 9980, iostat = ioval)
c$$$     $     (boxsiz(i), i = 1, 3), rcut
     $     (boxsiz(i), i = 1, 3)
      write(iutmp, *, err = 9980, iostat = ioval)
     $     nmol
      do 300 imol = 1, nmol
         imty = molty(imol)
         write(iutmp, *, err = 9980, iostat = ioval)
     $        imol, molid(imol), nframe(imty), ntag(imty)
         do 290 imind = 1, nframe(imty)
            iatom = iatmol(imol,imind)
            iaty = iatype(imty,imind)
            write(iutmp, *, err = 9980, iostat = ioval)
     $           imind, ident(iatom)
            write(iutmp, *, err = 9980, iostat = ioval)
     $           (pos(iatom,i), i = 1, 3)
            write(iutmp, *, err = 9980, iostat = ioval)
     $           (vel(iatom,i), i = 1, 3)
c$$$            write(iutmp, *, err = 9980, iostat = ioval)
c$$$     $           (vel(iatom,i), i = 1, 3), (acc(iatom,i), i = 1, 3)
            write(iutmp, *, err = 9980, iostat = ioval)
c$$$     $           q(iatom), qvel(iatom), qacc(iatom)
     $           q(iatom), qvel(iatom)
 290     continue
         do 291 itind = 1, ntag(imty)
            imind = itind + nframe(imty)
            iatom = iatmol(imol,imind)
            write(iutmp, *, err = 9980, iostat = ioval)
     $           itind, ident(iatom)
            write(iutmp, *, err = 9980, iostat = ioval)
c$$$     $           q(iatom), qvel(iatom), qacc(iatom)
     $           q(iatom), qvel(iatom)
 291     continue
 300  continue
      write(iutmp, *, err = 9980, iostat = ioval) iocbnv
cmov      write(iutmp, *, err = 9980, iostat = ioval) iobwfl
c$$$      write(iutmp, *, err = 9980, iostat = ioval) nciter
      write(iutmp, *, err = 9980, iostat = ioval) niters
c$$$      write(iutmp, *, err = 9980, iostat = ioval) moving, eqflag
      write(iutmp, *, err = 9980, iostat = ioval) eqflag
      write(iutmp, *, err = 9980, iostat = ioval) nscbas
      write(iutmp, *, err = 9980, iostat = ioval) wvirnx
      write(iutmp, *, err = 9980, iostat = ioval) maxbin
cmov      write(iutmp, *, err = 9980, iostat = ioval) nbinst
cmov      write(iutmp, *, err = 9980, iostat = ioval) apesum, apessq
cmov      write(iutmp, *, err = 9980, iostat = ioval) qpesum, qpessq
      write(iutmp, *, err = 9980, iostat = ioval) rktsum, rktssq
cmov      if (iobwfl) then
cmov         write(iutmp, *, err = 9980, iostat = ioval) 
cmov     $        (ang(i), i = 1, maxbin)
cmov         write(iutmp, *, err = 9980, iostat = ioval) 
cmov     $        (gkr(i), i = 1, maxbin)
cmov         write(iutmp, *, err = 9980, iostat = ioval) 
cmov     $        (nang(i), i = 1, maxbin)
cmov         write(iutmp, *, err = 9980, iostat = ioval) 
cmov     $        (nbondE(i), i = 1, maxbin)
cmov         write(iutmp, *, err = 9980, iostat = ioval) 
cmov     $        (npairE(i), i = 1, maxbin)
cmov         write(iutmp, *, err = 9980, iostat = ioval) 
cmov     $        (nwhb(i), i = 1, mxhbp1)
cmov      endif
cmov      write(iutmp, *, err = 9980, iostat = ioval) 
cmov     $     (((rg(i,j,l), l = 1, j), j = 1, maxaty),i = 1, maxbin)
cmov      write(iutmp, *, err = 9980, iostat = ioval) 
cmov     $     (((rg2(i,j,l), l = 1, j), j = 1, maxaty), i = 1, maxbin)
cmovc$$$      write(iutmp, *) (nHOH(i), i = 1, maxbin)
cmovc$$$      write(iutmp, *) (nrHH(i), i = 1, maxbin)
cmovc$$$      write(iutmp, *) (nrOH(i), i = 1, maxbin)
cmov      write(iutmp, *, err = 9980, iostat = ioval) 
cmov     $     ((nq(i,j), j = 1, maxaty), i = 1, maxbin)
      close(iutmp)
      return
 9980 write(iuout, *) 'wrteor: error writing to ', filnam
      call ioerr(ioval)
      stop
 9990 write(iuout, *) 'wrteor: error opening ', filnam
      call ioerr(ioval)
      end
c
c***********************************************************************
c     ioopen opens all the data files needed for qdyn.  Units 7 and 17 
c     are reserved for quick open/shut use.  Unit 11 is used in
c     Friesner's code (used in J()), and units in the 20s are reserved
c     for debugging....sjs 6/24/93
c***********************************************************************
c
      subroutine ioopen()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      character ctemp*30
c
c$$$      write(filnam,*) scrdir, '/anneal.', suffix
c$$$      open(8, file = filnam, form = 'unformatted', err = 9990,
c$$$     $     iostat = ioval)
c$$$      close(8, status = 'delete')
c$$$      if (qdflag .and. .not. qslvfl) then
c$$$         open(8, file = filnam, form = 'unformatted', status = 'new',
c$$$     $        err = 9990, iostat = ioval)
c$$$         write(iuout, *) 'writing annealing dynamics data to anneal.', 
c$$$     $        suffix
c$$$         if (qifile .eq. 'none') then
c$$$            itemp = int((nqstep + 0.1) / (maxsiz / iabyte))
c$$$            if (itemp .gt. nioint) then
c$$$               nioint = itemp
c$$$               write(iuout, *) 'ioopen: WARNING! had to set nioint to ',
c$$$     $              nioint, ' to avoid having an anneal.', suffix,
c$$$     $              ' larger than ', maxsiz, ' bytes'
c$$$            endif
c$$$            itemp = int(nqstep / maxrec) + 1
c$$$            if (itemp .gt. nioint) then
c$$$               nioint = itemp
c$$$               write(iuout, *) 'ioopen: WARNING! had to set nioint to ',
c$$$     $              nioint, ' to avoid exceeding ', maxrec,
c$$$     $              ' records in anneal.', suffix
c$$$            endif
c$$$            if (nioint .eq. 0) stop 'ioopen: nioint = 0 !'
c$$$            itemp = iabyte * (nqstep / nioint)
c$$$            if (itemp .gt. 2 ** 20) then
c$$$               temp = 10 * itemp / 2 ** 20 / 10
c$$$               write(ctemp, '(f5.1,a)') temp, ' megabytes'
c$$$            else
c$$$               temp = 10 * itemp / 2 ** 10 / 10
c$$$               write(ctemp, '(f5.1,a)') temp, ' kilobytes'
c$$$            endif
c$$$            write(iuout, *) 'anneal.', suffix,
c$$$     $           ' should be no larger than ', ctemp
c$$$         else
c$$$            write(iuout, *) 'anneal.', suffix,
c$$$     $           ' should be small unless the qifile is wrong'
c$$$         endif
c$$$      endif
      write(filnam, '(3a)') scrdir(1:index(scrdir, ' ')-1), 
     $     '/dyn.', suffix
      open(iudyn, file = filnam, form = 'unformatted', err = 9990,
     $     iostat = ioval)
      close(iudyn, status = 'delete')
c$$$      if (dynflg) then
      nitmax = ndstep + 2
c$$$         if (qdflag .and. qifile .eq. 'none' .and. .not. qslvfl)
c$$$     $        nitmax = nitmax + nqstep
      if (vscflg) then 
         nitmax = nitmax + mxsceq
      endif
      if (cbnflg) then
         nitmax = nitmax - niters
      endif
      open(iudyn, file = filnam, form = 'unformatted',
     $     status = 'new', err = 9990, iostat = ioval)
      write(iudyn) iodynv
      itemp = int((nitmax + 0.1d0) / (maxsiz / idbyte))
      if (itemp .gt. nioint) then
         nioint = itemp
         write(iuout, *) 'ioopen: WARNING! had to set nioint to ',
     $        nioint, ' to avoid having a dyn.', suffix, 
     $        ' larger than ', maxsiz, ' bytes'
      endif
      itemp = int(nitmax / maxrec) + 1
c      if (itemp .gt. nioint) then
c         nioint = itemp
c         write(iuout, *) 'ioopen: WARNING!  had to set nioint to ',
c     $        nioint, ' to avoid exceeding ', maxrec,
c     $        ' records in dyn.', suffix
c      endif
c      if (nioint .eq. 0) stop 'ioopen: nioint = 0!'
c      itemp = idbyte * (nitmax / nioint)
c      if (itemp .gt. 2 ** 20) then
c         temp = 10.0d0 * itemp / 2 ** 20 / 10
c         write(ctemp, '(f5.1,a)') temp, ' megabytes'
c      else
c         temp = 10.0d0 * itemp / 2 ** 10 / 10
c         write(ctemp, '(f5.1,a)') temp, ' kilobytes'
c      endif
c      write(iuout, *) 'dyn.', suffix,
c     $     ' should be no larger than ', ctemp
c$$$      write(filnam, *) datdir, '/sysdip.', suffix
c$$$      open(iudip, file = filnam, form = 'unformatted', 
c$$$     $     err = 9990, iostat = ioval)
c$$$      write(iudip) iodipv
c$$$      write(iuout, *) 'writing system dipole data to ', filnam
c$$$      endif
c$$$      if (dbflag) then
c$$$         write(filnam, *) scrdir, '/debug.', suffix
c$$$         open(20, file = filnam, err = 9990, iostat = ioval)
c$$$      endif
      return
 9990 write(iuout, *) 'qdyn: error opening ', filnam
      call ioerr(ioval)
      stop
      end
c***********************************************************************
c     rdeor reads in the eor file for the filname in filnam (should be
c     an eor.SUF file or *.sys file)  rdmodb needs to be called before
c     rdeor, to fill knwmid, nfrmid, and ntamid....sjs 7/22/94
c***********************************************************************
c
      subroutine rdeor(filenm, knwmid, nfrmid, ntamid, anaflg)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      character*80 filenm
      integer*4    nfrmid(maxmid), ntamid(maxmid)
      logical      anaflg, 
cmov     $     tiobwf,
     $             knwmid(maxmid)
c
      write(iuout, *) 'reading startup info from ', filenm
      open(iutmp, file = filenm, err = 9990, iostat = ioval)
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     iovers
      if (iovers .ne. iosysv) then
         write(iuout, *) 'rdeor:  incompatible system file'
         stop
      endif
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
c$$$     $     (boxsiz(i), i = 1, 3), rcut
     $     (boxsiz(i), i = 1, 3)
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     nmol
      if (nmol .gt. maxmol) then
         write(iuout, *) 'rdeor:  Too many molecules in system.'
         stop
      endif
      natoms = 0
      do 300 imol = 1, nmol
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        iimol, iid, nfrtmp, ntatmp
         if (iimol .ne. imol) then
            write(iuout, *) 'rdeor: lost count of molecules at ',
     $           'molecule ', imol
            stop
         endif
         if (iid .gt. maxmid) then
            write(iuout, *) 'rdeor: molecule id ', iid, 
     $           'exceeds maxmid of ', maxmid
            stop
         endif
         if (.not. knwmid(iid)) then
            write(iuout, '(a,i2)') 
     $           'rdeor: don''t recognize molecule id ', iid
            stop
         endif
         molid(imol) = iid
         hasmid(iid) = .true.
         if (nfrtmp .ne. nfrmid(iid)) then
            write(iuout, '(a,i4,a,i2,a,i1,a,i1)') 'rdeor: molecule ',
     $           imol, ' has id ', iid, ' and should have ',
     $           nfrmid(iid), ' frame atoms, not ', nfrtmp
            stop
         endif
         if (ntatmp .ne. ntamid(iid)) then
            write(iuout, '(a,i4,a,i2,a,i1,a,i1)') 'rdeor: molecule ',
     $           imol, ' has id ', iid, ' and should have ',
     $           ntamid(iid), ' tag atoms, not ', ntatmp
            stop
         endif
         if (natoms + nfrtmp + ntatmp .gt. maxatm) then
            write(iuout, '(a,i4,a)') 'rdeor: ', 
     $           natoms + nfrtmp + ntatmp, ' is too many atoms'
            stop
         endif
         do 290 ifr = 1, nfrtmp
            iatom = natoms + ifr
            molec(iatom) = imol
            iatmol(imol,ifr) = iatom
            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $           iifr, ident(iatom)
            if (iifr .ne. ifr) then
               write(iuout, *) 'rdeor: lost count of frame atoms ',
     $              'in molecule ', imol
               stop
            endif
            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $           (pos(iatom,i), i = 1, 3)
            if (stfile .ne. 'none') then
               read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $              (vel(iatom,i), i = 1, 3)
c$$$               read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
c$$$     $              (vel(iatom,i), i = 1, 3), (acc(iatom,i), i = 1, 3)
               read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $              q(iatom), qvel(iatom)
            endif
 290     continue
         natoms = natoms + nfrtmp
         do 291 itg = 1, ntatmp
            imind = nfrtmp + itg
            iatom = natoms + itg
            molec(iatom) = imol
            iatmol(imol,imind) = iatom
            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $           iitg, ident(iatom)
            if (iitg .ne. itg) then
               write(iuout, *) 'rdeor: lost count of tag atoms ',
     $              'in molecule ', imol
               stop
            endif
            if (stfile .ne. 'none') then
               read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
c$$$     $              q(iatom), qvel(iatom), qacc(iatom)
     $              q(iatom), qvel(iatom)
            endif
 291     continue
         natoms = natoms + ntatmp
         numatm(imol) = nfrtmp + ntatmp
 300  continue
      if (cbnflg .or. anaflg) then
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        iovers
         if (iovers .ne. iocbnv) then
            write(iuout, *) 'rdeor:  incompatible version of cbn file'
            stop
         endif
cmov         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmov     $        tiobwf
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
c$$$     $        nciter
     $        niters
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval) 
c$$$     $        moving, eqflag
     $        eqflag
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        nscbas
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        wvirnx
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        nbin
         if (nbin .ne. maxbin) then
            write(iuout, *) 'rdeor: cannot continue binning from the ',
     $           'startup file - maxbin has changed from ', nbin
            stop
         endif
cmov         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmov     $        nbinst
cmov         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmov     $        apesum, apessq
cmov 1       read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmov     $        qpesum, qpessq
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        rktsum, rktssq
cmov         if (tiobwf) then
cmov            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmov     $           (ang(i), i = 1, maxbin)
cmov            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmov     $           (gkr(i), i = 1, maxbin)
cmov            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmov     $           (nang(i), i = 1, maxbin)
cmov            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmov     $           (nbondE(i), i = 1, maxbin)
cmov            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmov     $           (npairE(i), i = 1, maxbin)
cmov            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmov     $           (nwhb(i), i = 1, mxhbp1)
cmov         endif
cmov         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmov     $        (((rg(i,j,l), l = 1, j), j = 1, maxaty), i = 1, maxbin)
cmov         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmov     $        (((rg2(i,j,l), l = 1, j), j = 1, maxaty), i = 1, maxbin)
cmovc$$$         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmovc$$$     $        (nHOH(i), i = 1, maxbin)
cmovc$$$         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmovc$$$     $        (nrHH(i), i = 1, maxbin)
cmovc$$$         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmovc$$$     $        (nrOH(i), i = 1, maxbin)
cmov         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
cmov     $        ((nq(i,j), j = 1, maxaty), i = 1, maxbin)
cmov      else
cmov         nbinst = 0
      endif
      close(iutmp)
      return
 9970 write(iuout, *) 'rdeor: error reading from ', filenm
      call ioerr(ioval)
      stop
 9980 write(iuout, *) 'rdeor: end of file in ', filenm
      call ioerr(ioval)
      stop
 9990 write(iuout, *) 'rdeor: error opening ', filenm
      call ioerr(ioval)
      stop
      end
c
c***********************************************************************
c     rdmodb reads in the info about molecule ids in models.dat, but
c     only remembers a few things.  The full details will be read in 
c     again in rdmodd, when it is known which models can be ignored.
c     ...sjs 7/28/94
c***********************************************************************
c
      subroutine rdmodb(knwmid, nfrmid, ntamid)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      character cjunk*30
      integer   nfrmid(maxmid), ntamid(maxmid)
      logical   ljunk, rdsome, 
     $          knwmid(maxmid)
c
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
c***********************************************************************
c     read from models.dat:
c***********************************************************************
      write(filnam, '(2a)') datdir(1:index(datdir, ' ')-1), 
     .     '/models.dat'
      open(iutmp, file = filnam, err = 9990, iostat = ioval)
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     iovers
      if (iovers .ne. iomodv) then
         write(iuout, *) 'rdmodb:  incompatible models file'
         stop
      endif
 200  continue
      read(iutmp, *, end = 300, err = 9970, iostat = ioval)
     $     imid, cjunk
      if (imid .gt. maxmid) then
         write(iuout, '(a,i2,a)') 'rdmodb: molecule id ', imid, 
     $        ' is too big'
         write(iuout, *) imid, maxmid
         stop
      endif
      knwmid(imid) = .true.
      rdsome = .true.
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     natyid
      do 205 iaty = 1, natyid
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        iiaty, ijunk
         if (iiaty .ne. iaty) then
            write(iuout, '(a, i2, a, i2)') 
     $           'rdmodb: lost count of atom types at number ', iaty, 
     $           ' in molecule id ', imid
            stop
         endif
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        ljunk
         if (ljunk) then
            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $           rjunk, rjunk
c$$$            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
c$$$     $           ljunk
c$$$            if (ljunk) then
c$$$               read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
c$$$     $              rjunk, rjunk, rjunk
c$$$            endif
         endif
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        ljunk
         if (ljunk) then
            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $           rjunk
         endif
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        ljunk
         if (ljunk) then
            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $           rjunk, rjunk
         else
            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $           ljunk
            if (ljunk) then
               read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $              rjunk
               read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $              ljunk
               if (ljunk) then
                  read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $                 rjunk
               endif
            endif
         endif
c$$$         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
c$$$     $        ljunk
c$$$         if (ljunk) then
c$$$            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
c$$$     $           rjunk, rjunk, rjunk, rjunk
c$$$         endif
 205  continue
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     nfrmid(imid), ntamid(imid)
      if (nfrmid(imid) .gt. 4 .or. 
     $     (nfrmid(imid) .lt. 2 .and. ntamid(imid) .ne. 0)) then
         write(iuout, '(a,i2,a,i1,a)') 'rdmodb: molecule id ', imid, 
     $        ' has ', nfrmid(imid), ' frame atoms'
         stop
      endif
      if (nfrmid(imid) + ntamid(imid) .gt. mxatml) then
         write(iuout, '(a,i2,a,i2,a)') 'rdmodb: molecule id ', imid,
     $        ' has ', nfrmid(imid) + ntamid(imid), 
     $        ' atoms, which is too many'
         stop
      endif
      do 210 ifr = 1, nfrmid(imid)
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        iifr, iaty
         if (iifr .ne. ifr) then
            write(iuout, *) 'rdmodb: lost count of frame atoms in ',
     $           'molecule id ', imid
            stop
         endif
 210  continue
      do 220 ita = 1, ntamid(imid)
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        iiat, iita, ijunk, (rjunk, i = 1, nfrmid(imid))
         if (iiat .ne. iita + nfrmid(imid) .or. iita .ne. ita) then
            write(iuout, *) 'rdmodb: lost count of tag atoms in ',
     $           'molecule id ', imid
            stop
         endif
 220  continue
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     nconid
      do 230 iconi = 1, nconid
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        ijunk, ijunk, rjunk
 230  continue
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     ndrdid
      do 233 idrdid = 1, ndrdid
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
c$$$     $        ijunk, ijunk, rjunk
     $        ijunk, idrspi, rjunk
 233  continue
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     rjunk
      go to 200
 300  continue
      if (.not. rdsome) go to 9980
      close(iutmp)
      return
 9970 write(iuout, *) 'rdmodb: error reading from ', filnam
      call ioerr(ioval)
      stop
 9980 write(iuout, *) 'rdmodb: unexpected end of file in ', filnam
      call ioerr(ioval)
      stop
 9990 write(iuout, *) 'rdmodb: error opening ', filnam
      call ioerr(ioval)
      stop
      end
c
c***********************************************************************
c     rdmodd reads in the parameters for all of the molecule ids in
c     models.dat, and remembers them for all molecule ids being used in
c     the simulation indexing them by molecule type....sjs 7/28/94
c***********************************************************************
c
      subroutine rdmodd()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      character tmodnm*30
      integer   iaida(maxaty),
     $          iadrd1(mxatml/2), iadrd2(mxatml/2),
     $          iatya(mxatml),
     $          iacon1(maxcon), iacon2(maxcon)
      logical   rdsome,
c$$$     $          hlja(maxaty), hmsa(maxaty), hfqa(maxaty), hr4a(maxaty),
     $          hlja(maxaty), hmsa(maxaty), hfqa(maxaty),
c$$$     $          hqda(maxaty), hrqa(maxaty)
     $          hrqa(maxaty), hztaa(maxaty)
      real*8    sprcon(mxatml/2),
     $          epsa(maxaty), siga(maxaty), dfqa(maxaty),
     $          rmsa(maxaty), chia(maxaty), ztaa(maxaty),
c$$$     $          r4ka(maxaty), gexa(maxaty), gcoa(maxaty), gpoa(maxaty),
c$$$     $          sg0a(maxaty), sg1a(maxaty), q0a(maxaty),
     $          rcons(maxcon),
     $          tgca(mxtgat,mxfrat)
c
c***********************************************************************
c     initialize some stuff:
c***********************************************************************
      naty = 0
      write(filnam, '(2a)') datdir(1:index(datdir, ' ')-1), 
     .     '/models.dat'
      open(iutmp, file = filnam, err = 9990, iostat = ioval)
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     iovers
c***********************************************************************
c     read in all of the data for a single molecule id:
c***********************************************************************
 200  continue
      read(iutmp, *, end = 300, err = 9970, iostat = ioval)
     $     imid, tmodnm
      rdsome = .true.
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     natyid
      do 205 iaty = 1, natyid
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        ijunk, iaida(iaty)
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        hlja(iaty)
         if (hlja(iaty)) then
            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $           epsa(iaty), siga(iaty)
c$$$            if (imid .eq. 6) then
c$$$               write(iuout, *) 'Drude LJ: ', epsa(iaty), siga(iaty)
c$$$            endif
c$$$            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
c$$$     $           hqda(iaty)
c$$$            if (hqda(iaty)) then
c$$$               read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
c$$$     $              sg0a(iaty), sg1a(iaty), q0a(iaty)
c$$$            endif
         endif
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        hmsa(iaty)
         if (hmsa(iaty)) then
            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $           rmsa(iaty)
c$$$            if (imid .eq. 6 .and. iaty .eq. 2) then
c$$$               write(iuout, *) 'Drude dm: ', rmsa(iaty)
c$$$            endif
         endif
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        hfqa(iaty)
         if (hfqa(iaty)) then
            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $           chia(iaty), ztaa(iaty)
            hrqa(iaty) = .false.
         else
            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $           hrqa(iaty)
            if (hrqa(iaty)) then
               read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $              dfqa(iaty)
               read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $              hztaa(iaty)
               if (hztaa(iaty)) then
                  read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $                 ztaa(iaty)
c$$$                  if (imid .eq. 6) then
c$$$                     write(iuout, *) 'Drude zeta: ', ztaa(iaty)
c$$$                  endif
               endif
            endif
         endif
c$$$         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
c$$$     $        hr4a(iaty)
c$$$         if (hr4a(iaty)) then
c$$$            read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
c$$$     $           r4ka(iaty), gexa(iaty), gcoa(iaty), gpoa(iaty)
c$$$         endif
 205  continue
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     nfrid, ntaid
      do 210 ifr = 1, nfrid
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        ijunk, iatya(ifr)
         if (.not. hmsa(iatya(ifr))) then
            write(iuout, '(a, i2, a)') 'rdmodd: molecule id ', 
     $           imid, ' has a massless frame atom'
            stop
         endif
 210  continue
      do 220 ita = 1, ntaid
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        ijunk, ijunk, iatya(nfrid + ita), 
     $        (tgca(ita,i), i = 1, nfrid)
         temp = 1.
         do 215 ifr = 2, nfrid
            temp = temp - tgca(ita,ifr)
 215     continue
         temp = temp - tgca(ita,1)
         if (temp * temp .gt. rattol) then
            write(iuout, *) 'rdmodd: tagcofs do not add to 1 in ',
     $           'molecule id ', imid
            stop
         endif
 220  continue
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     nconid
c$$$      if (nconid .ne. nfrid * (nfrid - 1) / 2) then
c$$$         write(iuout, *) 'rdmodd: flexible frames are not allowed'
c$$$         stop
c$$$      endif
      do 230 iconi = 1, nconid
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        iacon1(iconi), iacon2(iconi), rcons(iconi)
 230  continue
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     ndrdid
      do 233 idrdid = 1, ndrdid
         read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $        iadrd1(idrdid), iadrd2(idrdid), sprcon(idrdid)
c$$$         if (imid .eq. 6) then
c$$$            write(iuout, *) 'Drude mw^2:', sprcon(idrdid)
c$$$         endif
 233  continue
      read(iutmp, *, end = 9980, err = 9970, iostat = ioval)
     $     qmolt
c***********************************************************************
c     if it's a molecule that we care about, set some stuff, indexing
c     it by molecule type:
c***********************************************************************
      if (hasmid(imid)) then
         imty = mtymid(imid)
         modnam(imty) = tmodnm
         nframe(imty) = nfrid
         ntag(imty) = ntaid
         natmty(imty) = nfrid + ntaid
         hasljt(imty) = .false.
         nljmol(imty) = 0
         nmamol(imty) = 0
         do 235 iiaty = naty + 1, naty + natyid
            isfqat(iiaty) = .false.
            isqat(iiaty) = .false.
            isljat(iiaty) = .false.
            ismasv(iiaty) = .false.
            isztat(iiaty) = .false.
c$$$            isqdlj(iiaty) = .false.
 235     continue
         hasqt(imty) = .false.
         hasfqt(imty) = .false.
         hasrqt(imty) = .false.
         hasdrd(imty) = .false.
c$$$         hasr4t(imty) = .false.
         nqmol(imty) = 0
         nfqmol(imty) = 0
         nrqmol(imty) = 0
         ndrmol(imty) = 0
c$$$         nr4mol(imty) = 0
c$$$         ndrmol(imty) = ndrdid
         qmol(imty) = qmolt
         do 245 iaty = 1, natyid
            if (iaty .le. nfrid) then
               isfrat(naty + iaty) = .true.
            else
               isfrat(naty + iaty) = .false.
            endif
            ideaty(naty + iaty) = iaida(iaty)
            imtaty(naty + iaty) = imty
            if (hlja(iaty)) then
               trmlj(naty + iaty,naty + iaty,1) = epsa(iaty)
               trmlj(naty + iaty,naty + iaty,2) = 8. * epsa(iaty)
               trmlj(naty + iaty,naty + iaty,3) = 24. * epsa(iaty)
               trmlj(naty + iaty,naty + iaty,4) = siga(iaty)
c$$$               if (hqda(iaty)) then
c$$$                  isqdlj(naty + iaty) = .true.
c$$$                  trmlj(naty + iaty,naty + iaty,4) = sg0a(iaty)
c$$$                  trmlj(naty + iaty,naty + iaty,5) = sg1a(iaty)
c$$$                  trmlj(naty + iaty,naty + iaty,6) = q0a(iaty)
c$$$                  trmlj(naty + iaty,naty + iaty,11) = 
c$$$     $                 sg1a(iaty) * sg1a(iaty)
c$$$               endif
               eps6 = trmlj(naty + iaty,naty + iaty,4) ** 3
               eps6 = eps6 * eps6
               trmlj(naty + iaty,naty + iaty,7) = 
     $              eps6 * trmlj(naty + iaty,naty + iaty,2)
               trmlj(naty + iaty,naty + iaty,8) = 
     $              3. * trmlj(naty + iaty,naty + iaty,7)
               trmlj(naty + iaty,naty + iaty,9) = 
     $              eps6 * trmlj(naty + iaty,naty + iaty,7)
               trmlj(naty + iaty,naty + iaty,10) = 
     $              6. * trmlj(naty + iaty,naty + iaty,9)
c$$$c***********************************************************************
c$$$c     The mixed-LJ interactions are calculated using the additive
c$$$c     combining rule for sigma (geometric for epsilon):
c$$$c***********************************************************************
c$$$               do 240 iaty2 = 1, iaty - 1
c$$$                  if (hlja(iaty)) then
c$$$                     trmlj(naty + iaty,naty + iaty2,1) = 
c$$$     $                    sqrt(epsa(iaty) *
c$$$     $                    epsa(iaty2))
c$$$                     trmlj(naty + iaty2,naty + iaty,1) = 
c$$$     $                    trmlj(naty + iaty,naty + iaty2,1)
c$$$                     trmlj(naty + iaty,naty + iaty2,2) = 
c$$$     $                    8. * trmlj(naty + iaty,naty + iaty2,1)
c$$$                     trmlj(naty + iaty2,naty + iaty,2) = 
c$$$     $                    trmlj(naty + iaty,naty + iaty2,2)
c$$$                     trmlj(naty + iaty,naty + iaty2,3) = 
c$$$     $                    24. * trmlj(naty + iaty,naty + iaty2,1)
c$$$                     trmlj(naty + iaty2,naty + iaty,3) = 
c$$$     $                    trmlj(naty + iaty,naty + iaty2,3)
c$$$                     trmlj(naty + iaty,naty + iaty2,4) = 0.5 *
c$$$     $                    (siga(iaty) + siga(iaty2))
c$$$                     trmlj(naty + iaty2,naty + iaty,4) = 
c$$$     $                    trmlj(naty + iaty,naty + iaty2,4)
c$$$                     if (hqda(iaty2) .or. hqda(iaty1)) then
c$$$                        write(iuout, *) 'don''t know the trmlj ',
c$$$     $                       'factors for q-dep heterogeneous LJ'
c$$$                        stop
c$$$                     endif
c$$$                     eps6 = trmlj(naty + iaty,naty + iaty2,4) ** 3
c$$$                     eps6 = eps6 * eps6
c$$$                     trmlj(naty + iaty,naty + iaty2,7) = 
c$$$     $                    eps6 * trmlj(naty + iaty,naty + iaty2,2)
c$$$                     trmlj(naty + iaty2,naty + iaty,7) = 
c$$$     $                    trmlj(naty + iaty,naty + iaty2,7)
c$$$                     trmlj(naty + iaty,naty + iaty2,8) = 
c$$$     $                    3. * trmlj(naty + iaty,naty + iaty2,7)
c$$$                     trmlj(naty + iaty2,naty + iaty,8) = 
c$$$     $                    trmlj(naty + iaty,naty + iaty2,8)
c$$$                     trmlj(naty + iaty,naty + iaty2,9) = 
c$$$     $                    eps6 * trmlj(naty + iaty,naty + iaty2,7)
c$$$                     trmlj(naty + iaty2,naty + iaty,9) = 
c$$$     $                    trmlj(naty + iaty,naty + iaty2,9)
c$$$                     trmlj(naty + iaty,naty + iaty2,10) = 
c$$$     $                    6. * trmlj(naty + iaty,naty + iaty2,9)
c$$$                     trmlj(naty + iaty2,naty + iaty,10) = 
c$$$     $                    trmlj(naty + iaty,naty + iaty2,10)
c$$$                  endif
c$$$ 240           continue
               isljat(naty + iaty) = .true.
            endif
            if (hmsa(iaty)) then
               ismasv(naty + iaty) = .true.
               atmass(naty + iaty) = rmsa(iaty)
            endif
            if (hfqa(iaty)) then
               isfqat(naty + iaty) = .true.
               isqat(naty + iaty) = .true.
               isztat(naty + iaty) = .true.
               trmq(naty + iaty,1) = chia(iaty)
               trmq(naty + iaty,2) = ztaa(iaty)
            endif
            if (hrqa(iaty)) then
               defq(naty + iaty) = dfqa(iaty)
               isqat(naty + iaty) = .true.
               if (hztaa(iaty)) then
                  isztat(naty + iaty) = .true.
                  trmq(naty + iaty,2) = ztaa(iaty)
               endif
            endif
c$$$            if (hr4a(iaty)) then
c$$$               trmr4(naty + iaty,1) = gpoa(iaty)
c$$$               trmr4(naty + iaty,2) = gexa(iaty)
c$$$               trmr4(naty + iaty,3) = gcoa(iaty)
c$$$               trmr4(naty + iaty,4) = 2. * r4ka(iaty)
c$$$               trmr4(naty + iaty,5) = 2. * trmr4(naty + iaty,4)
c$$$            endif
 245     continue
         do 250 imind = 1, natmty(imty)
            iaty = iatya(imind)
            idemty(imty,imind) = iaida(iaty)
            iatype(imty,imind) = naty + iaty
            if (hlja(iaty)) then
               hasljt(imty) = .true.
               nljmol(imty) = nljmol(imty) + 1
               iljmol(imty,nljmol(imty)) = imind
            endif
            if (hmsa(iaty)) then
               nmamol(imty) = nmamol(imty) + 1
               imamol(imty,nmamol(imty)) = imind
            endif
            if (hfqa(iaty)) then
               hasfqt(imty) = .true.
               nfqmol(imty) = nfqmol(imty) + 1
               ifqmol(imty,nfqmol(imty)) = imind
            endif
            if (hrqa(iaty)) then
               hasrqt(imty) = .true.
               nrqmol(imty) = nrqmol(imty) + 1
               irqmol(imty,nrqmol(imty)) = imind
            endif
            if (hfqa(iaty) .or. hrqa(iaty)) then
               hasqt(imty) = .true.
               nqmol(imty) = nqmol(imty) + 1
               iqmol(imty,nqmol(imty)) = imind
            endif
c$$$            if (hr4a(iaty)) then
c$$$               hasr4t(imty) = .true.
c$$$               nr4mol(imty) = nr4mol(imty) + 1
c$$$               ir4mol(imty,nr4mol(imty)) = imind
c$$$            endif
 250     continue
         do 253 ifind = 1, nframe(imty)
            tagcof(imty,ifind,ifind) = 1.
            do 252 ifind2 = 1, ifind - 1
               tagcof(imty,ifind,ifind2) = 0.
               tagcof(imty,ifind2,ifind) = 0.
 252        continue
 253     continue
         do 255 itind = 1, ntag(imty)
            imind = nframe(imty) + itind
            do 254 ifind = 1, nframe(imty)
               tagcof(imty,imind,ifind) = tgca(itind,ifind)
 254        continue
 255     continue
         naty = naty + natyid
         if (naty .gt. maxaty) then
            write(iuout, '(a,i2,a)') 'rdmodd: ', naty, 
     $           ' is too many atom types'
            stop
         endif
         ncons(imty) = nconid
         do 260 iconi = 1, nconid
            icons(imty,iconi,1) = iacon1(iconi)
            icons(imty,iconi,2) = iacon2(iconi)
            r2cons(imty,iconi) = rcons(iconi) * rcons(iconi)
 260     continue
         ndrmol(imty) = ndrdid
         do 270 idrdid = 1, ndrdid
            hasdrd(imty) = .true.
            idrmol(imty,idrdid,1) = iadrd1(idrdid)
            idrmol(imty,idrdid,2) = iadrd2(idrdid)
            trmdrd(imty,idrdid) = sprcon(idrdid)
 270     continue
      endif
      go to 200
 300  continue
      if (.not. rdsome) go to 9980
      close(iutmp)
      return
 9970 write(iuout, *) 'rdmodd: error reading from file ', filnam
      call ioerr(iostat)
      stop
 9980 write(iuout, *) 'rdmodd: unexpected end of file in ', filnam
      call ioerr(iostat)
      stop
 9990 write(iuout, *) 'rdmodd: error opening ', filnam
      call ioerr(iostat)
      stop
      end
c
c***********************************************************************
c     rddet reads in the simulation details from the file passed to it
c     (usually details.dat, qdyn.in, or simdat.SUF)....sjs 8/9/94
c***********************************************************************
c     now passed an open filehandle.  file must be opened and closed
c     outside of this routine....sjs 9/8/96
c***********************************************************************
c
      subroutine rddet(ihandl, icjid)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
c$$$      character filenm*80
      integer*4 icjid(mxcstj,3)
c
c$$$      open(iutmp, file = filenm, err = 9990, iostat = ioval)
      read(ihandl, *, end = 9980, err = 9970, iostat = ioval)
     $     iovers
      if (iovers .ne. iosimv) then
         write(iuout, *) 'rddet:  incompatible simulation files'
         stop
      endif
      read(ihandl, *, end = 9980, err = 9970, iostat = ioval)
     $     datdir, scrdir
      read(ihandl, *, end = 9980, err = 9970, iostat = ioval) 
     $     stfile, syfile, qifile
      read(ihandl, *, end = 9980, err = 9970, iostat = ioval) rinfo
      read(ihandl, *, end = 9980, err = 9970, iostat = ioval) 
     $     ctflag, qslvfl, nqsolv, qmass
      read(ihandl, *, end = 9980, err = 9970, iostat = ioval) 
     $     ndstep, pbflag, dtbig, nl, nlf, nhf
      read(ihandl, *, end = 9980, err = 9970, iostat = ioval)
     $     dslvfl, ndsolv
      read(ihandl, *, end = 9980, err = 9970, iostat = ioval) 
     $     fcutr, fcutd, fcutb, fnearr, fneard, iplint, hvywt
      read(ihandl, *, end = 9980, err = 9970, iostat = ioval) 
     $     interJ, intraJ, ncint
      do 100 icint = 1, ncint
         read(ihandl, *, end = 9980, err = 9970, iostat = ioval)
     $        icjid(icint,1), icjid(icint,2), icjid(icint,3)
  100 continue
      read(ihandl, *, end = 9980, err = 9970, iostat = ioval)
     $     ewlflg, ewlkfl, ewsrfl, ewlkap, ewlkmx
      read(ihandl, *, end = 9980, err = 9970, iostat = ioval) T, Ttol,
     $     nscint, stint, vscflg, nsceql, mxsceq, qT, qTrang, qvscfl
      read(ihandl, *, end = 9980, err = 9970, iostat = ioval) cbnflg, 
     $     dbflag, iseed, nioint, hbonde, 
cmov     $     gstinf, 
     $     ocflag,
     $     Efield(1), Efield(2), Efield(3)
      read(ihandl, *, end = 9980, err = 9970, iostat = ioval) 
     $     ioflg, binflg, nbnint, 
cmov     $     iobwfl, 
     $     suffix
c$$$      close(iutmp)
      return
 9970 write(iuout, *) 'rddet: error reading from ', filnam
      call ioerr(ioval)
      stop
 9980 write(iuout, *) 'rddet: unexpected end of file in ', filnam
      call ioerr(ioval)
      stop
      end
c
c***********************************************************************
c     wrtfor is a debugging subroutine
c***********************************************************************
c
      subroutine wrtfor(idist, imass)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      write(filnam, '(2a)') datdir(1:index(datdir, ' ')-1), '/pos.xyz'
      open(iutmp, file = filnam)
      write(iutmp, *) natoms
      do 500 iatom = 1, natoms
         write(iutmp, *) (pos(iatom,i), i = 1, 3), q(iatom)
 500  continue
      close(iutmp)
      write(filnam, '(2a)') datdir(1:index(datdir, ' ')-1), 
     .     '/forces.xyz'
      open(iutmp, file = filnam)
      write(iutmp, *) natoms
      do 600 iatom = 1, natoms
         write(iutmp, *) (pos(iatom,i), i = 1, 3), 
     $        (1000.d0 * force(iatom,i,idist,imass), i = 1, 3)
 600  continue
      close(iutmp)
      return
      end
c
c***********************************************************************
c     wrttfr is a debugging subroutine
c***********************************************************************
c
      subroutine wrttfr()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      write(filnam, '(2a)') datdir(1:index(datdir, ' ')-1), '/pos.xyz'
      open(iutmp, file = filnam)
      write(iutmp, *) natoms
      do 500 iatom = 1, natoms
         write(iutmp, *) (pos(iatom,i), i = 1, 3), q(iatom)
 500  continue
      close(iutmp)
      write(filnam, '(2a)') datdir(1:index(datdir, ' ')-1), 
     .     '/forces.xyz'
      open(iutmp, file = filnam)
      write(iutmp, *) natoms
      do 600 iatom = 1, natoms
         write(iutmp, *) (pos(iatom,i), i = 1, 3), (1000.d0 * 
     $        (force(iatom,i,inear,ilight) + 
     $        force(iatom,i,inear,imixed) +
     $        force(iatom,i,inear,iheavy) +
     $        force(iatom,i,ifar,ilight) +
     $        force(iatom,i,ifar,imixed) +
     $        force(iatom,i,ifar,iheavy)), i = 1, 3)
 600  continue
      close(iutmp)
      return
      end

