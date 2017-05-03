      program analyz
c***********************************************************************
c     This program reads, analyzes, and rewrites the information in the
c     various data files spit out by qdyn.
c     Input:  bins.SUF, dyn.SUF (unformatted)
c             simdat.SUF, eor.SUF (formatted)
c     Output: fulle.SUF, mu.SUF, qe.SUF, q.SUF, p.SUF,
c             temps.SUF, trnord.SUF,
c             other files specified in bins.SUF
c     ...sjs 9/15/92
c***********************************************************************
c     It only needs the files from the most recent of a string of
c     cbn'ed runs.  The data (diffusion constant, dielectric constant)
c     that need total trajectory histories are analyzed in analz2.
c     ...sjs 9/27/94
c***********************************************************************
c     The new version of qdyn does no binning of its own, so analyz
c     now has to do everything....sjs 3/30/94
c***********************************************************************
c
      include 'implic'      
      include 'genpar'
      include 'qpar'
      include 'anapar'
c$$$      include 'commons'
      include 'anacom'
c     
c$$$      parameter (maxcbn = 100)
c$$$      parameter (tcmax = 200000.0d0)
c$$$c
      character    name*12, suff2*10, filenm*80, fildir*80
      character*20 suffl(mdatfl)
      integer*4    nfrmid(maxmid), ntamid(maxmid),
     $     ntres(maxcbn),
     $     nbrsnc(maxmol,maxmol), nbrout(maxmol,maxmol)
c$$$     $     numtc(maxrec)
      logical      accflg, coorfl, lspcok, nonzfl, rsqflg,
     $     knwmid(maxmid),
     $     nbrfl(maxmol,maxmol), xnbrfl(maxmol,maxmol)
      real*4       rio(maxdat), riovar(maxdat)
      real*8       dumvec(1), dumvc2(1),
     $     rltvec(3), tcsmu(3),
     $     smuave(3), smusum(3), smussq(3),
     $     gexpt(maxbin), xvalvc(maxbin), onevec(maxbin), xvec(maxbin),
     $     x2vec(maxbin), x3vec(maxbin), x4vec(maxbin),
     $     yvec(maxbin), xyvec(maxbin), xxyvec(maxbin),
c$$$     $     phi(maxrec), phissq(maxrec), rms(maxrec), rmsstd(maxrec),
c$$$     $     rint(maxrec),
     $     oldq(maxatm),
     $     dipole(maxmol,3),
     $     splcof(maxbin,4),
     $     smrio(maxdat),
c$$$     $     csmu(maxrec,3),
     $     oldpos(maxatm,3), tcpos(maxatm,3)
c$$$     $     posxt(maxatm,3,maxrec)
c     
      logical   anaflg
c
      parameter (linst = 10000)
c     
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
      lspcok = .true.
      nmu = 0
      npr = 0
      psum = 0.0
      pssq = 0.0
c$$$      do 100 i = 1, maxrec
c$$$         rint(i) = dble(i)
c$$$  100 continue
c$$$      symass = 0.
      do 110 imid = 1, maxmid
         knwmid(imid) = .false.
 110  continue
      do 130 imol2 = 1, maxmol
         do 120 imol1 = 1, maxmol
            nbrfl(imol1,imol2) = .false.
            xnbrfl(imol1,imol2) = .false.
 120     continue
 130  continue
      nntres = 0
      do 140 ibin = 1, maxcbn
         ntres(maxcbn) = 0
 140  continue
c***********************************************************************
c     Read in the input file analyz.in
c***********************************************************************
      write(filnam, *) 'analyz.in'
      open(iutmp, file = filnam, err = 9980, iostat = ioval)
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) fildir
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) nsuff 
      if (nsuff .gt. mdatfl) then
         write(iuout, *) 'analyz: maximum number of data files ',
     $        '(mdatfl) exceeded'
         stop
      endif
      do 200 idatfl = 1, nsuff
         read(iutmp, *, err = 9980, end = 9991, iostat = ioval)
     $        suffl(idatfl)
 200  continue
c$$$      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) nmuint
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) nsqint
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) coorfl
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) 
     $     rcoord, tout, ncbin
c$$$c***********************************************************************
c$$$c     Ask the user which data file to read (I wish FORTRASH took
c$$$c     command-line arguments...):
c$$$c***********************************************************************
c$$$      write(iuout, *) 'Which suffix?'
c$$$      read(iuout, *) suffix
c***********************************************************************
c     Read in the header stuff and initialize related things:
c***********************************************************************
      write(filenm,*) datdir,'simdat.',suffl(1)
      write(iuout, *) 'reading header info from simdat.',suffl(1)
      call rddet(filenm)
      call chkdet
      call setdet
      call rdmodb(knwmid, nfrmid, ntamid)
      write(stfile, *) 'eor.', suffl(1)
      write(filnam, *) datdir, 'eor.', suffl(1)
      anaflg = .true.
      call rdeor(filnam, knwmid, nfrmid, ntamid, anaflg)
      call setsys(anaflg)
      call rdmodd
      call setmod
      if (stfile .eq. 'none' .and. syfile .eq. 'none') then
         write(iuout, *)
     $     'analyz: WARNING! system.dat may have changed...'
      endif
      if (qdflag .and. ctflag) then
         write(iuout, *)
     $     'analyz: WARNING! dipole calcs may be nonsense'
      endif
c$$$      if (dynflg) then
         iterio = ndstep / 20
c$$$      else
c$$$         iterio = nqstep / 20
c$$$      endif
      iterio = (iterio / nioint) * nioint
      if (iterio .lt. 100) iterio = (100 / nioint + 1) * nioint
      rkTinv = rkinv / T
c$$$      if (dynflg) then
      if (.not. cubflg) write(iuout, *) 'analyz: WARNING! box is ', 
     $     'not cubic, translational order parameter may be useless'
c$$$      endif
      rktavg = rktsum / nbinst
      nout = int(tout / dt)
c***********************************************************************
c     Set up output files:
c***********************************************************************
      if (qdflag) then
         write(filnam, *) datdir, 'qe.', suffl(nsuff)
         open(12, file = filnam, err = 9980, iostat = ioval)
         write(iuout, *) 'writing charge energies to qe.', suffl(nsuff)
      endif
c$$$      if (qdflag .or. dynflg) then
      write(filnam,*) datdir,'temps.',suffl(nsuff)
      open(11, file = filnam, err = 9980, iostat = ioval)
      write(iuout, *) 'writing temperatures to temps.', suffl(nsuff)
c$$$      endif
c$$$c***********************************************************************
c$$$c     Analyze annealing data:
c$$$c***********************************************************************
      if (qdflag) then
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
         qtfac = 2.0d0 / (nqdof * k)
c$$$c***********************************************************************
c$$$c     Read in the data and spit it back out:
c$$$c***********************************************************************
c$$$         if (.not. qslvfl) then
c$$$            write(filnam,*) scrdir,'anneal.',suffix
c$$$            write(iuout,*) 'reading from anneal.',suffix
c$$$            open(8, file = filnam, status = 'old', form = 'unformatted',
c$$$     $           err = 9980, iostat = ioval)
c$$$  410       continue
c$$$            if (niters .eq. (niters / iterio) * iterio)
c$$$     $           write(iuout, *) 'Timestep ', niters
c$$$            ioval = iurdmp(8, rio, iabyte / nspbyt)
c$$$            if (ioval .eq. -1) then
c$$$               go to 480
c$$$            else if (ioval .ne. 0) then
c$$$               go to 9990
c$$$            endif
c$$$            do 420 iatom = 1,natoms
c$$$               q(iatom) = rio(iatom)
c$$$  420       continue
c$$$            do 430 iatom = 1,natoms
c$$$               qvel(iatom) = rio(natoms + iatom)
c$$$  430       continue
c$$$            qke = rio(2 * natoms + 1)
c$$$            qpe = rio(2 * natoms + 2)
c$$$            niters = rio(2 * natoms + 3)
c$$$            write(10, *) niters * dt, qke * fJm2kc / natoms,
c$$$     $           qpe * fJm2kc / natoms, (qke + qpe) * fJm2kc / natoms
c$$$            write(11, *) niters * dt, 0.0d0, qke * qtfac, T
c$$$            go to 410
c$$$  480       continue
c$$$            close(8)
c$$$         endif
c$$$      else
c$$$         write(iuout,*) 'No charge dynamics'
      endif
c***********************************************************************
c     Analyze the real dynamics data:
c***********************************************************************
c$$$      if (dynflg) then
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
      if (nmol .eq. 32) then 
         ncell = 2
      else if (nmol .eq. 108) then
         ncell = 3
      else if (nmol .eq. 256) then
         ncell = 4
      else
         write(iuout,*) 'not FCC - translational order parameter',
     $        ' may be meaningless'
         ncell = 1
      endif
      rtfac = 2.0d0 / (nrdof * k)
      vol = boxsiz(1) * boxsiz(2) * boxsiz(3)
      gkfac = 4.0d0 * pi / 9.0d0 * rkinv / vol * epsinv
      rltvec(1) = 1 * 2 * pi * ncell/ boxsiz(1) * (2 * nmol)**0
      rltvec(2) = -1 * 2 * pi * ncell / boxsiz(2)
      rltvec(3) = -1 * 2 * pi * ncell / boxsiz(3)
c$$$      do 600 irec = 1, maxrec
c$$$         numtc(irec) = 0
c$$$         phi(irec) = 0.0d0
c$$$         rms(irec) = 0.0d0
c$$$         rmsstd(irec) = 0.0d0
c$$$         csmu(irec,1) = 0.0d0
c$$$         csmu(irec,2) = 0.0d0
c$$$         csmu(irec,3) = 0.0d0
c$$$ 600  continue
c***********************************************************************
c     Open the data files:
c***********************************************************************
      write(filnam,*) datdir,'fulle.',suffl(nsuff)
      open(20, file = filnam, err = 9980, iostat = ioval)
      write(iuout, *) 'writing full energies to fulle.', suffl(nsuff)
      write(filnam,*) datdir,'trnord.',suffl(nsuff)
      open(21, file = filnam, err = 9980, iostat = ioval)
      write(iuout, *) 'writing equilibration data to trnord.', 
     $     suffl(nsuff)
c$$$      write(filnam, *) datdir, 'mu.', suffl(nsuff)
c$$$      open(22, file = filnam, err = 9980, iostat = ioval)
c$$$      write(iuout, *) 
c$$$     $     'writing total system dipole data to mu.', suffl(nsuff)
      write(filnam, *) datdir, 'adebug.', suffl(nsuff)
      open(23, file = filnam, err = 9980, iostat = ioval)
      write(iuout, *) 'writing miscellaneous stuff to adebug.', 
     $     suffl(nsuff)
      write(filnam, *) datdir, 'p.', suffl(nsuff)
      open(24, file = filnam, err = 9980, iostat = ioval)
      write(iuout, *) 'writing pressure data to p.', suffl(nsuff)
      if (pbflag) then
         write(filnam, *) datdir, 'mu.', suffl(nsuff)
         open(25, file = filnam, err = 9980, iostat = ioval)
         write(iuout, *) 'writing system dipole data to mu.',
     $        suffl(nsuff)
      endif
c***********************************************************************
c     Read in the first two records from the first data file, to get
c     the distance between timesteps:
c***********************************************************************
      write(filnam, *) fildir(1:index(fildir, ' ')-1), 'dyn.', suffl(1)
      open(iutmp, file = filnam, status = 'old',
     $     form = 'unformatted', err = 9980, iostat = ioval)
      read(iutmp) iovers
      ioval = iurdmp(iutmp, rio, idbyte / nspbyt)
      if (ioval .eq. -1) then
         go to 9991
      else if (ioval .ne. 0) then
         go to 9990
      endif
      itrone = rio(8 * natoms + 5)
      ioval = iurdmp(iutmp, rio, idbyte / nspbyt)
      if (ioval .eq. -1) then
         go to 9991
      else if (ioval .ne. 0) then
         go to 9990
      endif
      itemp2 = rio(8 * natoms + 5)
      itrgap = itemp2 - itrone
      if (itrgap .le. 0) then
         write(iuout, *) 'analyz: 1st two records are weird...'
         stop
      endif
      close(iutmp)
c***********************************************************************
c     loop through the data files:
c***********************************************************************
      nrec = 0
      i1rec = 1
      iwrite = 0
      do 1200 idatfl = 1, nsuff
         write(filnam,*) fildir(1:index(fildir, ' ')-1),
     $        'dyn.', suffl(idatfl)
         write(iuout,*) 'reading from dyn.',suffl(idatfl)
         open(iudyn, file = filnam, status = 'old',
     $        form = 'unformatted', err = 9980, iostat = ioval)
         read(iudyn) iovers
         if (iovers .ne. iodynv) then
            write(iuout, *) 'analyz: incompatible dyn file'
            stop
         endif
c***********************************************************************
c     read in the first record from the next data file, to see where it
c     starts, and avoid overlap:
c***********************************************************************
         if (idatfl .ne. nsuff) then
            write(filnam,*) fildir(1:index(fildir, ' ')-1), 
     $           'dyn.', suffl(idatfl+1)
            open(iutmp, file = filnam, status = 'old',
     $           form = 'unformatted', err = 9980, iostat = ioval)
            read(iutmp) iovers
            ioval = iurdmp(iutmp, rio, idbyte / nspbyt)
            if (ioval .eq. -1) then
               go to 9991
            else if (ioval .ne. 0) then
               go to 9990
            endif
            istop = rio(8 * natoms + 5)
            close(iutmp)
         else
            istop = ibigp
         endif
c$$$      itrone = rio(8 * natoms + 5)
c$$$      ioval = iurdmp(iudyn, rio, idbyte / nspbyt)
c$$$      if (ioval .eq. -1) then
c$$$         go to 1100
c$$$      else if (ioval .ne. 0) then
c$$$         go to 9990
c$$$      endif
c$$$      itemp2 = rio(8 * natoms + 5)
c$$$      itrgap = itemp2 - itrone
c$$$      if (itrgap .le. 0) stop 
c$$$     $     'analyz: 1st two records are for the same iteration'
c$$$      rewind(iudyn)
c$$$      read(iudyn) iovers
c$$$      nrec = 0
c$$$      i1rec = 1
c***********************************************************************
c     read records from this data file until the end, and do stuff with
c     the numbers:
c***********************************************************************
 720     continue
         trnord = 0.0d0
         nrec = nrec + 1
c$$$      if (nrec .ge. maxrec + 1) then
c$$$         write(iuout, *) 'only ', maxrec, ' records allowed'
c$$$         stop
c$$$      endif
         if (niters .eq. (niters / iterio) * iterio .and. niters .gt. 0)
     $        write(iuout, *) 'Timestep ', niters
         ioval = iurdmp(iudyn, rio, idbyte / nspbyt)
         if (ioval .eq. -1) then
            go to 1100
         else if (ioval .ne. 0) then
            go to 9990
         endif
         do 730 iatom = 1,natoms
            temp = rio(iatom)
            pos(iatom,1) = temp
c$$$         posxt(iatom,1,nrec) = temp
 730     continue
         do 735 iatom = 1,natoms
            temp = rio(natoms + iatom)
            pos(iatom,2) = temp
c$$$         posxt(iatom,2,nrec) = temp
 735     continue
         do 740 iatom = 1,natoms
            temp = rio(2 * natoms + iatom)
            pos(iatom,3) = temp
c$$$         posxt(iatom,3,nrec) = temp
 740     continue
         do 745 iatom = 1,natoms
            vel(iatom,1) = rio(3 * natoms + iatom)
 745     continue
         do 750 iatom = 1,natoms
            vel(iatom,2) = rio(4 * natoms + iatom)
 750     continue
         do 755 iatom = 1,natoms
            vel(iatom,3) = rio(5 * natoms + iatom)
 755     continue
         do 760 iatom = 1,natoms
            q(iatom) = rio(6 * natoms + iatom)
 760     continue
         do 765 iatom = 1,natoms
            qvel(iatom) = rio(7 * natoms + iatom)
 765     continue
         qke = rio(8 * natoms + 1)
         qpe = rio(8 * natoms + 2)
         ake = rio(8 * natoms + 3)
         ape = rio(8 * natoms + 4)
         niters = rio(8 * natoms + 5)
         rscale = rio(8 * natoms + 6)
         pressr = rio(8 * natoms + 7)
         dke = rio(8 * natoms + 8)
c***********************************************************************
c     Write out the charges and spring positions:
c***********************************************************************
         call brpair(2,1)
         xclo = -dxmat(1,1)
         yclo = -dymat(1,1)
         zclo = -dzmat(1,1)
         rclo = rmat(1,1) 
         rclh1 = rmat(2,1)
         rclh2 = rmat(3,1)
         call brmono(2)
         xom = dxmat(4,1)
         yom = dymat(4,1)
         zom = dzmat(4,1)
         rom = rmat(4,1)
         cosdip = cos((xclo * xom + yclo * yom + zclo * zom) / 
     $        rclo / rom)
         teste = q(4) / rclh1 / rclh1 + q(5) / rclh2 / rclh2 + 
     $        q(6) / rclo / rclo
         call brmono(1)
         write(23, *) niters * dt, dxmat(1,2), dymat(1,2), dzmat(1,2),
     $        rmat(1,2), rclo, rclh1, rclh2, q(4) / 100, q(5) / 100, 
     $        q(6) / 100, cosdip, teste / 10
c$$$     $        (q(iatmol(163,i)), i = 2, 4),
c$$$     $        (q(iatmol(29,i)), i = 2, 4), (q(iatmol(109,i)), i = 2, 4),
c$$$     $        (q(iatmol(10,i)), i = 2, 4), (q(iatmol(80,i)), i = 2, 4),
c$$$     $        (q(iatmol(78,i)), i = 2, 4), (q(iatmol(17,i)), i = 2, 4),
c$$$     $        (q(iatmol(28,i)), i = 2, 4), (q(iatmol(125,i)), i = 2, 4),
c$$$     $        (q(iatmol(3,i)), i = 2, 4)
c***********************************************************************
c     If niters exceeds istop, then this run crashed, and the next
c     data file should be used, starting with this timestep:
c************************************************** *********************
         if (niters .ge. istop) then
            nrec = nrec - 1
            write(iuout, *) suffl(idatfl)(1:index(suffl(idatfl),' ')),
     $           'crashed near t =', istop * dt
            go to 1100
         endif
c***********************************************************************
c     Do some pressure averaging:
c***********************************************************************
         if (rscale .gt. 0.0d0) then
            npr = 0
            psum = 0.0d0
            pave = 0.0d0
            pssq = 0.0d0
         else
            if (pressr .ne. 0.0) then
               npr = npr + 1
               psum = psum + pressr
               pave = psum / npr
               pssq = pssq + pressr * pressr
            endif
         endif
c***********************************************************************
c     Check the data spacing and T scaling:
c***********************************************************************
         nrec2 = (niters - itrone) / itrgap + 1
         if (lspcok .and. nrec2 .ne. nrec) then
            lspcok = .false.
            write(iuout, *) 'analyz: WARNING! irregular data spacing'
            write(iuout, *) 'dt = ', nrec, ' should be ', nrec2
         endif
         if (rscale .gt. 0.0d0) then
            i1rec = nrec
         endif
c***********************************************************************
c     Calculate the translation order parameter:
c***********************************************************************
         do 1000 imol = 1,nmol
            dotprd = pos(iatmol(imol, 1),1) * rltvec(1) 
     $           + pos(iatmol(imol, 1),2) * rltvec(2) 
     $           + pos(iatmol(imol, 1),3) * rltvec(3)
            trnord = trnord + dcos(dotprd)
 1000    continue
         trnord = trnord / nmol
c***********************************************************************
c     Calculate the temperature:
c***********************************************************************
         rketmp = gettmp()
c***********************************************************************
c     Calculate some residence time info:
c***********************************************************************
         if (coorfl) then
            do 1020 imol1 = 1, nmol
               iatom1 = iatmol(imol1,1)
               do 1015 imol2 = 1, imol - 1
                  iatom2 = iatmol(imol2,1)
                  dx = pos(iatom2,1) - pos(iatom1,1)
                  dy = pos(iatom2,2) - pos(iatom1,2)
                  dz = pos(iatom2,3) - pos(iatom1,3)
                  if (pbflag) then
                     dx = dx - dnint(dx / boxsiz(1)) * boxsiz(1)
                     dy = dy - dnint(dy / boxsiz(2)) * boxsiz(2)
                     dz = dz - dnint(dz / boxsiz(3)) * boxsiz(3)
                  endif
                  r2 = dx * dx + dy * dy + dz * dz
                  r = sqrt(r2)
                  if (r .lt. rcoord) then
                     if (.not. nbrfl(imol1,imol2)) then
                        nbrfl(imol1,imol2) = .true.
                        if (.not. xnbrfl(imol1,imol2)) then
                           nbrsnc(imol1,imol2) = niters
                        endif
                     endif
                  else
                     if (nbrfl(imol1,imol2)) then
                        nbrfl(imol1,imol2) = .false.
                        xnbrfl(imol1,imol2) = .true.
                        nbrout(imol1,imol2) = niters
                     else if (xnbrfl(imol1,imol2)) then
                        if (niters - nbrout(imol1,imol2) .gt. nout) then
                           xnbrfl(imol1,imol2) = .false.
                           tres = (nbrout(imol1,imol2) - 
     $                          nbrsnc(imol1,imol2)) * dt
                           ibin = tres / tcmax * ncbin + 1
                           if (ibin .lt. ncbin .and. ibin .gt. 0) then
                              ntres(ibin) = ntres(ibin) + 1
                           endif
                           nntres = nntres + 1
                        endif
                     endif
                  endif
 1015          continue
 1020       continue
         endif
c$$$c***********************************************************************
c$$$c     Calculate some dipole stuff:
c$$$c***********************************************************************
c$$$      do 1065 imol = 1, nmol
c$$$         imty = molty(imol)
c$$$         do 1060 iqind = 1, nqmol(imty)
c$$$            iatom = iatmol(imol,iqmol(imty,iqind))
c$$$            csmu(nrec,1) = csmu(nrec,1) + pos(iatom,1) * q(iatom)
c$$$            csmu(nrec,2) = csmu(nrec,2) + pos(iatom,2) * q(iatom)
c$$$            csmu(nrec,3) = csmu(nrec,3) + pos(iatom,3) * q(iatom)
c$$$ 1060    continue
c$$$ 1065 continue
c$$$      smusq = csmu(nrec,1) * csmu(nrec,1) +
c$$$     $     csmu(nrec,2) * csmu(nrec,2) +
c$$$     $     csmu(nrec,3) * csmu(nrec,3)
c***********************************************************************
c     Now write out the dynamics stats (sequential stuff only if it
c     has been nsqint timesteps since the last write)
c***********************************************************************
         if (niters - iwrite .ge. nsqint) then
            iwrite = niters
            if (qdflag) then
               write(12, *) niters * dt, qke * fJm2kc / nmol,
     $              qpe * fJm2kc / nmol, (qke + qpe) * fJm2kc / nmol
            endif
            write(11, *) niters * dt, qke * qtfac, rketmp,
     $           dke * 2. / (3 * ndatom * k)
c$$$     $           ake * rtfac
            write(20, *) niters * dt, ake * fJm2kc / nmol,
     $           ape * fJm2kc / nmol, (ake + ape) * fJm2kc / nmol
            write(21, *) niters * dt, trnord
c$$$      write(22, *) niters * dt, csmu(nrec,1) * feA2D,
c$$$     $     csmu(nrec,2) * feA2D, csmu(nrec,3) * feA2D
c$$$c$$$  $        csmu(nrec,2) * feA2D, csmu(nrec,3) * feA2D, 
c$$$c$$$  $        smum(1) * feA2D, smum(2) * feA2D, smum(3) * feA2D, 
c$$$c$$$  $        smusq * feA2D * feA2D, smumsq * feA2D * feA2D, 
c$$$c$$$  $        1.0d0 + 3.0d0 * (gkfac / rktavg * smumsq)
c$$$         write(23, *) niters * dt, (q(i), i = 1, 20)
            write(24, *) niters * dt, pressr * fpmd2a, pave * fpmd2a, 
     $           sqrt(pssq / npr - pave * pave) * fpmd2a
         endif
c***********************************************************************
c     go back to read another record:
c***********************************************************************
         go to 720
c***********************************************************************
c     show up here at end of the current data file:
c***********************************************************************
 1100    continue
         close(iudyn)
c***********************************************************************
c     read from the system dipole file, and fiddle with those numbers:
c***********************************************************************
         if (pbflag) then
            write(filnam, *) datdir, 'sysdip.', suffl(idatfl)
            write(iuout, *) 'reading from sysdip.', suffl(idatfl)
            open(iudip, file = filnam, status = 'old',
     $           form = 'unformatted', err = 9980, iostat = ioval)
            read(iudip) iovers
            if (iovers .ne. iodipv) then
               write(iuout, *) 'analyz:  incompatible sysdip file'
               stop
            endif
c***********************************************************************
c     read in the first record from the next sysdip file, to see where
c     it starts, and avoid overlap:
c***********************************************************************
            if (idatfl .ne. nsuff) then
               write(filnam, *) datdir, 'sysdip.', suffl(idatfl+1)
               open(iutmp, file = filnam, status = 'old',
     $              form = 'unformatted', err = 9980, iostat = ioval)
               read(iutmp) iovers
               ioval = iurdmp(iutmp, rio, 4)
               if (ioval .eq. -1) then
                  go to 9991
               else if (ioval .ne. 0) then
                  go to 9990
               endif
               istop = int(rio(1))
               close(iutmp)
            else
               rstop = bigpos
            endif
 1110       continue
            ioval = iurdmp(iudip, rio, 4)
            if (ioval .eq. -1) then
               go to 1190
            else if (ioval .ne. 0) then
               go to 9990
            endif
c***********************************************************************
c     if rio(1) exceeds istop, then this run crashed, and the next data
c     file should be used, starting with this timestep:
c***********************************************************************
            itime = int(rio(1))
            do 1120 idim = 1, 3
               smu(idim) = rio(idim+1)
 1120       continue
            if (itime .ge. istop) then
               write(iuout, *) 
     $              suffl(idatfl)(1:index(suffl(idatfl), ' ')),
     $              'crashed near t =', istop * dt
               go to 1190
            endif
c***********************************************************************
c     do some system dipole averaging:
c***********************************************************************
            if (itime .ge. nscbas) then
               nmu = nmu + 1
               do 1130 idim = 1, 3
                  smusum(idim) = smusum(idim) + smu(idim)
                  smuave(idim) = smusum(idim) / nmu
                  smussq(idim) = smussq(idim) + smu(idim) * smu(idim)
 1130          continue
            endif
            if (nmuint * (itime / nmuint) .eq. itime) then
               write(25, '(i8, 6f8.3, 3f12.3)') 
     $              int(rio(1) * dt), (smu(i) * feA2D, i = 1, 3), 
     $              (smuave(i) * feA2D, i = 1, 3), 
     $              (smussq(i) / nmu, i = 1, 3)
            endif
            go to 1110
c***********************************************************************
c     show up here at the end of the current sysdip file:
c***********************************************************************
 1190       continue
c$$$            close(iudyn)
            close(iudip)
         endif
         nrec = nrec - 1
 1200 continue
c***********************************************************************
c     all done reading the dyn and sysdip files, so close up the
c     associated output files:
c***********************************************************************
      if (qdflag) then
         close(12)
      endif
      close(11)
      close(20)
      close(21)
      close(22)
      close(23)
      close(24)
      if (pbflag) then
         close(25)
      endif
c***********************************************************************
c     write out the residence time curve:
c***********************************************************************
      if (coorfl) then
         write(filnam, *) datdir, 'tres.', suffl(nsuff)
         open(iutmp, file = filnam, err = 9980, iostat = ioval)
         do 1205 ibin = 1, ncbin
            write(7, *) (ibin - 0.5) * tcmax / ncbin, ntres(ibin),
     $           dble(ntres(ibin)) / nntres
 1205    continue
         write(iuout, *) nntres, 'total neighbor exchanges'
         close(7)
      endif
c$$$      nrec = nrec - 1
c$$$c***********************************************************************
c$$$c     Calculate the system dipole time correlation function:
c$$$c***********************************************************************
c$$$      do 1130 irec1 = 1, nrec - 1
c$$$         do 1120 irec2 = 0, nrec - irec1
c$$$            irec2p = irec2 + 1
c$$$            temp = csmu(irec1,1) * csmu(irec1 + irec2,1) +
c$$$     $           csmu(irec1,2) * csmu(irec1 + irec2,2) +
c$$$     $           csmu(irec1,3) * csmu(irec1 + irec2,3)
c$$$            phi(irec2p) = phi(irec2p) + temp
c$$$            phissq(irec2p) = phissq(irec2p) + temp * temp
c$$$ 1120    continue
c$$$ 1130 continue
c$$$      write(filnam, *) datdir, 'phi.', suffix
c$$$      open(iutmp, file = filnam, err = 9980, iostat = ioval)
c$$$      write(iuout, *) 'writing system dipole time correlation ',
c$$$     $     'function data to phi.', suffix
c$$$      do 1140 irec = 0, nrec - 2
c$$$         irecp1 = irec + 1
c$$$         itemp = nrec - irec
c$$$         phi(irecp1) = phi(irecp1) / itemp
c$$$         phissq(irecp1) = phissq(irecp1) / itemp
c$$$         temp = tconf(0.99d0, itemp - 1) * 
c$$$     $        sqrt((phissq(irecp1) - 
c$$$     $        phi(irecp1) * phi(irecp1)) / itemp)
c$$$         write(iutmp, *) irec * itrgap * dt, 
c$$$     $        phi(irecp1) / phi(1), temp / phi(1)
c$$$ 1140 continue
c$$$      close(iutmp)
c$$$c***********************************************************************
c$$$c     Calculate the mean square diffusion for all dts, along with the
c$$$c     standard deviations of those means:
c$$$c***********************************************************************
c$$$      write(iuout, *) 'calculating rms distance stats...'
c$$$      do 1145 irec = i1rec, nrec - 1
c$$$         do 1143 irec2 = irec + 1, nrec
c$$$            ndt = irec2 - irec
c$$$            numtc(ndt) = numtc(ndt) + 1
c$$$            do 1142 imol = 1, nmol
c$$$               iatom = iatmol(imol,1)
c$$$               dx = posxt(iatom,1,irec2) - posxt(iatom,1,irec)
c$$$               dy = posxt(iatom,2,irec2) - posxt(iatom,2,irec)
c$$$               dz = posxt(iatom,3,irec2) - posxt(iatom,3,irec)
c$$$               rmsd = dx * dx + dy * dy + dz * dz
c$$$               rms(ndt) = rms(ndt) + rmsd
c$$$               rmsstd(ndt) = rmsstd(ndt) + rmsd * rmsd
c$$$ 1142       continue
c$$$ 1143    continue
c$$$ 1145 continue
c$$$      do 1146 idt = 1, nrec - i1rec
c$$$         npts = numtc(idt) * nmol
c$$$         rms(idt) = rms(idt) / npts
c$$$         rmsstd(idt) = sqrt(rmsstd(idt) / npts - rms(idt) * rms(idt))
c$$$ 1146 continue
c$$$c***********************************************************************
c$$$c     Calculate the diffusion constant by doing a weighted linear
c$$$c     regression on all of the mean squared diffusion data after a
c$$$c     certain point where they are assumed to be linear:
c$$$c***********************************************************************
c$$$      linit = linst / (itrgap * dt) + 1
c$$$      nlinit = nrec - i1rec - linit + 1
c$$$      if (nlinit .gt. 1) then
c$$$         mhavwt = 1
c$$$         call fit(rint(linit), rms(linit), nlinit, rmsstd(linit),
c$$$     $        mhavwt, dintcp, dslope, sigdin, sigdsl, sumdev, qscore)
c$$$         conv = fA2cm * fA2cm * fs2smd / (6.0d0 * itrgap * dt)
c$$$         dslope = dslope * conv
c$$$         sigdsl = sigdsl * conv
c$$$         write(iuout, *) 'D = ', dslope, '+/-', sigdsl
c$$$         write(iuout, *) '   intcpt = ', dintcp, '+/-', sigdin
c$$$         write(iuout, *) '   chi^2 = ', sumdev
c$$$         write(iuout, *) '   q score = ', qscore
c$$$      else
c$$$         write(iuout,*) 'Run was too short to calculate a ',
c$$$     $        'diffusion constant'
c$$$      endif
c$$$      write(filnam, *) datdir,'rmsdis.',suffix
c$$$      open(iutmp, file = filnam, err = 9980, iostat = ioval)
c$$$      write(iuout, *) 'writing equilibration data to rmsdis.', suffix
c$$$      do 1150 irec = 1, nrec - 1
c$$$         if (numtc(irec) .ne. 0) then
c$$$            write(iutmp, *) 
c$$$     $           irec * itrgap * dt, rms(irec), rmsstd(irec)
c$$$         endif
c$$$ 1150 continue
c$$$      close(iutmp)
c$$$c***********************************************************************
c$$$c     read in and spit out the system dipole stuff:
c$$$c***********************************************************************
c$$$      write(filnam, *) datdir, 'sysdip.', suffix
c$$$      write(iuout, *) 'reading system dipole info from sysdip.', 
c$$$     $     suffix
c$$$      open(iudip, file = filnam, status = 'old',
c$$$     $     form = 'unformatted', err = 9980, iostat = ioval)
c$$$      read(iudip) iovers
c$$$      if (iovers .ne. iodipv) then
c$$$         write(iuout, *) 'analyz:  incompatible sysdip file'
c$$$         stop
c$$$      endif
c$$$      write(filnam, *) datdir, 'fsysdip.', suffix
c$$$      open(iutmp, file = filnam, err = 9980, iostat = ioval)
c$$$      write(iuout, *) 'writing system dipole data to fsysdip.',
c$$$     $     suffix
c$$$ 1160 continue
c$$$      ioval = iurdmp(iudip, rio, 4)
c$$$      if (ioval .eq. -1) then
c$$$         go to 1180
c$$$      else if (ioval .ne. 0) then
c$$$         go to 9990
c$$$      endif
c$$$      write(iutmp, *) rio(1) * dt, rio(2) * feA2D, rio(3) * feA2D,
c$$$     $     rio(4) * feA2D
c$$$      go to 1160
c$$$ 1180 continue
c$$$      close(iutmp)
c***********************************************************************
c     Now read in all of the binned info from only the last suffix, and
c     spit it back out in more readable form:
c***********************************************************************
      write(filnam, *) datdir, 'bins.', suffl(nsuff)
      open(iudyn, file = filnam, form = 'unformatted',
     $     status = 'old', err = 9980, iostat = ioval)
      write(iuout, *) 'reading bin info from bins.', suffl(nsuff)
      read(iudyn, err = 9990, end = 9991, iostat = ioval) iovers
      if (iovers .ne. iobinv) then
         write(iuout, *) 'analyz:  incompatible bins file'
         stop
      endif
      read(iudyn, err = 9990, end = 9991, iostat = ioval) nstat
      icount = 0
 1300 continue
      sum = 0
      write(filnam, *) datdir, 'bins.', suffl(nsuff)
      read(iudyn, err = 9990, end = 9991, iostat = ioval)
     $     name, nbins, binmin, binmax, accflg, rsqflg
      if (nbins .gt. maxdat) stop 'analyz: maxdat is too small'
      binw = (binmax - binmin) / nbins
      if (dbflag .and. name .eq. 'rprE') then
         ioval = iurdmp(iudyn, rio, nbins * mxrdiv)
      else
         ioval = iurdmp(iudyn, rio, nbins)
      endif
      if (ioval .eq. -1) then
         write(filnam, *) datdir, 'bins.', suffl(nsuff)
         go to 9991
      else if (ioval .ne. 0) then
         write(filnam, *) datdir, 'bins.', suffl(nsuff)
         go to 9990
      endif
 1210 continue
      write(filnam, *) datdir, name(1:index(name, ' ')-1), 
     $     '.', suffl(nsuff)
      open(iutmp, file = filnam, err = 9980, iostat = ioval)
      write(iuout, *) 'writing binned stats to ', 
     $     name(1:index(name, ' ')-1), '.', suffl(nsuff)
      if (dbflag .and. name .eq. 'rprE') then
         rbw = hbox / mxrdiv
         do 1212 irdiv = 1, mxrdiv
            do 1211 ibin = 1, nbins
               write(iutmp, *) binmin + (ibin - 0.5) * binw, 
     $              (irdiv - 0.5) * rbw, 
     $              rio((irdiv - 1) * maxbin + ibin)
               sum = sum + rio((irdiv - 1) * maxbin + ibin) * 
     $              binw * rbw
 1211       continue
            write(iutmp, *) 
 1212    continue
c***********************************************************************
c     Some stats are binned by r^2, and some need to be accumulated as
c     they are written out:
c***********************************************************************
      else
         if (name(1:1) .eq. 'g' .and. name(2:2) .ne. 'k' .and.
     $        name(2:2) .ne. 'q') then
c***********************************************************************
c     Read in the variances on the g(r) data:
c***********************************************************************
            ioval = iurdmp(iudyn, riovar, nbins)
            if (ioval .eq. -1) then
               write(filnam, *) datdir, 'bins.', suffl(nsuff)
               go to 9991
            else if (ioval .ne. 0) then
               write(filnam, *) datdir, 'bins.', suffl(nsuff)
               go to 9990
            endif
c$$$c***********************************************************************
c$$$c     LS fit parabolae to windows of the g(r) and determine where the
c$$$c     peaks and valleys are:
c$$$c***********************************************************************
c$$$            write(filnam, *) datdir, 'LSfitdebug',
c$$$     $           name(1:index(name, ' ')-1), '.', suffl(nsuff)
c$$$            open(iutmp2, file = filnam, err = 9980, iostat = ioval)
c$$$            write(filnam, *) datdir, 'LSfit', 
c$$$     $           name(1:index(name, ' ')-1), '.', suffl(nsuff)
c$$$            open(iutmp3, file = filnam, err = 9980, iostat = ioval)
c$$$c$$$c***********************************************************************
c$$$c$$$c     First check for maxima and minima, smoothing until there are fewer
c$$$c$$$c     than one per 2 A.
c$$$c$$$c***********************************************************************
c$$$c$$$               do 1215 npts = 2, nbins
c$$$c$$$                  call smooth(rio, smrio, nbins, npts)
c$$$c$$$                  nex = xcount(smrio, nbins + 1 - npts)
c$$$c$$$                  write(iuout, *) npts, 'point binning found', nex, 
c$$$c$$$     $                 'extrema'
c$$$c$$$ 1215          continue
c$$$            write(iutmp2, *) 'Starting LS fit'
c$$$            do 1216 ibin = 1, nbins
c$$$               xvalvc(ibin) = binmin + (ibin - 0.5d0) * binw
c$$$               onevec(ibin) = 1.0d0 / riovar(ibin)
c$$$               xvec(ibin) = onevec(ibin) * xvalvc(ibin)
c$$$               x2vec(ibin) = xvec(ibin) * xvalvc(ibin)
c$$$               x3vec(ibin) = x2vec(ibin) * xvalvc(ibin)
c$$$               x4vec(ibin) = x3vec(ibin) * xvalvc(ibin)
c$$$               yvec(ibin) = rio(ibin) * onevec(ibin)
c$$$               xyvec(ibin) = yvec(ibin) * xvalvc(ibin)
c$$$               xxyvec(ibin) = xyvec(ibin) * xvalvc(ibin)
c$$$ 1216       continue
c$$$            nwinpt = int(0.5 / binw) + 1
c$$$            if (nwinpt .le. 4) then
c$$$               nwinpt = 4
c$$$               write(iuout, *) 'analyz: WARNING! 4 g(r) points ',
c$$$     $              'gives a window of ', 2 * binw, ' A'
c$$$            endif
c$$$            write(iutmp2, *) binw, 'A between data points'
c$$$            write(iutmp2, *) nwinpt, 'points per window'
c$$$            sum1 = 0.0d0
c$$$            sumx = 0.0d0
c$$$            sumx2 = 0.0d0
c$$$            sumx3 = 0.0d0
c$$$            sumx4 = 0.0d0
c$$$            sumy = 0.0d0
c$$$            sumxy = 0.0d0
c$$$            sumx2y = 0.0d0
c$$$            iibin = 0
c$$$ 1217       iibin = iibin + 1
c$$$            if (iibin .gt. nbins - nwinpt + 1) then 
c$$$               write(iutmp2, *) 'no data to fit'
c$$$               go to 1223
c$$$            endif
c$$$            if (rio(iibin) .eq. 0.0d0) go to 1217
c$$$            do 1218 ibin = iibin, iibin + nwinpt - 1
c$$$               sum1 = sum1 + onevec(ibin)
c$$$               sumx = sumx + xvec(ibin)
c$$$               sumx2 = sumx2 + x2vec(ibin)
c$$$               sumx3 = sumx3 + x3vec(ibin)
c$$$               sumx4 = sumx4 + x4vec(ibin)
c$$$               sumy = sumy + yvec(ibin)
c$$$               sumxy = sumxy + xyvec(ibin)
c$$$               sumx2y = sumx2y + xxyvec(ibin)
c$$$ 1218       continue
c$$$            xmin = xvalvc(iibin)
c$$$            xmax = xvalvc(iibin + nwinpt - 1)
c$$$            rmata = sumx2 * sumx4 - sumx3 * sumx3
c$$$            rmatd = sumx2 * sumx3 - sumx * sumx4
c$$$            temp = sumx2 * sumx2
c$$$            rmate = sum1 * sumx4 - temp
c$$$            rmatg = sumx * sumx3 - temp
c$$$            rmath = sumx * sumx2 - sum1 * sumx3
c$$$            rmati = sum1 * sumx2 - sumx * sumx
c$$$            denom = 1.0d0 / 
c$$$     $           (sum1 * rmata + sumx * rmatd + sumx2 * rmatg)
c$$$            c0 = (rmata * sumy + rmatd * sumxy + rmatg * sumx2y) *
c$$$     $           denom
c$$$            c1 = (rmatd * sumy + rmate * sumxy + rmath * sumx2y) *
c$$$     $           denom
c$$$            c2 = (rmatg * sumy + rmath * sumxy + rmati * sumx2y) *
c$$$     $           denom
c$$$            sig20 = rmata * denom
c$$$            sig21 = rmate * denom
c$$$            sig22 = rmati * denom
c$$$            sig10 = sqrt(sig20)
c$$$            sig11 = sqrt(sig21)
c$$$            sig12 = sqrt(sig22)
c$$$            rel20 = sig20 / (c0 * c0)
c$$$            rel21 = sig21 / (c1 * c1)
c$$$            rel22 = sig22 / (c2 * c2)
c$$$            write(iutmp2, *) 'window 1:'
c$$$            write(iutmp2, *) '  ', xmin, ' A to ', xmax, ' A'
c$$$            write(iutmp2, *) '  ', c0, ' + ', c1, ' x + ', c2, 'x^2'
c$$$            write(iutmp2, *) 
c$$$     $           '  ', sig10, '   ', sig11, '     ', sig12
c$$$            if (c2. ne. 0.0d0) then
c$$$               extrx = -0.5 * c1 / c2
c$$$               extrxe = sqrt(rel20 + rel22) * extrx
c$$$               temp = c0 * c2
c$$$               temp2 = 4.0d0 * temp - c1 * c1
c$$$               extry = temp2 * 0.25d0 / c2
c$$$               extrye = sqrt(((rel20 + rel22) * temp * temp + 
c$$$     $              2.0d0 * sig21) / (temp2 * temp2) + rel22) * extry
c$$$               write(iutmp2, *) '  extremum at ', extrx, extry
c$$$               write(iutmp2, *) '  errors of   ', extrxe, extrye
c$$$               write(iutmp2, *) '  strength of ', abs(c2)
c$$$               chisq = 0.0d0
c$$$               do 1219 jbin = 1, nwinpt
c$$$                  temp = rio(jbin) - c0 - c1 * xvalvc(jbin) - 
c$$$     $                 c2 * xvalvc(jbin) * xvalvc(jbin)
c$$$                  chisq = chisq + temp * temp
c$$$ 1219          continue
c$$$               chisq = chisq / nwinpt
c$$$               if (extrx .lt. xmax .and. extrx .gt. xmin) then
c$$$                  write(iutmp3, *) extrx, extry, chisq
c$$$                  if (c2 .gt. 0.0d0) then
c$$$                     write(iuout, *) 'Found a min at ', extrx
c$$$                  else if (c2 .lt. 0.0d0) then
c$$$                     write(iuout, *) 'Found a max at ', extrx
c$$$                  endif
c$$$               endif
c$$$            endif
c$$$ 1220       do 1222 isbin = iibin + 1, nbins - nwinpt + 1
c$$$               ipbin = isbin - 1
c$$$               ifbin = isbin + nwinpt - 1
c$$$               if (rio(ifbin) .eq. 0.0d0) go to 1223
c$$$               sum1 = sum1 - onevec(ipbin)
c$$$               sumx = sumx - xvec(ipbin)
c$$$               sumx2 = sumx2 - x2vec(ipbin)
c$$$               sumx3 = sumx3 - x3vec(ipbin)
c$$$               sumx4 = sumx4 - x4vec(ipbin)
c$$$               sumy = sumy - yvec(ipbin)
c$$$               sumxy = sumxy - xyvec(ipbin)
c$$$               sumx2y = sumx2y - xxyvec(ipbin)
c$$$               sum1 = sum1 + onevec(ifbin)
c$$$               sumx = sumx + xvec(ifbin)
c$$$               sumx2 = sumx2 + x2vec(ifbin)
c$$$               sumx3 = sumx3 + x3vec(ifbin)
c$$$               sumx4 = sumx4 + x4vec(ifbin)
c$$$               sumy = sumy + yvec(ifbin)
c$$$               sumxy = sumxy + xyvec(ifbin)
c$$$               sumx2y = sumx2y + xxyvec(ifbin)
c$$$               xmin = xvalvc(isbin)
c$$$               xmax = xvalvc(ifbin)
c$$$               rmata = sumx2 * sumx4 - sumx3 * sumx3
c$$$               rmatd = sumx2 * sumx3 - sumx * sumx4
c$$$               temp = sumx2 * sumx2
c$$$               rmate = sum1 * sumx4 - temp
c$$$               rmatg = sumx * sumx3 - temp
c$$$               rmath = sumx * sumx2 - sum1 * sumx3
c$$$               rmati = sum1 * sumx2 - sumx * sumx
c$$$               denom = 1.0d0 / 
c$$$     $              (sum1 * rmata + sumx * rmatd + sumx2 * rmatg)
c$$$               c0 = (rmata * sumy + rmatd * sumxy + rmatg * sumx2y) *
c$$$     $              denom
c$$$               c1 = (rmatd * sumy + rmate * sumxy + rmath * sumx2y) *
c$$$     $              denom
c$$$               c2 = (rmatg * sumy + rmath * sumxy + rmati * sumx2y) *
c$$$     $              denom
c$$$               sig20 = rmata * denom
c$$$               sig21 = rmate * denom
c$$$               sig22 = rmati * denom
c$$$               sig10 = sqrt(sig20)
c$$$               sig11 = sqrt(sig21)
c$$$               sig12 = sqrt(sig22)
c$$$               rel20 = sig20 / (c0 * c0)
c$$$               rel21 = sig21 / (c1 * c1)
c$$$               rel22 = sig22 / (c2 * c2)
c$$$               write(iutmp2, *) 'window ', isbin, ':'
c$$$               write(iutmp2, *) '  ', xmin, ' A to ', xmax, ' A'
c$$$               write(iutmp2, *) 
c$$$     $              '  ', c0, ' + ', c1, ' x + ', c2, 'x^2'
c$$$               write(iutmp2, *) 
c$$$     $              '  ', sig10, '   ', sig11, '     ', sig12
c$$$               if (c2 .ne. 0.0d0) then
c$$$                  extrx = -0.5 * c1 / c2
c$$$                  extrxe = sqrt(rel20 + rel22) * extrx
c$$$                  temp = c0 * c2
c$$$                  temp2 = 4.0d0 * temp - c1 * c1
c$$$                  extry = temp2 * 0.25d0 / c2
c$$$                  extrye = sqrt(((rel20 + rel22) * temp * temp + 
c$$$     $                 2.0d0 * sig21) / (temp2 * temp2) + rel22) * 
c$$$     $                 extry
c$$$                  write(iutmp2, *) '  extremum at ', extrx, extry
c$$$                  write(iutmp2, *) '  errors of   ', extrxe, extrye
c$$$                  write(iutmp2, *) '  strength of ', abs(c2)
c$$$                  chisq = 0.0d0
c$$$                  do 1221 jbin = isbin, ifbin
c$$$                     temp = rio(jbin) - c0 - c1 * xvalvc(jbin) -
c$$$     $                    c2 * xvalvc(jbin) * xvalvc(jbin)
c$$$                     chisq = chisq + temp * temp
c$$$ 1221             continue
c$$$                  chisq = chisq / nwinpt
c$$$                  write(iutmp2, *) 'chisq = ', chisq
c$$$                  if (extrx .lt. xmax .and. extrx .gt. xmin) then
c$$$                     write(iutmp3, *) extrx, extry, chisq
c$$$                     if (c2 .gt. 0.0d0) then
c$$$                        write(iuout, *) 'Found a min at ', extrx
c$$$                     else if (c2 .lt. 0.0d0) then
c$$$                        write(iuout, *) 'Found a max at ', extrx
c$$$                     endif
c$$$                  endif
c$$$               endif
c$$$ 1222       continue
c$$$ 1223       continue
c$$$            close(iutmp2)
c$$$            close(iutmp3)
c***********************************************************************
c     Read in the Soper g(r) data:
c***********************************************************************
            if (iobwfl) then
               gdiff = 0.0d0
               numg = 0
               if (name(2:5) .eq. '0808') then
                  write(filnam, *) datdir, 'gOO_soper.dat'
                  write(iuout, *) 
     $                 'reading experimental data from gOO_soper.dat'
               else if (name(2:5) .eq. '0801' .or.
     $                 name(2:5) .eq. '0108') then
                  write(filnam, *) datdir, 'gOH_soper.dat'
                  write(iuout, *)
     $                 'reading experimental data from gOH_soper.dat'
               else if (name(2:5) .eq. '0101') then
                  write(filnam, *) datdir, 'gHH_soper.dat' 
                  write(iuout, *) 
     $                 'reading experimental data from gHH_soper.dat'
               else
                  stop 'analyz: unknown g(r)'
               endif
               open(iutmp2, file = filnam, err = 9980, iostat = ioval)
               nlines = 1
               read(iutmp2, *, err = 9990, end = 1240, iostat = ioval) 
     $              tmpwid, gexpt(1)
               tmpwid = tmpwid * 2
               nlines = 2
 1230          continue
               read(iutmp2, *, err = 9990, end = 1240, iostat = ioval) 
     $              rjunk, gexpt(nlines)
               nlines = nlines + 1
               go to 1230
 1240          continue
               nlines = nlines - 1
               close(iutmp2)
            endif
         endif
         do 1245 ibin = 1, nbins
            if (accflg) then
               stat = stat + rio(ibin)
            else
               stat = rio(ibin)
            endif
            xval = binmin + (ibin - 0.5d0) * binw
            sum = sum + rio(ibin) * binw
            if (rsqflg)
     $           xval = sqrt(xval)
            if (name(1:1) .eq. 'g' .and. name(2:2) .ne. 'k' .and.
     $           name(2:2) .ne. 'q') then
               if (iobwfl) then
                  ixpbin = int(xval / tmpwid) + 1
                  if (ixpbin .ge. 2 .and. ixpbin .le. nlines - 1) then
                     xpxval = (xval + tmpwid * (0.5d0 - ixpbin)) /
     $                    tmpwid
                     xpyval = 0.5d0 * 
     $                    (xpxval * (xpxval - 1) * gexpt(ixpbin - 1) +
     $                    2 * (1 - xpxval * xpxval) * gexpt(ixpbin) +
     $                    xpxval * (xpxval + 1) * gexpt(ixpbin + 1))
                     xpdiff = stat - xpyval
                     gdiff = gdiff + xpdiff * xpdiff * binw
                     numg = numg + 1
                  else
                     xpdiff = 0.0d0
                  endif
                  write(iutmp, *) 
     $                 xval, stat, sqrt(riovar(ibin)), xpdiff
               else
                  write(iutmp, *) xval, stat, sqrt(riovar(ibin))
               endif
            else
               write(iutmp, *) xval, stat
            endif
 1245    continue
         if (name(1:1) .eq. 'g' .and. name(2:2) .ne. 'k' .and. 
     $        name(2:2) .ne. 'q' .and. iobwfl) then
            write(iuout, '(f15.5,a,i3,a)') gdiff,
     $           ' difference from Soper, for ', numg, 
     $           ' data points'
         else
            write(iuout, '(a, f5.1, a)') '(', 100 * sum, '% caught)'
         endif
      endif
      close(iutmp)
      icount = icount + 1
      if (icount .lt. nstat) go to 1300
c$$$      else
c$$$         write(iuout,*) 'No real dynamics'
c$$$      endif
      apeave = apesum / nbinst
      qpeave = qpesum / nbinst
      pave = psum / nrec
      write(iuout, '(a,f9.3,a,f8.3)') ' <PE> = ', apeave,
     $     ' +/- ', sqrt(apessq / nbinst - apeave * apeave)
      write(iuout, '(a,f9.3,a,f8.3)') '<qPE> = ', qpeave,
     $     ' +/- ', sqrt(qpessq / nbinst - qpeave * qpeave)
      write(iuout, '()') '<P> = ', pave * fpmd2a,
     $     ' +/- ', sqrt(pssq / nrec - pave * pave) * fpmd2a
      write(iuout,*) 'Finished reading.'
      stop
c***********************************************************************
c     Various error messages:
c***********************************************************************
c$$$ 9970 write(iuout, *) 'analyz: atom identities in simdat.', suffix,
c$$$     $     ' and ', filnam, ' don''t match'
c$$$      stop
c$$$ 9971 write(iuout, *) 'analyz: molecule assignments in simdat.', suffix,
c$$$     $     ' and ', filnam, ' don''t match'
c$$$      stop
 9980 write(iuout, *) 'analyz: error opening file ',filnam
      call ioerr(ioval)
      stop
 9990 write(iuout, *) 'analyz: error reading file ',filnam
      call ioerr(ioval)
      stop
 9991 write(iuout, *) 'analyz: unexpected end-of-file in',filnam
      stop
c$$$ 9992 write(iuout, *) 'analyz: data from dyn.',filnam,
c$$$     $           ' doesn''t match data from anneal.',filnam
c$$$      stop
      end
