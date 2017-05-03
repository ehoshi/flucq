      program qdyn
c***********************************************************************
c     This program will do dynamics for a system of molecules, while
c     letting their individual partial charges vary.
c     The charges are calculated by equalizing electronegativities, much
c     like Rappe and Goddard.  This is done in conjunction with the 
c     dynamics in an extended Lagrangian type algorithm.
c***********************************************************************
c     Units:
c
c     mass:        1 atomic mass unit (amu)
c     time:        1 decifemtosecond (1 dfs = 1e-16 s)
c     length:      1 Angstrom (A)
c     charge:      1 electron (e)
c     temperature: 1 Kelvin (K)
c
c     Everything else (forces, energies) is derived in terms of these.
c
c     energy:      1 amu A^2 / dfs^2
c
c     The charge properties are a little weird, since their motions
c     occur in 1-D charge space, not in 3-D distance space:
c
c     qvel:        1 e / dfs
c     qmass:       1 amu A^2 / e^2
c***********************************************************************
c     RESPA has been hardwired in, with a far/near/far split on
c     distances nested inside a light/heavy/light split on masses.  The
c     loop counts (nl, nlf, nhf) can be set to recover velocity Verlet.
c***********************************************************************
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c     
      logical   tscflg, getnr, getfar, getlt, gethvy, gete
      real*4    rio(maxdat)
      real*8    qinit(maxatm)
c
      integer*4 rdsim
c     
c***********************************************************************
c     Initialize stuff:
c***********************************************************************
c     bdeflg is a hack to count solvation (bonding) energies
c***********************************************************************
      bdeflg = .true.
      neede = .true.
      tscflg = .false.
      call qdinit
      nqinit = 0
      nqscbs = 0
      sqreps = dsqrt(pi4eps)
c***********************************************************************
c     Read in the simulation details, and initialize all of the
c     associated stuff:
c***********************************************************************
      itemp = rdsim()
      write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), '/debug.', 
     $     suffix
      open(99, file = filnam, err = 9990, iostat = ioval)
      if (cbnflg) then
         first = .false.
      else
         first = .true.
      endif
      if (ewlkfl) call getkvc
      if (itemp .ne. 0) then
         call ioerr(itemp)
         stop
      endif
      rclbw = sqrt(3.d0) * boxby2(1) / nrclbn
      do 450 ibin = 1, nrclbn
         bdercl(ibin) = 0.d0
         bdrssq(ibin) = 0.d0
         nbdr(ibin) = 0
 450  continue
c      iterio = (ndstep - niters) / 20
      iterio = 10000
      if (iterio .lt. 1) iterio = 1
      Tinv = 1.0d0 / T
      if (dbflag) rkTinv = rkinv / T
      if (nfqaty .ne. 0) then
         qtfac = 2.d0 / (nqdof * k)
      endif
c$$$      rtfac = 2 / (nrdof * k)
      do 460 iatom = 1, natoms
         qinit(iatom) = q(iatom)
  460 continue
      if (ioflg) then
         call ioopen()
      endif
c$$$      if (dbflag) then
c$$$         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
c$$$     $        '/drudepos.', suffix
c$$$         open(21, file = filnam, err = 9990, iostat = ioval)
c$$$         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
c$$$     $        '/cmvel.', suffix
c$$$         open(22, file = filnam, err = 9990, iostat = ioval)
c$$$      endif
c***********************************************************************
c     Write out a preliminary version of simdat.SUF, that will allow
c     analyze to be run before qdyn finishes, in order to check on the
c     simulation as it is running:
c***********************************************************************
      call wrtsim()
c***********************************************************************
c     Dynamics starts here; initialize stuff needed each iteration:
c***********************************************************************
      write(iuout, *) 'Starting dynamics at step ', niters
c     write once at the start, so we can analyze properties even at the
c     early stages of a run
      call wrteor
 1000 continue
      call itinit
c***********************************************************************
c     The RESPA breakup is first into light/heavy/light, then both the 
c     heavy and light propagators are broken up into slow/fast/slow
c     forces.  
c***********************************************************************
c     Some calls
c     are commented out b/c the light particles currently need no charge
c     dynamics, and vverlv can be called instead of rattlv because the
c     light particles currently need no RATTLEing.
c***********************************************************************
      if (moving) then
         if (nlatom .ne. 0) then
            getnr = .true.
            gete = .false.
            getlt = .true.
            gethvy = .false.
            do 1100 idt = 1, nl / 2
               call vverlv(dtl, ifar, ilight)
               call qvvrlv(dtl, ifar, ilight)
               do 1090 iidt = 1, nlf
                  getfar = (iidt + 1 .gt. nlf)
                  gethvy = (idt + 1 .gt. nl / 2 .and. iidt + 1 .gt. nlf)
                  call vverlv(dtlf, inear, ilight)
                  call qvvrlv(dtlf, inear, ilight)
                  call vverlx(dtlf, ilight)
                  call qvvrlx(dtlf, ilight)
                  call getf(getnr, getfar, getlt, gethvy, gete)
                  call qvvrlv(dtlf, inear, ilight)
                  call vverlv(dtlf, inear, ilight)
 1090          continue
               call qvvrlv(dtl, ifar, ilight)
               call vverlv(dtl, ifar, ilight)
 1100       continue
         endif
         getnr = .true.
c$$$         getlt = .false.
         gethvy = .true.
         call rattlv(dths, ifar, iheavy)
         call qvvrlv(dths, ifar, iheavy)
         do 1110 idt = 1, nhf
            getlt = (idt .eq. nhf)
            getfar = (idt .eq. nhf)
            gete = (idt .eq. nhf)
            call qvvrlv(dthf, inear, iheavy)
            call rattlx(dthf, inear, iheavy)
            call qvvrlx(dthf, iheavy)
            call getf(getnr, getfar, getlt, gethvy, gete)
            call rattlv(dthf, inear, iheavy)
            call qvvrlv(dthf, inear, iheavy)
 1110    continue
         call qvvrlv(dths, ifar, iheavy)
         call rattlv(dths, ifar, iheavy)
         if (nlatom .ne. 0) then
            getnr = .true.
            getlt = .true.
            gethvy = .false.
            do 1120 idt = 1, nl / 2
               call vverlv(dtl, ifar, ilight)
c               call qvvrlv(dtl, ifar, ilight)
               do 1115 iidt = 1, nlf
                  getfar = (iidt + 1 .gt. nlf)
                  gete = (ioiter .and. idt + 1 .gt. nl / 2 .and.
     $                 iidt .eq. nlf)
c$$$                  gethvy = (idt + 1 .gt. nl / 2 .and. iidt + 1 .gt. nlf)
c$$$                  gete = (idt + 1 .gt. nl / 2 .and. iidt + 1 .gt. nlf)
                  call vverlv(dtlf, inear, ilight)
c                  call qvvrlv(dtlf, inear, ilight)
                  call vverlx(dtlf, ilight)
c                  call qvvrlx(dtlf, ilight)
                  call getf(getnr, getfar, getlt, gethvy, gete)
c                  call qvvrlv(dtlf, inear, ilight)
                  call vverlv(dtlf, inear, ilight)
 1115          continue
c               call qvvrlv(dtl, ifar, ilight)
               call vverlv(dtl, ifar, ilight)
 1120       continue
         endif
      else
         getnr = .true.
         getfar = .true.
         getlt = .true.
         gethvy = .true.
         gete = .true.
         movlt = .true.
         movhvy = .true.
         call getf(getnr, getfar, getlt, gethvy, gete)
      endif
      do 1530 imol = 1, nmol
         imty = molty(imol)
         do 1525 ifind = 1, nframe(imty)
            iatom = iatmol(imol,ifind)
            iaty = iatype(imty,ifind)
            ake = ake + atmass(iaty) * (vel(iatom,1) * vel(iatom,1) +
     $           vel(iatom,2) * vel(iatom,2) +
     $           vel(iatom,3) * vel(iatom,3))
            qke = qke + qmass * qvel(iatom) ** 2
 1525    continue
 1530 continue
      ake = ake * 0.5
c$$$      write(78, '(f9.1,8f13.5)') niters * dtbig, (q(i), i = 1,8)
c***********************************************************************
c     if it's a Drude resolving iteration, then resolve.
c***********************************************************************
      if (dsiter) then
         write(iuout, *) 'resolving for Drude oscillator positions ',
     $        'at timestep ', niters, ':'
         call dsolv
      endif
c***********************************************************************
c     If it's a charge resolving iteration, then resolve.  The charge 
c     velocities need to be reset here, and the qke recalculated:
c***********************************************************************
      if (qsiter) then
         call qsolv2
         write(iuout, *) 'resolved for charges at timestep ', 
     $        niters, ':'
         qke = 0.0d0
         rms = 0.d0
         do 1555 imol = 1, nmol
            imty = molty(imol)
            do 1550 iqind = 1, nfqmol(imty)
               iatom = iatmol(imol,ifqmol(imty,iqind))
               dq = q(iatom) - qans(iatom)
               rms = rms + dq * dq
               q(iatom) = qans(iatom)
               qvel(iatom) = 0.0d0
               qke = qke + qmass * qvel(iatom) * qvel(iatom)
 1550       continue
 1555    continue
         rms = sqrt(rms / nfqatm)
         write(iuout, *) '   rms charge deviation ', rms
         if (rms .gt. qrmsct .and. dslvfl) then
            qwrong = .true.
            if (dslvfl) then
               dwrong = .true.
            endif
            write(iuout, *) '  charges not converged...'
         else
            qwrong = .false.
            write(iuout, *) '  charges converged okay'
         endif
         if (first) then
            write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1),
     $           '/charges.', suffix
            open(iutmp, file = filnam, err = 9990, iostat = ioval)
            do 1560 iatom = 1, natoms
               write(iutmp,'(3i4,f12.8)') iatom, ident(iatom),
     $              molec(iatom), q(iatom)
 1560       continue
            close(iutmp)
            write(iuout,*) 'The charges are written to charges.',
     $           suffix
         endif
      endif
      qke = qke * 0.5
      ake = ake + qke
      if (ioiter) then
         qtote = qke + qpe(ilight) + qpe(imixed) + qpe(iheavy)
      endif
      qtemp = qke * qtfac
      if (ioiter) then
         atote = ake + ape(ilight) + ape(imixed) + ape(iheavy)
      endif
c***********************************************************************
c     Calculate the temperature, and keep stats:
c***********************************************************************
c$$$      rketmp = gettmp()
      if (ndatom .gt. 0) then
         call getdke()
      endif
      call gettmp()
      ntstep = ntstep + 1
      rktsum = rktsum + rketmp
      rktssq = rktssq + rketmp * rketmp
c***********************************************************************
c     quick hack to get stats on solvation energies -- it assumes that
c     a chloride ion (or some other important molecule) is molecule #1
c***********************************************************************
c$$$      if (bdeflg) then
c$$$         rjunk = 0.d0
c$$$         rjunk = rjunk + ebit(1,ilight) + ebit(1,imixed) 
c$$$     .        + ebit(1,iheavy)
c$$$         rjunk = rjunk 
c$$$     .        + ebit(2,ilight) + ebit(2,imixed) + ebit(2,iheavy)
c$$$         rjunk = rjunk 
c$$$     .        + ebit(3,ilight) + ebit(3,imixed) + ebit(3,iheavy)
c$$$         rjunk = rjunk
c$$$     .        + ebit(4,ilight) + ebit(4,imixed) + ebit(4,iheavy)
c$$$         rjunk = rjunk 
c$$$     .        + ebit(5,ilight) + ebit(5,imixed) + ebit(5,iheavy)
c$$$         rjunk = rjunk 
c$$$     .        + ebit(6,ilight) + ebit(6,imixed) + ebit(6,iheavy)
c$$$         rjunk = rjunk 
c$$$     .        + ebit(7,ilight) + ebit(7,imixed) + ebit(7,iheavy)
c$$$         rjunk = rjunk
c$$$     .        + ebit(8,ilight) + ebit(8,imixed) + ebit(8,iheavy)
c$$$         rjunk = rjunk 
c$$$     .        + ebit(9,ilight) + ebit(9,imixed) + ebit(9,iheavy)
c$$$         rjunk = rjunk
c$$$     .        + ebit(10,ilight) + ebit(10,imixed) + ebit(10,iheavy)
c$$$         write(iuout,*) rjunk
c$$$         rjunk = 0.d0
c$$$         do 6969 imol = 1, nmol
c$$$            rjunk = rjunk + coule(imol,ilight)
c$$$     .           + coule(imol,imixed) + coule(imol,iheavy)
c$$$            rjunk = rjunk + rlje(imol,ilight)
c$$$     .           + rlje(imol,imixed) + rlje(imol,iheavy)
c$$$            rjunk = rjunk + (selfe(imol,ilight) + selfe(imol,imixed)
c$$$     .           + selfe(imol,iheavy)) * 2.d0
c$$$ 6969    continue
c$$$         rjunk = rjunk * 0.5d0
c$$$         write(iuout, *) rjunk
c$$$      endif
      if (bdeflg .and. nbde * (niters / nbde) .eq. niters) then
         do 7988 imol = 1, nmol
            sdmre(imol) = 0.d0
            do 7987 jmol = 1, nmol
               sdmre(imol) = sdmre(imol) + dimere(imol,jmol,ilight)
     .              + dimere(imol,jmol,imixed)
     .              + dimere(imol,jmol,iheavy)
 7987       continue
c$$$            if (sdmre(imol) .eq. 0.d0) then
c$$$               write(iuout, *) niters, imol, sdmre(imol)
c$$$            endif
 7988    continue
         do 8990 imol = 1, nmol
            solve(imol) = 0.d0
            cemine = coule(imol,ilight) + coule(imol,imixed)
     .           + coule(imol,iheavy)
            solve(imol) = solve(imol) + cemine
            solve(imol) = solve(imol) + rlje(imol,ilight)
     .           + rlje(imol,imixed) + rlje(imol,iheavy)
            solve(imol) = solve(imol) + selfe(imol,ilight)
     .           + selfe(imol,imixed) + selfe(imol,iheavy)
c$$$            solve(imol) = solve(imol) + selfe(imol,ilight)
c$$$     .           + selfe(imol,imixed) + selfe(imol,iheavy)
            do 8988 jmol = 1, nmol
               djunk = dimere(imol,jmol,ilight) 
     .              + dimere(imol,jmol,imixed) 
     .              + dimere(imol,jmol,iheavy)
               solve(imol) = solve(imol) + djunk
     .              / sdmre(jmol) * (selfe(jmol,ilight)
     .              + selfe(jmol,imixed) + selfe(jmol,iheavy))
 8988       continue
            solve(imol) = solve(imol) / 2.d0
 8990    continue
         do 1590 imol = 2, nmol
            call brhead(imol, 1, dxinc, dyinc, dzinc)
            ibin = int(rmat(1,1) / rclbw) + 1
c           to prevent problems when rmat is NaN
            if (ibin .lt. 1) then
               ibin = 1
            else if (ibin .gt. nrclbn) then
               ibin = nrclbn
            endif
            bdercl(ibin) = bdercl(ibin) + solve(imol)
            bdrssq(ibin) = bdrssq(ibin) + solve(imol)**2
            nbdr(ibin) = nbdr(ibin) + 1
 1590    continue
      endif	
c***********************************************************************
c     Rescale the real velocities if necessary.  If the system has
c     equilibrated sufficiently, then unset the velocity-scaling flag.
c     All velocities are scaled; if just the cm velocities should be
c     scaled, then changes need to be made:
c***********************************************************************
      if (frzflg .and. moving) then
         do 1605 imol = 1, nmol
            imty = molty(imol)
            do 1600 ifind = 1, nframe(imty)
               iatom = iatmol(imol,ifind)
               vel(iatom,1) = 0.0d0
               vel(iatom,2) = 0.0d0
               vel(iatom,3) = 0.0d0
 1600       continue
 1605    continue
         ake = 0.0d0
         if (ioiter) then
            atote = ake + ape(ilight) + ape(imixed) + ape(iheavy)
         endif
      else if (vscflg .and. .not. eqflag .and. moving
     $        .and. nscint * (niters / nscint) .eq. niters) then
         if (ntstep .gt. 1) then
            rktavg = rktsum / ntstep
            write(iuout, *) 'rktavg = ', rktavg
            temp = dabs(rktavg - T)
            if (temp .gt. Ttol) then
               neffbn = ntstep * dtbig / stint
               write(iuout, *) 'neffbn = ', neffbn
               if (neffbn .gt. 1) then
                  rktvar = (rktssq / ntstep - rktavg ** 2) / 
     $                 dble(neffbn)
                  write(iuout, *) 'rktstd = ', sqrt(rktvar)
                  ttest = temp / sqrt(rktvar)
                  temp = tfunc(ttest, neffbn - 1)
                  write(iuout, *) 'F(t) = ', temp
c***********************************************************************
c     tfunc returns the two-sided t function, and we want the one-sided
c     version.  For there to be a 10% chance that the mean is as
c     far away as ttest std deviations in one specific direction 
c     (there's not much reason to be more strict than this), we need
c     F(ttest) > 0.8:
c***********************************************************************
                  if (temp .gt. 0.8) then
                     tscflg = .true.
                     write(iuout, *) 'Scaled temperature from ',
     $                    rktavg, ' K after ', niters - nscbas, ' steps'
c***********************************************************************
c     This is better for small clusters or other systems with large
c     variations in the temperature away from the average.  It attempts
c     to change the KE by an absolute amount, instead of by the fraction
c     by which the average is off, since for small systems the
c     instantaneous temp can differ substantially from the average.
c***********************************************************************
c     It would seem to make sense to change the KE by double the amount
c     needed, so that half can bleed into (or in from) the PE.  But I
c     tried this and it tends to overcorrect.
c***********************************************************************
                     write(iuout, *) 'rketmp = ', rketmp
                     temp = (rketmp + T - rktavg) / rketmp
                     if (temp .gt. 2.d0) then
                        temp = 2.d0
                        write(iuout, *) '(only doubling instantaneous ',
     $                       'temperature)'
                     else if (temp .lt. 0.0d0) then
                        temp = 0.0d0
                        write(iuout, *) '(can only freeze ',
     $                       'velocities)'
                     endif
                     write(iuout, *) 'KE scale factor = ', temp
                     temp = dsqrt(temp)
                     do 1615 imol = 1, nmol
                        imty = molty(imol)
                        do 1608 idrind = 1, ndrmol(imty)
                           ibsind = idrmol(imty,idrind,1)
                           ibatom = iatmol(imol,ibsind)
                           ispind = idrmol(imty,idrind,2)
                           isatom = iatmol(imol,ispind)
                           vel(isatom,1) = vel(isatom,1) + 
     $                          (temp - 1) * vel(ibatom,1)
                           vel(isatom,2) = vel(isatom,2) +
     $                          (temp - 1) * vel(ibatom,2)
                           vel(isatom,3) = vel(isatom,3) +
     $                          (temp - 1) * vel(ibatom,3)
                           vel(ibatom,1) = vel(ibatom,1) * temp
                           vel(ibatom,2) = vel(ibatom,2) * temp
                           vel(ibatom,3) = vel(ibatom,3) * temp
 1608                   continue
                        do 1610 ifind = 1, nframe(imty)
                           iatom = iatmol(imol,ifind)
                           iaty = iatype(imty,ifind)
                           if (isdrat(iaty)) then
                              go to 1610
                           endif
                           vel(iatom,1) = vel(iatom,1) * temp
                           vel(iatom,2) = vel(iatom,2) * temp
                           vel(iatom,3) = vel(iatom,3) * temp
 1610                   continue
 1615                continue
                     ake = ake * temp * temp
                     if (ioiter) then
                        atote = ake + ape(ilight) + ape(imixed) + 
     $                       ape(iheavy)
                     endif
                     rketmp = rketmp * temp * temp
                     write(iuout, *) 'new T = ', rketmp
                     nscbas = niters
                     rktsum = 0.d0
                     rktssq = 0.d0
                     ntstep = 0
                  endif
               endif
            endif
         endif
         if (niters - nscbas .ge. nsceql) then
            eqflag = .true.
            write(iuout, *) 'system has equilibrated at step ',
     $           niters, '(no recent T-scaling)'
         else if (niters .ge. mxsceq) then
            eqflag = .true.
            write(iuout, *) 'system has equilibrated at step ',
     $           niters, '(it''s been running long enough)'
         endif
      endif
      if (qvscfl .and. nscint * (niters / nscint) .eq. niters) then
         qtemp = qke * qtfac
         if (dabs((qtemp - qt) / qt) .gt. qtrang) then
            qwrong = .true.
            write(iuout, *) 'Charge temp is ', qtemp, 
     $           ' K at timestep ', niters
            write(iuout, *) 'Charges will be bottomed next timestep.'
         endif
      endif
c$$$      call getdke()
c$$$      if (dbflag) then
c$$$         call gtcomv(cmvel)
c$$$         write(22, *) niters * dtbig, (cmvel(i), i = 1, 3)
c$$$         do 1535 ietype = 1, 9
c$$$            write(iuout, *) ' ', ebit(ietype) / nmol * fJm2kc
c$$$ 1535    continue
c$$$         write(iuout, *) '--------------------------------'
c$$$      endif
c***********************************************************************
c     We're done an iteration, so dump out lots of stuff if necessary:
c***********************************************************************
c     This is a short "progress report" to show that we're still 
c     running:
c***********************************************************************
      if (iterio * (niters / iterio) .eq. niters .and.
     $     niters .ne. 0) then
         write(iuout, *) niters, ' steps:'
         if (moving) then
            if (.not. eqflag) then
               write(iuout, *) 'Currently equilibrating.'
            else
               write(iuout, *) 'Currently doing dynamics.'
            endif
         else
            write(iuout, *) 'Currently annealing.'
         endif
         write(iuout,*) '   qT = ',qtemp
         if (moving) then
            write(iuout, *) 'KE T = ', rketmp
         endif
         call wrteor
      endif
      if (dbflag) then
         do 1690 ii = 1, 10
            write(iuout, *) ii, ebit(ii,ilight) * fjm2kc / nmol, 
     $           ebit(ii,imixed) * fjm2kc / nmol, 
     $           ebit(ii,iheavy) * fjm2kc / nmol
 1690    continue
      endif
c***********************************************************************
c     This is the real I/O:
c***********************************************************************
      if (ioiter) then
         if (moving) then
            do 1720 iatom = 1,natoms
               rio(iatom) = pos(iatom,1)
 1720       continue
            do 1730 iatom = 1,natoms
               rio(natoms + iatom) = pos(iatom,2)
 1730       continue
            do 1740 iatom = 1,natoms
               rio(2 * natoms + iatom) = pos(iatom,3)
 1740       continue
            do 1750 iatom = 1,natoms
               rio(3 * natoms + iatom) = vel(iatom,1)
 1750       continue
            do 1760 iatom = 1,natoms
               rio(4 * natoms + iatom) = vel(iatom,2)
 1760       continue
            do 1770 iatom = 1,natoms
               rio(5 * natoms + iatom) = vel(iatom,3)
 1770       continue
            do 1780 iatom = 1,natoms
               rio(6 * natoms + iatom) = q(iatom)
 1780       continue
            do 1790 iatom = 1,natoms
               rio(7 * natoms + iatom) = qvel(iatom)
 1790       continue
            rio(8 * natoms + 1) = qke
            rio(8 * natoms + 2) = qpe(ilight) + qpe(imixed) +
     $           qpe(iheavy)
            rio(8 * natoms + 3) = ake
            rio(8 * natoms + 4) = ape(ilight) + ape(imixed) + 
     $           ape(iheavy)
            rio(8 * natoms + 5) = niters
            if (tscflg) then
               rio(8 * natoms + 6) = 1.
               tscflg = .false.
            else
               rio(8 * natoms + 6) = 0.
            endif
c$$$            wvir = wvirlj + wvirqt + wvirqa + wvirk + wvirkd
            wvir = wvirlj(1) + wvirlj(2) + wvirlj(3) +
     $           wvirqt(1) + wvirqt(2) + wvirqt(3) +
     $           wvirqa(1) + wvirqa(2) + wvirqa(3) +
     $           wvirk + wvirkd
c$$$            pressr = rnkbv * rketmp + wvir / 3. / vol
            pressr = rnkbv * rcmtmp + wvir / 3. / vol
c$$$            write(6, '(7(g10.4,1x))') rnkbv * rketmp * fpmd2a, 
c$$$     $           wvirlj / 3. / vol * fpmd2a, 
c$$$     $           wvirqt / 3. / vol * fpmd2a, 
c$$$     $           wvirqa / 3. / vol * fpmd2a,
c$$$     $           wvirk  / 3. / vol * fpmd2a,
c$$$     $           wvirkd / 3. / vol * fpmd2a,
c$$$     $           pressr * fpmd2a
            rio(8 * natoms + 7) = pressr
            rio(8 * natoms + 8) = dke
            call irdump(iudyn, rio, idbyte / nspbyt)
            write(99,*) rnkbv * rketmp * fpmd2a, 
     $           (wvirlj(1) + wvirlj(2) + wvirlj(3)) /
c$$$     $           (wvirlj(1) + wvirqt(1) + wvirqa(1)) / 
     $           3. / vol * fpmd2a, 
     $           (wvirqt(1) + wvirqt(2) + wvirqt(3)) /
c$$$     $           (wvirlj(2) + wvirqt(2) + wvirqa(2)) / 
     $           3. / vol * fpmd2a, 
     $           (wvirqa(1) + wvirqa(2) + wvirqa(3)) / 
c$$$     $           (wvirlj(3) + wvirqt(3) + wvirqa(3)) / 
     $           3. / vol * fpmd2a, 
     $           wvirk / 3. / vol * fpmd2a, 
     $           wvirkd / 3. / vol * fpmd2a,
     $           pressr * fpmd2a
c$$$            write(99,*) niters * dtbig, rketmp, rcmtmp
         endif
      endif
      if (moving .and. .not. qsiter .and. .not. dsiter .and. first) then
         first = .false.
      endif
c***********************************************************************
c     ndstep tells how many dynamics steps to run *after* the charge 
c     dynamics and equilibration phases are over:
c***********************************************************************
 9900 continue
      if (.not. eqflag .or. 
     $     ((vscflg .and. niters .lt. nscbas + ndstep) .or.
     $     (.not. vscflg .and. niters .lt. ndstep))) then
         go to 1000
      endif
      write(iuout, *) 'Dynamics finished at step ', niters
c***********************************************************************
c     We are now out of the dynamics loop.  Write out all of the 
c     pertinent simulation details so that analyze can read them later:
c***********************************************************************
      call wrtsim()
      call wrteor()
c***********************************************************************
c     write out hacked-together solvation energy data:
c***********************************************************************
      if (bdeflg) then
         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1),
     .         '/bonde.', suffix
         open (iutmp, file = filnam, err = 9990, iostat = ioval)
         do 8000 ibin = 1, nrclbn
            rwatcl = (ibin + 0.5) * rclbw
            if (nbdr(ibin) .gt. 1) then
               bdravg = bdercl(ibin) / nbdr(ibin)
c this is wrong.  it should be std = ssq / N - avg**2
c but the where_solv routine is counting on this mistake
c when it is fixed here, make sure to fix where_solv also
               bdrstd = bdrssq(ibin) / nbdr(ibin) - bdravg
               bdrstd = sqrt(bdrstd)
               write(iutmp, *) rwatcl, bdravg * fJm2kc, 
     .              bdrstd * fJm2kc, nbdr(ibin)
            else if (nbdr(ibin) .eq. 1) then
               write(iutmp, *) rwatcl, bdercl(ibin) * fJm2kc, 0.d0, 
     .              nbdr(ibin)
            else
               write(iutmp, *) rwatcl, 0.d0, 0.d0, nbdr(ibin)
            endif
 8000    continue
         close(iutmp)
      endif
c***********************************************************************
c     Print out all closing messages and shut down:
c***********************************************************************
      if (ioflg) then
            close(iudyn)
      endif
c$$$      if (dbflag) then
c$$$         close(21)
c$$$         close(22)
c$$$      endif
      close(99)
 9100 format(1x,'Charges not equilibrated after ',i6,' charge steps')
 9200 format(1x,'Charges equilibrated after ',i6,' charge steps')
 9300 format(1x,'RMS ',a,' energy conservation was ',g8.2)
      go to 9999
 9980 write(iuout, *) 'qdyn: error writing to ', filnam
      stop
 9990 write(iuout, *) 'qdyn: error opening ', filnam
      call ioerr(ioval)
      stop
 9999 continue
      end
