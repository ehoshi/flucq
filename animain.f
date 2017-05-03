      program animat
c***********************************************************************
c     This program reads in a qdyn trajectory and writes it out in
c     xyz format, for viewing with xmol....sjs 11/30/94
c***********************************************************************
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'anapar'
      include 'commons'
c
      character    symbol*2, filenm*80, fildir*80
      character*20 suffl(mdatfl)
      integer*4    nfrmid(maxmid), ntamid(maxmid),
     $             icjid(mxcstj,3)
      logical      anaflg, lspcok, inbox,
     $             knwmid(maxmid)
      real*4       rio(maxdat)
c
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
      neede = .false.
      lspcok = .true.
c***********************************************************************
c     Read in the input file animat.in
c***********************************************************************
c$$$      write(filnam, *) 'animat.in'
      write(filnam, '(a)') 'animat.in'
      open(iutmp, file = filnam, err = 9980, iostat = ioval)
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) fildir
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) datdir
      read(iutmp, *, err = 9980, end = 9991, iostat = ival) nsuff
      if (nsuff .gt. mdatfl) then
         write(iuout, *) 'animat: maximum number of data files ',
     $        '(mdatfl) exceeded'
      endif
      do 200 idatfl = 1, nsuff
         read(iutmp, *, err = 9980, end = 9991, iostat = ioval) 
     $        suffl(idatfl)
 200  continue
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) nsqint
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) startt
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) stopt
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) bfrac
      close(iutmp)
      bfrac = bfrac ** (1.0d0 / 3)
c***********************************************************************
c     Read in the header stuff and initialize related things:
c***********************************************************************
      write(filenm,'(3a)') datdir(1:index(datdir, ' ')-1), '/simdat.',
     $     suffl(1)
      write(iuout, *) 'reading header info from simdat.',suffl(1)
      open(iutmp2, file = filenm, err = 9980, iostat = ioval)
      call rddet(iutmp2, icjid)
      close(iutmp2)
      call chkdet(icjid)
      call setdet
      call rdmodb(knwmid, nfrmid, ntamid)
      write(stfile, *) 'eor.', suffl(1)
      write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
     $     '/eor.', suffl(1)
      anaflg = .true.
      call rdeor(filnam, knwmid, nfrmid, ntamid, anaflg)
      call setsys(anaflg)
      call rdmodd
      call setmod(icjid)
      if (stfile .eq. 'none' .and. syfile .eq. 'none') then
         write(iuout, *)
     $     'analyz: WARNING! system.dat may have changed...'
      endif
c***********************************************************************
c     Open the animation data file:
c***********************************************************************
      write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), '/animat.', 
     $     suffl(nsuff)
      open(20, file = filnam, err = 9980, iostat = ioval)
c***********************************************************************
c     loop through the data files:
c***********************************************************************
      iwrite = 0
      do 1200 idatfl = 1, nsuff
         write(filnam,'(3a)') fildir(1:index(fildir, ' ')-1),
     $        '/dyn.', suffl(idatfl)
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
            write(filnam,'(3a)') fildir(1:index(fildir, ' ')-1), 
     $           '/dyn.', suffl(idatfl+1)
            open(iutmp, file = filnam, status = 'old',
     $           form = 'unformatted', err = 9980, iostat = ioval)
            read(iutmp) iovers
            ioval = iurdmp(iutmp, rio, idbyte / nspbyt)
            if (ioval .eq. -1) then
               go to 9991
            endif
            istop = rio(8 * natoms + 5)
            close(iutmp)
         else
            istop = ibigp
         endif
c***********************************************************************
c     read records from this data file until the end, and output the
c     xyz coordinates:
c***********************************************************************
 720     continue
         ioval = iurdmp(iudyn, rio, idbyte / nspbyt)
         if (ioval .eq. -1) then
            go to 1100
         endif
         do 730 iatom = 1,natoms
            pos(iatom,1) = rio(iatom)
 730     continue
         do 735 iatom = 1,natoms
            pos(iatom,2) = rio(natoms + iatom)
 735     continue
         do 740 iatom = 1,natoms
            pos(iatom,3) = rio(2 * natoms + iatom)
 740     continue
         do 760 iatom = 1,natoms
            q(iatom) = rio(6 * natoms + iatom)
 760     continue
         niters = rio(8 * natoms + 5)
         if (pbflag) then 
            call foldr(imixed)
         else
            do 770 iatom = 1, natoms
               fldpos(iatom,1) = pos(iatom,1)
               fldpos(iatom,2) = pos(iatom,2)
               fldpos(iatom,3) = pos(iatom,3)
 770        continue
         endif
c***********************************************************************
c     If niters exceeds istop, then this run crashed, and the next
c     data file should be used, starting with this timestep:
c***********************************************************************
         if (niters .ge. istop) then
            write(iuout, *) suffl(idatfl)(1:index(suffl(idatfl),' ')),
     $           'crashed near t =', istop * dt
            go to 1100
         endif
c***********************************************************************
c     Now write out the xyz animation frame, if it has been nsqint
c     timesteps since the last write:
c***********************************************************************
c     Be careful: the calculation of whether a molecule is in the 
c     fraction of the box that we care about depends on how the unit
c     cell is defined, and this has been known to change.  It assumes
c     (-L/2,L/2) right now.
c***********************************************************************
         if (niters - iwrite .ge. nsqint .and.
     $        niters * dtbig .ge. startt .and.
     $        niters * dtbig .le. stopt) then
            iwrite = niters
            inatom = 0
            do 1010 imol = 1, nmol
               imty = molty(imol)
               inbox = .false.
               iatom = iatmol(imol,1)
               if (bfrac .eq. 1.0d0) then
                  inbox = .true.
               else if ((fldpos(iatom,1) + boxby2(1)) .le. 
     $                 bfrac * boxsiz(1) .and.
     $                 (fldpos(iatom,2) + boxby2(2)) .le. 
     $                 bfrac * boxsiz(2) .and.
     $                 (fldpos(iatom,3) + boxby2(3)) .le. 
     $                 bfrac * boxsiz(3)) then
                  inbox = .true.
               endif
               if (inbox) then
                  do 1000 imind = 1, natmty(imty)
                     iatom = iatmol(imol,imind)
                     if (ismasv(iatype(imty,imind))) then
                        inatom = inatom + 1
                     endif
 1000             continue
               endif
 1010       continue
            write(20, *) inatom
            write(20, *) ' '
            do 1020 imol = 1, nmol
               iatom = iatmol(imol,1)
               if (bfrac .eq. 1.d0) then
                  inbox = .true.
               else if ((fldpos(iatom,1) + boxby2(1)) .le. 
     $                 bfrac * boxsiz(1) .and.
     $                 (fldpos(iatom,2) + boxby2(2)) .le. 
     $                 bfrac * boxsiz(2) .and.
     $                 (fldpos(iatom,3) + boxby2(3)) .le. 
     $                 bfrac * boxsiz(3)) then
                  inbox = .true.
               endif
               if (inbox) then
                  do 1015 imind = 1, numatm(imol)
                     iatom = iatmol(imol,imind)
                     if (ident(iatom) .eq. 1) then
                        write(symbol, '(a2)') 'H'
                     else if (ident(iatom) .eq. 8) then
                        write(symbol, '(a2)') 'O'
                     else if (ident(iatom) .eq. 17) then
                        write(symbol, '(a2)') 'Cl'
                     else if (ident(iatom) .eq. 35) then
                        write(symbol, '(a2)') 'Br'
                     else
                        write(symbol, '(a2)') '?'
                     endif
                     if (ismasv(iatype(molty(imol),imind))) then
                        write(20, '(a2,1x,4f9.2)') symbol, 
     $                       (fldpos(iatom, i), i = 1, 3), q(iatom)
                     endif
 1015             continue
               endif
 1020       continue
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
 1200 continue
c***********************************************************************
c     close the output file:
c***********************************************************************
      close(20)
      write(iuout, *) 'All done.'
      stop
c***********************************************************************
c     Various error messages:
c***********************************************************************
 9980 write(iuout, *) 'analyz: error opening file ',filnam
      call ioerr(ioval)
      stop
 9990 write(iuout, *) 'analyz: error reading file ',filnam
      call ioerr(ioval)
      stop
 9991 write(iuout, *) 'analyz: unexpected end-of-file in',filnam
      stop
      end
