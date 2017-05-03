      program props
c***********************************************************************
c     The new version of qdyn does no binning of its own, so props does
c     all of the binning, collecting of statistics, output formatting,
c     etc....sjs 3/30/94
c***********************************************************************
c     
      include 'implic'      
      include 'genpar'
      include 'qpar'
      include 'propar'
      include 'commons'
c     
      character    name*80, 
c$$$     .     filenm*80, 
     .     fildir*80, eopfil*80,
     $     alpha*3
      character*20 suffl(mdatfl)
      integer*4    nfold(3),
     $     nfrmid(maxmid), ntamid(maxmid),
     $     icjid(mxcstj,3),
     $     iclmol(maxmol),
     $     nang(maxbin), nwthet(nangbn),
     $     nq(maxbin,maxaty),
     $     mxstnb(maxtbn,maxmty,maxmty), 
     $     nnbrt(maxtbn,maxmty,maxmty),
     $     nstnbr(maxtbn,maxmty,maxmty),
     $     nmz(maxbin,maxmty),
     $     ncg(maxbin,maxaty),
     $     ngz(maxbin,maxaty),
     $     ntg(maxbin,maxaty,maxaty),
     $     nnbr(maxmol,maxtbn),
     $     nbrcel(maxmol,maxmol,3),
     $     nbrlst(maxmol,mxnbr,maxtbn),
     $     nblscl(maxmol,mxnbr,maxtbn,3),
     $     nitbas(0:maxtim), nmsd(0:maxtim),
     $     ngrbst
      logical      anaflg, angflg, gorflg, gqrflg, 
     $     gcmflg, gzflg, mzflg, snpflg,
     $     lspcok,
     $     lnwbas, msdflg, mbasyt,
     $     qbnflg, trcflg, trsflg, clsflg, eopflg,
     $     rgyflg, corflg, pesflg, smuflg,
     $     gorwfl, gqrwfl, drflg,
     $     flag1, flag2, flag3, flag4, flag5, flag6, 
     $     flag7, flag8, flag9, flag10, flag11, flag12, flag13,
     $     flag14,
     $     getnr, getfar, getlt, gethvy, gete,
     $     cladd,
     $     knwmid(maxmid),
     $     incl(maxmol),
     $     isnbr(maxmol,maxmol)
      real*4       com(3), psmu(3), smusum(3), smussq(3),
     $     rnbr(maxmty,maxmty),
     $     qbwi(maxaty),
c$$$  $             gscal2(maxaty,maxaty),
     $     rmumin(maxmol), 
     $     rmu(maxmol,3),
     $     rio(maxdat),
     $     rsep(maxmol,maxmol)
      real*8       ang(maxbin), ang2(maxbin), 
     $     rthet(nangbn), rthet2(nangbn),
     $     rcg(maxbin,maxaty), rcg2(maxbin,maxaty),
     $     rgz(maxbin,maxaty), rgz2(maxbin,maxaty),
     $     rmz(maxbin,maxmty,3), rmz2(maxbin,maxmty,3),
     $     rg(maxbin,maxaty,maxaty), rg2(maxbin,maxaty,maxaty),
     $     rgq(maxbin,maxaty,maxaty), rgq2(maxbin,maxaty,maxaty),
     $     rmsd(0:maxtim), rmssd(0:maxtim),
     $     basmsd(maxatm,maxbas,3)
c     
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
      neede = .false.
      lspcok = .true.
      do 110 imid = 1, maxmid
         knwmid(imid) = .false.
  110 continue
c***********************************************************************
c     These are approximate correlation times for various statistics,
c     used to estimate errors:
c     angct = angular correlation time (for ang. corr. fxn.)
c     strct = structural correlation time (for g(r)s and related props)
c     genct = generic correlation time (for traces, etc.)
c***********************************************************************
      angct = 20000.d0
      strct = 2000.d0
      genct = 3000.d0
c***********************************************************************
c     Read in the input file props.in
c***********************************************************************
      write(filnam, '(a)') 'props.in'
      open(iutmp, file = filnam, err = 9980, iostat = ioval)
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) fildir
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) datdir
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) eopflg
      if (eopflg) then
         read(iutmp, *, err = 9980, end = 9991, iostat = ioval)
     $        eopfil
      endif
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) nsuff 
      if (nsuff .gt. mdatfl) then
         write(iuout, *) 'props: maximum number of data files ',
     $        '(mdatfl) exceeded'
         stop
      endif
      do 200 idatfl = 1, nsuff
         read(iutmp, *, err = 9980, end = 9991, iostat = ioval)
     $        suffl(idatfl)
  200 continue
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) nsqint
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) trcflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) qbnflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) angflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) gorflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) gqrflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) gcmflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) gzflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) mzflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) trsflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) clsflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) rgyflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) corflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) pesflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) smuflg
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) 
     $     snpflg, snptim
      read(iutmp, *, err = 9980, end = 9991, iostat = ioval) msdflg
      close(iutmp)
c***********************************************************************
c     Open the output files:
c***********************************************************************
      if (trcflg) then
         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
     .        '/temps.',
     $        suffl(nsuff)
         open(11, file = filnam, err = 9980, iostat = ioval)
         write(iuout, *) 'writing temperatures to temps.', suffl(nsuff)
c$$$  write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
c$$$  $        'qe.', suffl(nsuff)
c$$$  open(12, file = filnam, err = 9980, iostat = ioval)
c$$$  write(iuout, *) 'writing charge energies to qe.', suffl(nsuff)
         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
     $        '/fulle.', suffl(nsuff)
         open(20, file = filnam, err = 9980, iostat = ioval)
         write(iuout, *) 'writing full energies to fulle.', suffl(nsuff)
         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1),
     $        '/press.', suffl(nsuff)
         open(21, file = filnam, err = 9980, iostat = ioval)
         write(iuout, *) 'writing pressure to press.', suffl(nsuff)
c$$$         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1),
c$$$     $        '/special.', suffl(nsuff)
c$$$         open(30, file = filnam, err = 9980, iostat = ioval)
c$$$         write(iuout, *) 'writing com to special.', suffl(nsuff)
c$$$  write(iuout, *) 'writing q(1) to special.', 
c$$$  $        suffl(nsuff)
c$$$  write(iuout, *) 'writing Drude positions to special.', 
c$$$  $        suffl(nsuff)
c$$$  write(iuout, *) 'writing cluster sizes to special.',
c$$$  $        suffl(nsuff)
      endif
c***********************************************************************
c     Read in the header stuff and initialize related things:
c***********************************************************************
      if (pesflg) then
         neede = .true.
      endif
c$$$  write(filenm, '(3a)') datdir(1:index(datdir, ' ')-1), 
c$$$  $     '/simdat.', suffl(1)
      write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1),
     $     '/simdat.', suffl(1)
      write(iuout, *) 'reading header info from simdat.',suffl(1)
      open(iutmp, file = filnam, err = 9980, iostat = ioval)
c$$$  call rddet(filenm, icjid)
      call rddet(iutmp, icjid)
      close(iutmp)
      call chkdet(icjid)
      call setdet
      call rdmodb(knwmid, nfrmid, ntamid)
      write(stfile, *) '/eor.', suffl(1)
      write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
     $     '/eor.', suffl(1)
      anaflg = .true.
      call rdeor(filnam, knwmid, nfrmid, ntamid, anaflg)
      call setsys(anaflg)
      call chksys(icjid)
      call rdmodd
      call setmod(icjid)
      if (pesflg) then
         call tagup(imixed)
         call conchk
         call trmjst
         do 310 imty = 1, nmty
            call mpeset(imty)
  310    continue
         basepe = 0.
         do 320 imol = 1, nmol
            basepe = basepe + rmonpe(molty(imol))
  320    continue
         basepe = 2. * basepe
         if (stfile .eq. 'none' .and. syfile .eq. 'none') then
            write(iuout, *)
     $           'props: WARNING! system.dat may have changed...'
         endif
      endif
c***********************************************************************
c     set some variables:
c***********************************************************************
      iterio = ndstep / 20
      iterio = (iterio / nioint) * nioint
      niters = 0
c$$$  if (iterio .lt. 100) then
c$$$  iterio = (100 / nioint) * nioint
c$$$  if (iterio .lt. 100) then
c$$$  iterio = iterio + nioint
c$$$  endif
c$$$  endif
      if (iterio .lt. 1) then
         iterio = 1
      endif
      qtfac = 2.0d0 / (nqdof * k)
      Escale = fJm2kc / nmol
      fpby3 = 4.d0 * pi / 3.d0
      fpby3v = fpby3 / vol
c$$$  fpby3v = 4. * pi / (3. * vol)
      rbmn = 0.d0
      rbmx = dsqrt(3.0d0) * hbox
      rbw = (rbmx - rbmn) / maxbin
      rbwi = 1.d0 / rbw
      zbmx = boxsiz(3) * 0.5d0
      zbmn = -zbmx
      zbw = (zbmx - zbmn) / maxbin
      zbwi = 1.d0 / zbw
      if (qbnflg) then
         do 400 ityp1 = 1, naty
            qbwi(ityp1) = 1.0d0 / ((qhighg(ideaty(ityp1)) -
     $           qlowg(ideaty(ityp1))) / maxbin)
  400    continue
      endif
      if (gorflg) then
         do 430 iaty1 = 1, naty
            do 420 iaty2 = 1, naty
               do 410 ibin = 1, maxbin
                  ntg(ibin,iaty2,iaty1) = 0
  410          continue
  420       continue
  430    continue
      endif
      if (gcmflg) then
         do 450 iaty = 1, naty
            do 440 ibin = 1, maxbin
               ncg(ibin,iaty) = 0
  440       continue
  450    continue
      endif
      if (gzflg) then
         do 455 iaty = 1, naty
            do 452 ibin = 1, maxbin
               ngz(ibin,iaty) = 0
  452       continue
  455    continue
      endif
      if (hasdrd(molty(1))) then
         drflg = .true.
         drsum = 0.d0
         drssq = 0.d0
         dtsum = 0.d0
         dtssq = 0.d0
      endif
      if (corflg) then
         abmn = 0.d0
         abmx = pi
         abw = (abmx - abmn) / nangbn
         abwi = 1.d0 / abw
         do 460 ibin = 1, nangbn
            nwthet(ibin) = 0
  460    continue
      endif
      if (msdflg) then
         ioldbs = -9999999
         mbasyt = .false.
         nmsdba = 0
         do 470 itime = 0, maxtim
            rmsd(itime) = 0.d0
            rmssd(itime) = 0.d0
            nmsd(itime) = 0
  470    continue
      endif
      if (pesflg) then
         moving = .false.
         getnr = .true.
         getfar = .true.
         getlt = .true.
         gethvy = .true.
         gete = .true.
         pliter = .true.
      endif
c***********************************************************************
c     set the various counting and binning stats, either by reading the
c     eop file or zeroing them:
c***********************************************************************
      if (eopflg) then
         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
     $        '/', eopfil(1:index(eopfil, ' ')-1)
         write(iuout, *) 'reading prior binning info from ', filnam
         open(iutmp, file = filnam, status = 'old',
     $        form = 'unformatted', err = 9980, iostat = ioval)
         read(iutmp) iovers
         if (iovers .ne. ioeopv) then
            write(iuout, *) 'props: outdated eop file'
            stop
         endif
         read(iutmp) flag1, flag2, flag3, flag4, flag5, flag6, 
     $        flag7, flag8, flag9, flag10, flag11, flag12, flag13,
     $        flag14
         if ((flag1 .neqv. trcflg) .or. (flag2 .neqv. qbnflg) .or. 
     $        (flag3 .neqv. angflg) .or. (flag4 .neqv. gorflg) .or.
     $        (flag5 .neqv. gqrflg) .or. (flag6 .neqv. gcmflg) .or.
     $        (flag7 .neqv. gzflg) .or. (flag8 .neqv. mzflg) .or.
     $        (flag9 .neqv. trsflg) .or. (flag10 .neqv. clsflg) .or.
     $        (flag11 .neqv. rgyflg) .or. (flag12 .neqv. corflg) .or.
     $        (flag13 .neqv. pesflg) .or. (flag14 .neqv. smuflg)) then
            write(iuout, *) 'props: inconsistent input flags.  ',
     $           'flags should be:'
            write(iuout, *) 'trcflg: ', flag1
            write(iuout, *) 'qbnflg: ', flag2
            write(iuout, *) 'angflg: ', flag3
            write(iuout, *) 'gorflg: ', flag4
            write(iuout, *) 'gqrflg: ', flag5
            write(iuout, *) 'gcmflg: ', flag6
            write(iuout, *) 'gzflg:  ', flag7
            write(iuout, *) 'mzflg:  ', flag8
            write(iuout, *) 'trsflg: ', flag9
            write(iuout, *) 'clsflg: ', flag10
            write(iuout, *) 'rgyflg: ', flag11
            write(iuout, *) 'corflg: ', flag12
            write(iuout, *) 'pesflg: ', flag13
            write(iuout, *) 'smuflg: ', flag14
            stop
         endif
         read(iutmp) neopbs
         read(iutmp) ibin
         if (ibin .ne. maxbin) then
            write(iuout, *) 'props: maxbin has changed'
            stop
         endif
         if (qbnflg) then
            narrsz = maxbin * maxaty
            ioval = iuidmp(iutmp, nq, narrsz)
         endif
         if (angflg) then
            narrsz = maxbin
            ioval = iuidmp(iutmp, nang, narrsz)
            ioval = iuddmp(iutmp, ang, narrsz)
            ioval = iuddmp(iutmp, ang2, narrsz)
         endif
         if (gorflg) then
            narrsz = maxbin * maxaty * maxaty
            ioval = iuddmp(iutmp, rg, narrsz)
            ioval = iuddmp(iutmp, rg2, narrsz)
         endif
         if (gqrflg) then
            narrsz = maxbin * maxaty * maxaty
            ioval = iuddmp(iutmp, rgq, narrsz)
            ioval = iuddmp(iutmp, rgq2, narrsz)
         endif
         if (gcmflg) then
            narrsz = maxbin * maxaty
            ioval = iuddmp(iutmp, rcg, narrsz)
            ioval = iuddmp(iutmp, rcg2, narrsz)
         endif
         if (gzflg) then
            narrsz = maxbin * maxaty
            ioval = iuddmp(iutmp, rgz, narrsz)
            ioval = iuddmp(iutmp, rgz2, narrsz)
         endif
         if (mzflg) then
            narrsz = maxbin * maxmty * 3
            ioval = iuddmp(iutmp, rmz, narrsz)
            ioval = iuddmp(iutmp, rmz2, narrsz)
            narrsz = maxbin * maxmty
            ioval = iuidmp(iutmp, nmz, narrsz)
         endif
         if (trsflg) then
            narrsz = maxtbn * maxmty * maxmty
            ioval = iuidmp(iutmp, nstnbr, narrsz)
            ioval = iuidmp(iutmp, mxstnb, narrsz)
         endif
         if (corflg) then
            narrsz = nangbn
            ioval = iuddmp(iutmp, rthet, narrsz)
            ioval = iuddmp(iutmp, rthet2, narrsz)
            read(iutmp) ntthet
         endif
         read(iutmp) apesum, apessq, akesum, akessq, atesum, atessq,
     $        rktsum, rktssq, rttsum, rttssq, qtsum, qtssq, 
     $        dtsum, dtssq, psum, pssq
         if (clsflg) then
            read(iutmp) clssum, clsssq
         endif
         if (drflg) then
            read(iutmp) drsum, drssq
         endif
         if (rgyflg) then
            read(iutmp) rgysum, rgyssq
         endif
         if (smuflg) then
            read(iutmp) (smusum(i), i = 1,3), (smussq(i), i = 1,3)
         endif
         close(iutmp)
      else
         neopbs = 0
         if (trcflg) then
            apesum = 0.d0
            apessq = 0.d0
            akesum = 0.d0
            akessq = 0.d0
            atesum = 0.d0
            atessq = 0.d0
            rktsum = 0.d0
            rktssq = 0.d0
            rttsum = 0.d0
            rttssq = 0.d0
            qtsum = 0.d0
            qtssq = 0.d0
            psum = 0.d0
            pssq = 0.d0
            if (clsflg) then
               clssum = 0.d0
               clsssq = 0.d0
            endif
            if (drflg) then
               drsum = 0.d0
               drssq = 0.d0
            endif
            if (rgyflg) then
               rgysum = 0.d0
               rgyssq = 0.d0
            endif
         endif
         if (qbnflg) then
            do 490 ityp1 = 1, naty
               do 480 ibin = 1, maxbin
                  nq(ibin,ityp1) = 0
  480          continue
  490       continue
         endif
         if (angflg) then
            do 500 ibin = 1, maxbin
               ang(ibin) = 0.d0
               ang2(ibin) = 0.d0
               nang(ibin) = 0
  500       continue
         endif
         if (gorflg) then
            ngrbst = 0
            do 530 iaty1 = 1, naty
               do 520 iaty2 = 1, naty
                  do 510 ibin = 1, maxbin
                     rg(ibin,iaty2,iaty1) = 0.d0
                     rg2(ibin,iaty2,iaty1) = 0.d0
                     ntg(ibin,iaty2,iaty1) = 0
  510             continue
  520          continue
  530       continue
            if (gqrflg) then
               do 533 iaty1 = 1, naty
                  do 532 iaty2 = 1, naty
                     do 531 ibin = 1, maxbin
                        rgq(ibin,iaty2,iaty1) = 0.d0
                        rgq2(ibin,iaty2,iaty1) = 0.d0
  531                continue
  532             continue
  533          continue
            endif
         endif
         if (gcmflg) then
            do 535 iaty = 1, naty
               do 534 ibin = 1, maxbin
                  rcg(ibin,iaty) = 0.d0
                  rcg2(ibin,iaty) = 0.d0
                  ncg(ibin,iaty) = 0
  534          continue
  535       continue
         endif
         if (gzflg) then
            do 537 iaty = 1, naty
               do 536 ibin = 1, maxbin
                  rgz(ibin,iaty) = 0.d0
                  rgz2(ibin,iaty) = 0.d0
                  ngz(ibin,iaty) = 0
  536          continue
  537       continue
         endif
         if (mzflg) then
            if (.not. gzflg) then
               write(iuout, *) 'props: gzflg must be T to use mzflg'
               stop
            endif
            do 539 imty = 1, nmty
               do 538 ibin = 1, maxbin
                  rmz(ibin,imty,1) = 0.d0
                  rmz2(ibin,imty,1) = 0.d0
                  rmz(ibin,imty,2) = 0.d0
                  rmz2(ibin,imty,2) = 0.d0
                  rmz(ibin,imty,3) = 0.d0
                  rmz2(ibin,imty,3) = 0.d0
                  nmz(ibin,imty) = 0
  538          continue
  539       continue
         endif
         if (trsflg) then
            itrsbn = 0
            do 543 imty = 1, nmty
               do 542 imty2 = 1, imty
                  do 541 ibin = 1, maxtbn
                     mxstnb(ibin,imty,imty2) = 0
                     nstnbr(ibin,imty,imty2) = 0
  541             continue
  542          continue
  543       continue
         endif
         if (corflg) then
            do 544 ibin = 1, nangbn
               rthet(ibin) = 0.d0
               rthet2(ibin) = 0.d0
  544       continue
            ntthet = 0
         endif
         if (smuflg) then
            smusum(1) = 0.d0
            smusum(2) = 0.d0
            smusum(3) = 0.d0
            smussq(1) = 0.d0
            smussq(2) = 0.d0
            smussq(3) = 0.d0
         endif
      endif
c***********************************************************************
c     To calculate residence times or cluster sizes, we need to have a 
c     definition for how distant a pair can be and still be neighbors, 
c     or in the cluster.  These data are stored in gdata.in.  Read them
c     in, possibly scaling for the cluster size neighbor definition:
c***********************************************************************
      if (trsflg .or. clsflg .or. corflg .or. gzflg) then
         if (clsflg .and. trsflg) then
            write(iuout, *) 'props: WARNING! tres pair distance will ',
     $           'be used for cluster sizes.  cluster sizes will ',
     $           'appear too small'
         endif
         do 546 imty = 1, nmty
            do 545 imty2 = 1, imty
               rnbr(imty2,imty) = 0.
  545       continue
  546    continue
         write(filnam, '(a)') 'gdata.in'
         open(iutmp, file = filnam, err = 9980, iostat = ioval)
         read(iutmp, *, err = 9980, end = 9991, iostat = ioval) iovers
         if (iovers .ne. iomodv) then
            write(iuout, *) 'rdmodb:  gdata and models versions are ',
     $           'incompatible'
            stop
         endif
  547    continue
         read(iutmp, *, err = 9980, end = 550, iostat = ioval)
     $        imid1, imid2, trnbr
         if (hasmid(imid1) .and. hasmid(imid2)) then
            if (.not. trsflg) then
               rnbr(mtymid(imid1),mtymid(imid2)) = trnbr * 1.
               rnbr(mtymid(imid2),mtymid(imid1)) = trnbr * 1.
            else
               rnbr(mtymid(imid1),mtymid(imid2)) = trnbr
               rnbr(mtymid(imid2),mtymid(imid1)) = trnbr
            endif
         endif
         go to 547
  550    continue
         close(iutmp)
         do 570 imty = 1, nmty
            do 560 imty2 = 1, imty
               if (rnbr(imty2,imty) .eq. 0.) then
                  write(iuout, '(a,i2,a,i2,a,a)') 
     $                 'props: WARNING! molecules ', midmty(imty), 
     $                 ' and ', midmty(imty2), ' do not have a pair ',
     $                 'distance defined in gdata.in'
               endif
  560       continue
  570    continue
      endif
c***********************************************************************
c     Read in the first two records from the first data file, to get
c     the distance between timesteps:
c***********************************************************************
      write(filnam, '(3a)') fildir(1:index(fildir, ' ')-1), 
     $     '/dyn.', suffl(1)
      open(iutmp, file = filnam, status = 'old',
     $     form = 'unformatted', err = 9980, iostat = ioval)
      read(iutmp) iovers
      ioval = iurdmp(iutmp, rio, idbyte / nspbyt)
      if (ioval .eq. -1) then
         go to 9991
      endif
      itrone = rio(8 * natoms + 5)
      ioval = iurdmp(iutmp, rio, idbyte / nspbyt)
      if (ioval .eq. -1) then
         go to 9991
      endif
      itemp2 = rio(8 * natoms + 5)
      itrgap = itemp2 - itrone
      if (itrgap .le. 0) then
         write(iuout, *) 'props: 1st two records are weird...'
         stop
      endif
      close(iutmp)
c***********************************************************************
c     loop through the data files:
c***********************************************************************
      nrec = 0
      i1rec = 1
      iwrite = 0
      do 2310 idatfl = 1, nsuff
         write(filnam, '(3a)') fildir(1:index(fildir, ' ')-1),
     $        '/dyn.', suffl(idatfl)
         write(iuout,*) 'reading from dyn.',suffl(idatfl)
         open(iudyn, file = filnam, status = 'old',
     $        form = 'unformatted', err = 9980, iostat = ioval)
         read(iudyn) iovers
         if (iovers .ne. iodynv) then
            write(iuout, *) 'props: incompatible dyn file'
            stop
         endif
c***********************************************************************
c     read in the first record from the next data file, to see where it
c     starts, and avoid overlap:
c***********************************************************************
         if (idatfl .ne. nsuff) then
            write(filnam, '(3a)') fildir(1:index(fildir, ' ')-1), 
     $           '/dyn.', suffl(idatfl+1)
            open(iutmp, file = filnam, status = 'old',
     $           form = 'unformatted', err = 9980, iostat = ioval)
            read(iutmp) iovers
            ioval = iurdmp(iutmp, rio, idbyte / nspbyt)
            if (ioval .eq. -1) then
               go to 9991
            endif
            write(iuout, *) 'read istop'
            istop = rio(8 * natoms + 5)
            close(iutmp)
         else
            istop = ibigp
         endif
c***********************************************************************
c     read records from this data file until the end, and do stuff with
c     the numbers:
c***********************************************************************
  720    continue
         nrec = nrec + 1
         if (niters .eq. (niters / iterio) * iterio .and. niters .gt. 0)
     $        write(iuout, *) 'Timestep ', niters
         ioval = iurdmp(iudyn, rio, idbyte / nspbyt)
         if (ioval .eq. -1) then
            nrec = nrec - 1
            go to 2300
         endif
         do 730 iatom = 1,natoms
            temp = rio(iatom)
            pos(iatom,1) = temp
  730    continue
         do 735 iatom = 1,natoms
            temp = rio(natoms + iatom)
            pos(iatom,2) = temp
  735    continue
         do 740 iatom = 1,natoms
            temp = rio(2 * natoms + iatom)
            pos(iatom,3) = temp
  740    continue
         do 745 iatom = 1,natoms
            vel(iatom,1) = rio(3 * natoms + iatom)
  745    continue
         do 750 iatom = 1,natoms
            vel(iatom,2) = rio(4 * natoms + iatom)
  750    continue
         do 755 iatom = 1,natoms
            vel(iatom,3) = rio(5 * natoms + iatom)
  755    continue
c$$$  if (.not. pesflg) then
         do 760 iatom = 1,natoms
            q(iatom) = rio(6 * natoms + iatom)
  760    continue
c$$$  endif
         do 765 iatom = 1,natoms
            qvel(iatom) = rio(7 * natoms + iatom)
  765    continue
         qke = rio(8 * natoms + 1)
         qpe(imixed) = rio(8 * natoms + 2)
         ake = rio(8 * natoms + 3)
         ape(imixed) = rio(8 * natoms + 4)
         niters = rio(8 * natoms + 5)
         rscale = rio(8 * natoms + 6)
         pressr = rio(8 * natoms + 7)
         dke = rio(8 * natoms + 8)
c***********************************************************************
c     If niters exceeds istop, then this run crashed, and the next
c     data file should be used, starting with this timestep:
c***********************************************************************
         if (niters .ge. istop) then
            nrec = nrec - 1
            write(iuout, *) suffl(idatfl)(1:index(suffl(idatfl),' ')),
     $           'crashed near t =', istop * dtbig
            go to 2300
         endif
c***********************************************************************
c     Check the data spacing:
c***********************************************************************
         nrec2 = (niters - itrone) / itrgap + 1
         if (lspcok .and. nrec2 .ne. nrec) then
            lspcok = .false.
            write(iuout, *) 'props: WARNING! irregular data spacing'
            write(iuout, *) 't = ', nrec, ' should be ', nrec2
         endif
c***********************************************************************
c     take a snapshot and write it to an eor file if desired:
c***********************************************************************
         if (snpflg .and. niters * dtbig .ge. snptim) then
            write(iuout, *) 'taking a snapshot at t = ', niters * dtbig
            call wrtsnp
            snpflg = .false.
         endif
c***********************************************************************
c     If we're photoexciting anybody, recalculate the charges and the
c     system energy:
c***********************************************************************
         if (pesflg) then
c$$$  if (qslvfl) then
c$$$  qwrong = .true.
c$$$  endif
c$$$  if (dslvfl) then
c$$$  dwrong = .true.
c$$$  endif
c$$$  1100       continue
            movlt = .true.
            movhvy = .true.
            q(1) = 5.0d0
c$$$  if (qwrong .and. dwrong) then
c$$$  if (qsiter) then
c$$$  if (dslvfl) then
c$$$  dsiter = .true.
c$$$  endif
c$$$  qsiter = .false.
c$$$  else
c$$$  dsiter = .false.
c$$$  qsiter = .true.
c$$$  endif
c$$$  else if (dwrong) then
c$$$  dsiter = .true.
c$$$  qsiter = .false.
c$$$  else if (qwrong) then
c$$$  dsiter = .false.
c$$$  qsiter = .true.
c$$$  else
c$$$  dsiter = .false.
c$$$  qsiter = .false.
c$$$  endif
            call getf(getnr, getfar, getlt, gethvy, gete)
c$$$  if (qsiter) then
c$$$  call qsolv2
c$$$  rms = 0.d0
c$$$  do 1120 imol = 1, nmol
c$$$  imty = molty(imol)
c$$$  do 1110 iqind = 1, nfqmol(imty)
c$$$  iatom = iatmol(imol,ifqmol(imty,iqind))
c$$$  dq = q(iatom) - qans(iatom)
c$$$  rms = rms + dq * dq
c$$$  q(iatom) = qans(iatom)
c$$$  1110             continue
c$$$  1120          continue
c$$$  rms = sqrt(rms / nfqatm)
c$$$  if (rms .gt. qrmsct .and. dslvfl) then
c$$$  qwrong = .true.
c$$$  dwrong = .true.
c$$$  else
c$$$  qwrong = .false.
c$$$  endif
c$$$  endif
c$$$  if (dsiter) then
c$$$  call dsolv
c$$$  endif
c$$$  if (qsiter .or. dsiter) then
c$$$  go to 1100
c$$$  endif
         endif
         ape(imixed) = ape(iheavy) + ape(imixed) + ape(ilight)
         qpe(imixed) = qpe(iheavy) + qpe(imixed) + qpe(ilight)
c***********************************************************************
c     If temps were rescaled in this interval, reset all accumulation
c     and binning stats:
c***********************************************************************
         if (rscale .gt. 0.0d0) then
            if (niters - itrgap .ge. nscbas) then
               write(iuout, *) 'props: WARNING! temps were rescaled ',
     $              'after nscbas'
            endif
            i1rec = nrec
            if (trcflg) then
               apesum = 0.d0
               apessq = 0.d0
               akesum = 0.d0
               akessq = 0.d0
               atesum = 0.d0
               atessq = 0.d0
               rktsum = 0.d0
               rktssq = 0.d0
               rttsum = 0.d0
               rttssq = 0.d0
               qtsum = 0.d0
               qtssq = 0.d0
               dtsum = 0.d0
               dtssq = 0.d0
               psum = 0.d0
               pssq = 0.d0
               if (clsflg) then
                  clssum = 0.d0
                  clsssq = 0.d0
               endif
               if (drflg) then
                  drsum = 0.d0
                  drssq = 0.d0
               endif
               if (rgyflg) then
                  rgysum = 0.d0
                  rgyssq = 0.d0
               endif
            endif
            if (smuflg) then
               smusum(1) = 0.d0
               smusum(2) = 0.d0
               smusum(3) = 0.d0
               smussq(1) = 0.d0
               smussq(2) = 0.d0
               smussq(3) = 0.d0
            endif
         endif
c***********************************************************************
c     Calculate mean square displacement
c***********************************************************************
         if (msdflg) then
            if (niters - ioldbs .ge. 10 .and. nmsdba .lt. maxbas) then
               lnwbas = .true.
               ioldbs = niters
               nmsdba = nmsdba + 1
               nitbas(nmsdba) = niters
            else
               lnwbas = .false.
            endif
            if (niters - nitbas(1) .gt. maxtim) then
               write(iuout, *) 'step range = ', niters - nitbas(1),
     $            ' exceeds maxtim = ', maxtim
               write(iuout, *) 'increase maxtim and recompile'
               stop
            endif
            do 820 imol =  1, nmol
               imty = molty(imol)
c              Calculate the center of mass of the molecule
               com(1) = 0.d0
               com(2) = 0.d0
               com(3) = 0.d0
               rmlmas = 0.d0
               do 800 imaind = 1, nmamol(imty)
                  imind = imamol(imty,imaind)
                  iatom = iatmol(imol,imind)
                  iaty = iatype(imty,imind)
                  com(1) = com(1) + atmass(iaty) * pos(iatom,1)
                  com(2) = com(2) + atmass(iaty) * pos(iatom,2)
                  com(3) = com(3) + atmass(iaty) * pos(iatom,3)
                  rmlmas = rmlmas + atmass(iaty)
 800           continue
               com(1) = com(1) / rmlmas
               com(2) = com(2) / rmlmas
               com(3) = com(3) / rmlmas
               if (lnwbas) then
                  basmsd(imol,nmsdba,1) = com(1)
                  basmsd(imol,nmsdba,2) = com(2)
                  basmsd(imol,nmsdba,3) = com(3)
               endif
               do 810 imsdba = 1, nmsdba
                  dx = com(1) - basmsd(imol,imsdba,1)
                  dy = com(2) - basmsd(imol,imsdba,2)
                  dz = com(3) - basmsd(imol,imsdba,3)
                  sd = dx * dx + dy * dy + dz * dz
                  itrmsd = niters - nitbas(imsdba)
                  rmsd(itrmsd) = rmsd(itrmsd) + sd
                  rmssd(itrmsd) = rmssd(itrmsd) + sd * sd
                  nmsd(itrmsd) = nmsd(itrmsd) + 1
 810           continue
 820        continue
         endif
c***********************************************************************
c     Fold the positions back into the central box if necessary:
c***********************************************************************
         if (pbflag) then
            call foldr(imixed)
         endif
         if (niters .gt. nscbas) then
            ngrbst = ngrbst + 1
c***********************************************************************
c     Loop over molecule pairs, calculate intermolecular separations, 
c     and decide whether each molecule pair is a neighbor:
c***********************************************************************
            if (angflg .or. gorflg .or. trsflg .or. clsflg .or.
     $           corflg .or. gzflg) then
               do 1330 imol1 = 1, nmol
                  imty1 = molty(imol1)
                  do 1320 imol2 = 1, imol1 - 1
                     imty2 = molty(imol2)
                     call pbrhd(imol1, imol2, dxinc, dyinc, 
     $                    dzinc, nfold)
                     rsep(imol1,imol2) = rmat(1,1)
                     if (trsflg .or. clsflg .or. corflg 
     $                    .or. gzflg) then
                        if (rmat(1,1) .lt. rnbr(imty1,imty2)) then
                           isnbr(imol1,imol2) = .true.
                        else
                           isnbr(imol1,imol2) = .false.
                        endif
                        nbrcel(imol1,imol2,1) = nfold(1)
                        nbrcel(imol1,imol2,2) = nfold(2)
                        nbrcel(imol1,imol2,3) = nfold(3)
                     endif
c***********************************************************************
c     loop over atom pairs in this molecule pair and count the number in
c     each distance bin for g(r) purposes:
c***********************************************************************
                     if (gorflg) then
                        call pbrpr(imol1, imol2, 
     $                       dxinc, dyinc, dzinc, gqrflg)
                        do 1310 imind1 = 1, natmty(imty1)
c$$$  do 1310 irlin1 = 1, nrlmol(imty1)
c$$$  imind1 = irlmol(imty1,irlin1)
                           iatyp1 = iatype(imty1,imind1)
                           do 1300 imind2 = 1, natmty(imty2)
c$$$  do 1300 irlin2 = 1, nrlmol(imty2)
c$$$  imind2 = irlmol(imty2,irlin2)
                              iatyp2 = iatype(imty2,imind2)
                              if ((isrl(iatyp1) .and. 
     $                             isrl(iatyp2)) .or.
     $                             (gqrflg .and. 
     $                             (isrl(iatyp1) .and. 
     $                             isfqat(iatyp2)) .or.
     $                             (isfqat(iatyp1) .and. isrl(iatyp2)))) 
     $                             then
                                 ibin = int((rmat(imind1,imind2) - 
     $                                rbmn) *
     $                                rbwi) + 1
                                 if (ibin .ge. 1 .and. 
     $                                ibin .le. maxbin) then
                                    if (iatyp1 .le. iatyp2) then
                                       ntg(ibin,iatyp1,iatyp2) = 
     $                                      ntg(ibin,iatyp1,iatyp2) + 2
                                    else
                                       ntg(ibin,iatyp2,iatyp1) =
     $                                      ntg(ibin,iatyp2,iatyp1) + 2
                                    endif
                                    if (gqrflg) then
                                       if (isfqat(iatyp1)) then
                                          rgq(ibin,iatyp1,iatyp2) = 
     $                                         rgq(ibin,iatyp1,iatyp2) +
     $                                         q(iatmol(imol1,imind1))
                                          rgq2(ibin,iatyp1,iatyp2) = 
     $                                         rgq2(ibin,iatyp1,
     $                                         iatyp2) +
     $                                         q(iatmol(imol1,
     $                                         imind1)) ** 2
                                       endif
                                       if (isfqat(iatyp2)) then
                                          rgq(ibin,iatyp2,iatyp1) =
     $                                         rgq(ibin,iatyp2,iatyp1) +
     $                                         q(iatmol(imol2,imind2))
                                          rgq2(ibin,iatyp2,iatyp1) =
     $                                         rgq2(ibin,iatyp2,
     $                                         iatyp1) +
     $                                         q(iatmol(imol2,
     $                                         imind2)) ** 2
                                       endif
                                    endif
                                 endif
                              endif
 1300                      continue
 1310                   continue
                     endif
 1320             continue
 1330          continue
c***********************************************************************
c     now accumulate this step's g(r) data into the full tally:
c***********************************************************************
               if (gorflg) then
                  do 1420 iaty1 = 1, naty
                     do 1410 iaty2 = 1, iaty1
                        do 1400 ibin = 1, maxbin
                           rg(ibin,iaty2,iaty1) = rg(ibin,iaty2,iaty1) +
     $                          ntg(ibin,iaty2,iaty1)
                           rg2(ibin,iaty2,iaty1) = 
     $                          rg2(ibin,iaty2,iaty1) +
     $                          ntg(ibin,iaty2,iaty1) ** 2
                           ntg(ibin,iaty2,iaty1) = 0
 1400                   continue
 1410                continue
 1420             continue
               endif
c***********************************************************************
c     loop over maxtbn previous steps (less if there aren't that many
c     since equilibration), and count how many of the molecule pairs
c     that were neighbors in that step still are (and are still in the
c     same image cell):
c***********************************************************************
               if (trsflg) then
                  itrsbn = itrsbn + 1
 1480             continue
                  if (itrsbn .gt. maxtbn) then
                     itrsbn = itrsbn - maxtbn
                  else
                     go to 1490
                  endif
 1490             continue
                  nstep = nrec - i1rec
                  if (nstep .ge. maxtbn) then
                     nstep = maxtbn
                  endif
                  do 1540 istep = 1, nstep
                     ibin = itrsbn - istep
                     if (ibin .le. 0) then
                        ibin = ibin + maxtbn
                     endif
                     do 1510 imty = 1, nmty
                        do 1500 imty2 = 1, imty
                           mxstnb(istep,imty,imty2) = 
     $                          mxstnb(istep,imty,imty2) +
     $                          nnbrt(ibin,imty,imty2)
 1500                   continue
 1510                continue
                     do 1530 imol = 1, nmol
                        imty = molty(imol)
                        do 1520 inbri = 1, nnbr(imol,ibin)
                           inbr = nbrlst(imol,inbri,ibin)
                           imty2 = molty(inbr)
                           if (isnbr(imol,inbr)) then
c$$$  call pbrhd(imol, inbr, nfold)
                              if (nbrcel(imol,inbr,1) .eq. 
     $                             nblscl(imol,inbri,ibin,1) .and.
     $                             nbrcel(imol,inbr,2) .eq.
     $                             nblscl(imol,inbri,ibin,2) .and.
     $                             nbrcel(imol,inbr,3) .eq. 
     $                             nblscl(imol,inbri,ibin,3)) then
c$$$  if (isnbr(imol,inbr) .and. 
c$$$  $                       nblscl(imol,inbri,ibin,1) .eq. 
c$$$  $                       (nfold(imol,1) - nfold(inbr,1)) .and. 
c$$$  $                       nblscl(imol,inbri,ibin,2) .eq.
c$$$  $                       (nfold(imol,2) - nfold(inbr,2)) .and.
c$$$  $                       nblscl(imol,inbri,ibin,3) .eq.
c$$$  $                       (nfold(imol,3) - nfold(inbr,3))) then
                                 if (imty .ge. imty2) then
                                    nstnbr(istep,imty,imty2) =
     $                                   nstnbr(istep,imty,imty2) + 1
                                 else
                                    nstnbr(istep,imty2,imty) =
     $                                   nstnbr(istep,imty2,imty) + 1
                                 endif
                              endif
                           endif
 1520                   continue
 1530                continue
 1540             continue
c***********************************************************************
c     put this step's neighbor info into the neighbor list for use in
c     later steps:
c***********************************************************************
                  do 1610 imty = 1, nmty
                     do 1600 imty2 = 1, imty
                        nnbrt(itrsbn,imty,imty2) = 0
 1600                continue
 1610             continue
                  do 1630 imol = 1, nmol
                     imty = molty(imol)
                     nnbr(imol,itrsbn) = 0
                     do 1620 imol2 = 1, imol - 1
                        imty2 = molty(imol2)
                        if (isnbr(imol,imol2)) then
                           nnbr(imol,itrsbn) = nnbr(imol,itrsbn) + 1
                           if (nnbr(imol,itrsbn) .gt. mxnbr) then
                              write(iuout, '(a,i5,a,i2,a)') 
     $                             'props: molecule ', imol, 
     $                             ' has more than ', mxnbr, 
     $                             ' neighbors.  Increase mxnbr.'
                              stop
                           endif
c$$$  call pbrhd(imol, imol2, nfold)
                           nbrlst(imol,nnbr(imol,itrsbn),itrsbn) = imol2
                           nblscl(imol,nnbr(imol,itrsbn),itrsbn,1) =
     $                          nbrcel(imol,imol2,1)
c$$$  $                       nfold(1)
                           nblscl(imol,nnbr(imol,itrsbn),itrsbn,2) =
     $                          nbrcel(imol,imol2,2)
c$$$  $                       nfold(2)
                           nblscl(imol,nnbr(imol,itrsbn),itrsbn,3) =
     $                          nbrcel(imol,imol2,3)
c$$$  $                       nfold(3)
c$$$  nblscl(imol,nnbr(imol,itrsbn),itrsbn,1) =
c$$$  $                       nfold(imol,1) - nfold(imol2,1)
c$$$  nblscl(imol,nnbr(imol,itrsbn),itrsbn,2) =
c$$$  $                       nfold(imol,2) - nfold(imol2,3)
c$$$  nblscl(imol,nnbr(imol,itrsbn),itrsbn,3) =
c$$$  $                       nfold(imol,2) - nfold(imol2,3)
                           if (imty .ge. imty2) then
                              nnbrt(itrsbn,imty,imty2) =
     $                             nnbrt(itrsbn,imty,imty2) + 1
                           else
                              nnbrt(itrsbn,imty2,imty) =
     $                             nnbrt(itrsbn,imty2,imty) + 1
                           endif
                        endif
 1620                continue
 1630             continue
               endif
c***********************************************************************
c     Calculate the size of the cluster that contains the target 
c     molecule (hopefully a Cl- ion).  Do this by looping over molecules
c     to see which are neighbors of the target, adding these to the 
c     cluster.  Then loop again, adding everyone who is a neighbor of 
c     the bigger cluster.  Repeat until cluster stops growing.
c***********************************************************************
               if (clsflg .or. corflg .or. gzflg) then
                  do 1690 imol = 1, nmol
                     incl(imol) = .false.
 1690             continue
                  ntarg = 1
                  nclmol = 1
                  iclmol(nclmol) = ntarg
                  incl(ntarg) = .true.
 1700             continue
                  cladd = .false.
                  do 1720 itmol = 1, nmol
                     if (incl(itmol)) then
                        go to 1720
                     endif
                     do 1710 icli = 1, nclmol
                        if (incl(itmol)) then
                           go to 1710
                        endif
                        icmol = iclmol(icli)
                        if (itmol .gt. icmol) then
                           imol1 = itmol
                           imol2 = icmol
                        else
                           imol1 = icmol
                           imol2 = itmol
                        endif
                        if (isnbr(imol1,imol2)) then
                           nclmol = nclmol + 1
                           iclmol(nclmol) = itmol
                           incl(itmol) = .true.
                           cladd = .true.
                        endif
 1710                continue
 1720             continue
                  if (cladd) then
                     go to 1700
                  endif
                  clssum = clssum + nclmol
                  clsssq = clsssq + nclmol ** 2
               endif
            endif
c***********************************************************************
c     calculate the center of mass of the water-only cluster:
c***********************************************************************
            if (((gcmflg .or. corflg .or. rgyflg) .and. 
     $           .not. pbflag) .or. gzflg) then
               com(1) = 0.d0
               com(2) = 0.d0
               com(3) = 0.d0
               clmass = 0.d0
               do 1740 imol = 2, nmol
                  imty = molty(imol)
                  if (.not. incl(imol)) then
                     go to 1740
                  endif
                  do 1730 imaind = 1, nmamol(imty)
                     imind = imamol(imty,imaind)
                     iatom = iatmol(imol,imind)
                     iaty = iatype(imty,imind)
                     com(1) = com(1) + atmass(iaty) * pos(iatom,1)
                     com(2) = com(2) + atmass(iaty) * pos(iatom,2)
                     com(3) = com(3) + atmass(iaty) * pos(iatom,3)
                     clmass = clmass + atmass(iaty)
 1730             continue
 1740          continue
               com(1) = com(1) / clmass
               com(2) = com(2) / clmass
               com(3) = com(3) / clmass
            endif
c***********************************************************************
c     now get the g(r) of all cluster atoms to this center of mass:
c***********************************************************************
            if (gcmflg .and. .not. pbflag) then
               do 1760 imol = 1, nmol
                  imty = molty(imol)
                  if (.not. incl(imol)) then
                     go to 1760
                  endif
                  do 1750 irlind = 1, nrlmol(imty)
                     imind = irlmol(imty,irlind)
                     iaty = iatype(imty,imind)
                     iatom = iatmol(imol,imind)
                     dx = pos(iatom,1) - com(1)
                     dy = pos(iatom,2) - com(2)
                     dz = pos(iatom,3) - com(3)
                     r = sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                     ibin = int((r - rbmn) * rbwi) + 1
                     if (ibin .ge. 1 .and. ibin .le. maxbin) then
                        ncg(ibin,iaty) = ncg(ibin,iaty) + 2
                     endif
 1750             continue
 1760          continue
               do 1775 iaty = 1, naty
                  do 1770 ibin = 1, maxbin
                     rcg(ibin,iaty) = rcg(ibin,iaty) + ncg(ibin,iaty)
                     rcg2(ibin,iaty) = rcg2(ibin,iaty) +
     $                    ncg(ibin,iaty) ** 2
                     ncg(ibin,iaty) = 0
 1770             continue
 1775          continue
            endif
c***********************************************************************
c     and get the angular distribution of water oxygens about the 
c     Cl->c.o.m. vector:
c***********************************************************************
            if (corflg .and. .not. pbflag) then
               iclat = iatmol(1,1)
               cvecx = pos(iclat,1) - com(1)
               cvecy = pos(iclat,2) - com(2)
               cvecz = pos(iclat,3) - com(3)
               cvecm = sqrt(cvecx ** 2 + cvecy ** 2 + cvecz ** 2)
               do 1785 imol = 2, nmol
                  imty = molty(imol)
                  if (.not. incl(imol)) then
                     go to 1785
                  endif
                  iaty = iatype(imty,1)
                  iatom = iatmol(imol,1)
c$$$  dx = com(1) - pos(iatom,1)
c$$$  dy = com(2) - pos(iatom,2)
c$$$  dz = com(3) - pos(iatom,3)
                  dx = pos(iclat,1) - pos(iatom,1)
                  dy = pos(iclat,2) - pos(iatom,2)
                  dz = pos(iclat,3) - pos(iatom,3)
                  r = sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                  theta = acos((cvecx * dx + cvecy * dy + cvecz * dz) / 
     $                 cvecm / r)
                  ibin = int((theta - abmn) * abwi) + 1
                  if (ibin .ge. 1 .and. ibin .le. nangbn) then
                     nwthet(ibin) = nwthet(ibin) + 1
                     ntthet = ntthet + 1
                  endif
 1785          continue
               do 1790 ibin = 1, nangbn
                  rthet(ibin) = rthet(ibin) + nwthet(ibin)
                  rthet2(ibin) = rthet2(ibin) + nwthet(ibin) ** 2
                  nwthet(ibin) = 0
 1790          continue
            endif
c***********************************************************************
c     for the radius of gyration and density profile, use the real 
c     center of mass:
c***********************************************************************
            if ((rgyflg .and. .not. pbflag) .or. gzflg) then
               com(1) = com(1) * clmass
               com(2) = com(2) * clmass
               com(3) = com(3) * clmass
               imol = 1
               imty = molty(imol)
               do 1792 imaind = 1, nmamol(imty)
                  imind = imamol(imty,imaind)
                  iatom = iatmol(imol,imind)
                  iaty = iatype(imty,imind)
                  com(1) = com(1) + atmass(iaty) * pos(iatom,1)
                  com(2) = com(2) + atmass(iaty) * pos(iatom,2)
                  com(3) = com(3) + atmass(iaty) * pos(iatom,3)
                  clmass = clmass + atmass(iaty)
 1792          continue
               com(1) = com(1) / clmass
               com(2) = com(2) / clmass
               com(3) = com(3) / clmass
c$$$               write(30, *) niters * dtbig, (com(i), i = 1,3)
            endif
c***********************************************************************
c     get the density profile along the z-axis, with respect to the 
c     center of mass of the system.
c***********************************************************************
            if (gzflg) then
               do 1794 imol = 1, nmol
                  imty = molty(imol)
                  do 1793 irlind = 1, nrlmol(imty)
                     imind = irlmol(imty,irlind)
                     iaty = iatype(imty,imind)
                     iatom = iatmol(imol,imind)
                     ibin = int((fldpos(iatom,3) - com(3) - zbmn) 
     $                    * zbwi) + 1
                     if (ibin .ge. 1 .and. ibin .le. maxbin) then
                        ngz(ibin,iaty) = ngz(ibin,iaty) + 1
                     endif
 1793             continue
 1794          continue
               do 1796 iaty = 1, naty
                  do 1795 ibin = 1, maxbin
                     rgz(ibin,iaty) = rgz(ibin,iaty) + ngz(ibin,iaty)
                     rgz2(ibin,iaty) = rgz2(ibin,iaty) +
     $                    ngz(ibin,iaty) ** 2
                     ngz(ibin,iaty) = 0
 1795             continue
 1796          continue
            endif
c***********************************************************************
c     get the radius of gyration:
c***********************************************************************
            if (rgyflg .and. .not. pbflag) then
               nmaat = 0
               do 1798 imol = 1, nmol
                  imty = molty(imol)
                  if (.not. incl(imol)) then
                     go to 1798
                  endif
                  do 1797 imaind = 1, nmamol(imty)
                     nmaat = nmaat + 1
                     imind = imamol(imty,imaind)
                     iaty = iatype(imty,imind)
                     iatom = iatmol(imol,imind)
                     dx = pos(iatom,1) - com(1)
                     dy = pos(iatom,2) - com(2)
                     dz = pos(iatom,3) - com(3)
                     r2 = dx ** 2 + dy ** 2 + dz ** 2
                     rgy = rgy + r2
 1797             continue
 1798          continue
               rgy = sqrt(rgy / nmaat)
               rgysum = rgysum + rgy
               rgyssq = rgyssq + rgy ** 2
            endif
         endif
c***********************************************************************
c     Calculate the temperature:
c***********************************************************************
         call gettmp
c***********************************************************************
c     Calculate molecular dipoles:
c***********************************************************************
         if ((angflg .or. smuflg .or. mzflg) 
     $        .and. niters .ge. nscbas) then
            do 1810 imol = 1, nmol
               rmu(imol,1) = 0.d0
               rmu(imol,2) = 0.d0
               rmu(imol,3) = 0.d0
               imty = molty(imol)
               do 1800 iqind = 1, nqmol(imty)
                  iatom = iatmol(imol,iqmol(imty,iqind))
                  rmu(imol,1) = rmu(imol,1) + q(iatom) * fldpos(iatom,1)
                  rmu(imol,2) = rmu(imol,2) + q(iatom) * fldpos(iatom,2)
                  rmu(imol,3) = rmu(imol,3) + q(iatom) * fldpos(iatom,3)
 1800          continue
               ihead = iatmol(imol,iqmol(imty,1))
               rmu(imol,1) = rmu(imol,1) - qmol(imty) * fldpos(ihead,1)
               rmu(imol,2) = rmu(imol,2) - qmol(imty) * fldpos(ihead,2)
               rmu(imol,3) = rmu(imol,3) - qmol(imty) * fldpos(ihead,3)
               rmumin(imol) = sqrt(rmu(imol,1) ** 2 +
     $              rmu(imol,2) ** 2 + rmu(imol,3) ** 2)
 1810       continue
         endif
c***********************************************************************
c     Calculate system dipole:
c***********************************************************************
         psmu(1) = 0.d0
         psmu(2) = 0.d0
         psmu(3) = 0.d0
         do 1820 imol = 1, nmol
            psmu(1) = psmu(1) + rmu(imol,1)
            psmu(2) = psmu(2) + rmu(imol,2)
            psmu(3) = psmu(3) + rmu(imol,3)
 1820    continue
c***********************************************************************
c     calculate the distribution of dipoles along the z-axis:
c***********************************************************************
         if (mzflg) then
            do 1840 imol = 1, nmol
               imty = molty(imol)
               iatom = iatmol(imol,1)
               ibin = int((fldpos(iatom,3) - com(3) - zbmn) * zbwi) + 1
               if (ibin .ge. 1 .and. ibin .le. maxbin) then
                  rmz(ibin,imty,1) = rmz(ibin,imty,1) + rmu(imol,1)
                  rmz2(ibin,imty,1) = rmz2(ibin,imty,1) +
     $                 rmu(imol,1) ** 2
                  rmz(ibin,imty,2) = rmz(ibin,imty,2) + rmu(imol,2)
                  rmz2(ibin,imty,2) = rmz2(ibin,imty,2) + 
     $                 rmu(imol,2) ** 2
                  rmz(ibin,imty,3) = rmz(ibin,imty,3) + rmu(imol,3)
                  rmz2(ibin,imty,3) = rmz2(ibin,imty,3) +
     $                 rmu(imol,3) ** 2
                  nmz(ibin,imty) = nmz(ibin,imty) + 1
               endif
 1840       continue
         endif
c***********************************************************************
c     calculate charge distribution:
c***********************************************************************
         if (qbnflg .and. niters .ge. nscbas) then
            do 1860 imol = 1, nmol
               imty = molty(imol)
               do 1850 iqind = 1, nfqmol(imty)
                  imind = ifqmol(imty,iqind)
                  iatom = iatmol(imol,imind)
                  iatid = ident(iatom)
                  iattyp = iatype(imty,imind)
                  ibin = int((q(iatom) - qlowg(iatid)) * 
     $                 qbwi(iattyp)) + 1
                  if (ibin .ge. 1 .and. ibin .le. maxbin) then
                     nq(ibin,iattyp) = nq(ibin,iattyp) + 1
                  endif
 1850          continue
 1860       continue
         endif
c***********************************************************************
c     Calculate angular correlation functions:
c***********************************************************************
         if (angflg .and. niters .ge. nscbas) then
            do 1910 imol1 = 1, nmol
               do 1900 imol2 = 1, imol1 - 1
                  ibin = int((rsep(imol1,imol2) - rbmn) * rbwi) + 1
                  if (ibin .ge. 1 .and. ibin .le. maxbin) then
                     cosang = (rmu(imol1,1) * rmu(imol2,1) +
     $                    rmu(imol1,2) * rmu(imol2,2) +
     $                    rmu(imol1,3) * rmu(imol2,3)) / 
     $                    (rmumin(imol1) * rmumin(imol2))
                     ang(ibin) = ang(ibin) + cosang
                     ang2(ibin) = ang2(ibin) + cosang * cosang
                     nang(ibin) = nang(ibin) + 1
                  endif
 1900          continue
 1910       continue
         endif
c***********************************************************************
c     accumulate various stats:
c***********************************************************************
         apesum = apesum + ape(imixed)
         apessq = apessq + ape(imixed) ** 2
         akesum = akesum + ake
         akessq = akessq + ake ** 2
         atote = ake + ape(imixed)
         atesum = atesum + atote
         atessq = atessq + atote ** 2
         rktsum = rktsum + rketmp
         rktssq = rktssq + rketmp ** 2
         rttsum = rttsum + rcmtmp
         rttssq = rttssq + rcmtmp ** 2
         qtemp = qke * qtfac
         qtsum = qtsum + qtemp
         qtssq = qtssq + qtemp ** 2
         dtemp = dke * 2. / (3 * ndatom * k)
         dtsum = dtsum + dtemp
         dtssq = dtssq + dtemp ** 2
         psum = psum + pressr
         pssq = pssq + pressr ** 2
         if (drflg) then
            dx = pos(1,1) - pos(2,1)
            dy = pos(1,2) - pos(2,2)
            dz = pos(1,3) - pos(2,3)
            dr = sqrt(dx ** 2 + dy ** 2 + dz ** 2)
            drsum = drsum + dr
            drssq = drssq + dr ** 2
         endif
         if (smuflg) then
            smusum(1) = smusum(1) + psmu(1)
            smusum(2) = smusum(2) + psmu(2)
            smusum(3) = smusum(3) + psmu(3)
            smussq(1) = smussq(1) + psmu(1) ** 2
            smussq(2) = smussq(2) + psmu(2) ** 2
            smussq(3) = smussq(3) + psmu(3) ** 2
         endif
c***********************************************************************
c     write out traces:
c***********************************************************************
         if (trcflg .and. niters - iwrite .ge. nsqint) then
            iwrite = niters
            nbinst = nrec - i1rec + 1 + neopbs
c            nbinst = nrec - i1rec + neopbs
            stinef = genct / (itrgap * dtbig)
            if (stinef .lt. 1.d0) then
               stinef = 1.d0
            endif
            effbin = ceil(nbinst / stinef)
            stinef = nbinst / effbin
            rmean1 = qtsum / nbinst
            rmnsq1 = qtssq / nbinst
            stdev1 = rmnsq1 - rmean1 ** 2
            rmean2 = rktsum / nbinst
            rmnsq2 = rktssq / nbinst
            stdev2 = rmnsq2 - rmean2 ** 2
            rmean3 = rttsum / nbinst
            rmnsq3 = rttssq / nbinst
            stdev3 = rmnsq3 - rmean3 ** 2
            rmean4 = dtsum / nbinst
            rmnsq4 = dtssq / nbinst
            stdev4 = rmnsq3 - rmean3 ** 2
            errb1 = sqrt(stdev1 / (nbinst / stinef))
            errb2 = sqrt(stdev2 / (nbinst / stinef))
            errb3 = sqrt(stdev3 / (nbinst / stinef))
            errb4 = sqrt(stdev4 / (nbinst / stinef))
            write(11, '(13f15.6)') niters * dtbig, qtemp, rmean1, errb1,
     $           rketmp, rmean2, errb2, 
     $           rcmtmp, rmean3, errb3,
     $           dtemp,  rmean4, errb4
c$$$  write(12, *) niters * dtbig, qke * Escale,
c$$$  $           qpe(imixed) * Escale,
c$$$  $           (qke + qpe(imixed)) * Escale
            rmean1 = akesum / nbinst
            rmean2 = apesum / nbinst
            rmean3 = atesum / nbinst
            rmnsq1 = akessq / nbinst
            rmnsq2 = apessq / nbinst
            rmnsq3 = atessq / nbinst
            stdev1 = rmnsq1 - rmean1 ** 2
            stdev2 = rmnsq2 - rmean2 ** 2
            stdev3 = rmnsq3 - rmean3 ** 2
            errb1 = sqrt(stdev1 / (nbinst / stinef))
            errb2 = sqrt(stdev2 / (nbinst / stinef))
            errb3 = sqrt(stdev3 / (nbinst / stinef))
            write(20, '(10e21.13)') niters * dtbig, ake * Escale, 
     $           rmean1 * Escale, errb1 * Escale, 
     $           ape(imixed) * Escale, rmean2 * Escale, 
     $           errb1 * Escale, (ake + ape(imixed)) * Escale,
     $           rmean3 * Escale, errb1 * Escale
            rmean1 = psum / nbinst
            rmnsq1 = pssq / nbinst
            stdev1 = rmnsq1 - rmean1 ** 2
            errb1 = sqrt(stdev1 / (nbinst / stinef))
            write(21, '(4e21.13)') niters * dtbig, pressr * fpmd2a, 
     $           rmean1 * fpmd2a, errb1 * fpmd2a
c$$$  dx = pos(1,1) - pos(2,1)
c$$$  dy = pos(1,2) - pos(2,2)
c$$$  dz = pos(1,3) - pos(2,3)
c$$$  write(30, *) niters * dtbig, dx, dy, dz, 
c$$$  $           sqrt(dx ** 2 + dy ** 2 + dz ** 2)
c$$$  write(30, *) niters * dtbig, nclmol
c$$$  write(30, *) niters * dtbig, (q(4 * i), i = 1, 10)
         endif
c***********************************************************************
c     go back to read another record:
c***********************************************************************
         go to 720
c***********************************************************************
c     show up here at end of the current data file:
c***********************************************************************
 2300    continue
         close(iudyn)
 2310 continue
c***********************************************************************
c     write out correlation functions, averages, etc:
c***********************************************************************
      nbinst = nrec - i1rec + 1 + neopbs
c       nbinst = nrec - i1rec + neopbs
c***********************************************************************
c     write out the charge distributions:
c***********************************************************************
      if (qbnflg) then
         do 2350 iaty = 1, naty
            if (isfqat(iaty)) then
               iaid = ideaty(iaty)
               write(name, '(''q'', i2.2, ''.'')') iaid
               write(filnam, '(4a)') datdir(1:index(datdir, ' ')-1),
     $              '/', name(1:4), suffl(nsuff)
               open(iutmp, file = filnam, err = 9980, iostat = ioval)
               temp = qbwi(iaty) / (ngrbst * nwaty(iaty))
               do 2345 ibin = 1, maxbin
                  write(iutmp, *) qlowg(iaid) + 
     $                 (ibin - 0.5) / qbwi(iaty), nq(ibin,iaty) * temp
 2345          continue
               close(iutmp)
            endif
 2350    continue
      endif
c***********************************************************************
c     write out the angular correlation function:
c***********************************************************************
      if (angflg) then
         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
     $        '/angcor.', suffl(nsuff)
         open(iutmp, file = filnam, err = 9980, iostat = ioval)
         write(iuout, *) 'writing angular correlation function to ',
     $        'angcor.', suffl(nsuff)
         do 2400 ibin = 1, maxbin
            if (nang(ibin) .gt. 0.) then
               rmean = ang(ibin) / nang(ibin)
               rmnsq = ang2(ibin) / nang(ibin)
               stdev = rmnsq - rmean * rmean
               stinef = angct / (itrgap * dtbig)
               if (stinef .lt. 1.d0) then
                  stinef = 1.d0
               endif
               effbin = ceil(ngrbst / stinef)
               stinef = ngrbst / effbin
               errbar = sqrt(stdev / (nang(ibin) / stinef))
               write(iutmp, *) rbmn + (ibin - 0.5) * rbw, rmean, errbar
            endif
 2400    continue
         close(iutmp)
      endif
c$$$  c***********************************************************************
c$$$  c     calculate gscal2, a counting factor for the number of
c$$$  c     atom-type pairs:
c$$$  c***********************************************************************
      if (gorflg .or. gqrflg) then
c$$$  rideal = 4.0d0 * pi * (rndens / 3.d0) / 3.d0
c$$$  do 2530 iaty1 = 1, naty
c$$$  if (isrl(iaty1)) then
c$$$  imty1 = imtaty(iaty1)
c$$$  iaid1 = ideaty(iaty1)
c$$$  n1 = 0
c$$$  do 2500 imind1 = 1, natmty(imty1)
c$$$  if (iatype(imty1,imind1) .eq. iaty1) then
c$$$  n1 = n1 + 1
c$$$  endif
c$$$  2500          continue
c$$$  gscal2(iaty1,iaty1) = 1. / (n1 * n1)
c$$$  gscal2(iaty1,iaty1) = 1. / (nwaty(iaty1) ** 2)
c$$$  do 2520 iaty2 = 1, iaty1 - 1
c$$$  if (isrl(iaty2)) then
c$$$  imty2 = imtaty(iaty2)
c$$$  iaiad2 = ideaty(iaty2)
c$$$  n2 = 0
c$$$  do 2510 imind2 = 1, natmty(imty2)
c$$$  if (iatype(imty2,imind2) .eq. iaty2) then
c$$$  n2 = n2 + 1
c$$$  endif
c$$$  2510                continue
c$$$  gscal2(iaty1,iaty2) = 0.5 / (n1 * n2)
c$$$  gscal2(iaty1,iaty2) = 0.5 / 
c$$$  $                    (nwaty(iaty1) * nwaty(iaty2))
c$$$  endif
c$$$  2520          continue
c$$$  endif
c$$$  2530    continue
         do 2560 iaty1 = 1, naty
            if (isrl(iaty1) .or. isfqat(iaty1)) then
               do 2550 iaty2 = 1, iaty1
                  if (isrl(iaty2) .or. isfqat(iaty2)) then
                     if (.not. isrl(iaty1) .and. .not. isrl(iaty2)) then
                        go to 2550
                     endif
                     if (gorflg .and. isrl(iaty1) .and. isrl(iaty2)) 
     $                    then
                        gorwfl = .true.
                        if (ideaty(iaty1) .ge. ideaty(iaty2)) then
                           iaid1 = ideaty(iaty1)
                           iaid2 = ideaty(iaty2)
                        else
                           iaid1 = ideaty(iaty2)
                           iaid2 = ideaty(iaty1)
                        endif
                        write(name, '(''g'', 2i2.2, ''.'')') 
     $                       iaid1, iaid2
                        write(filnam, '(4a)') 
     $                       datdir(1:index(datdir, ' ')-1), '/',
     $                       name(1:6), suffl(nsuff)
                        open(iutmp, file = filnam, err = 9980, 
     $                       iostat = ioval)
                        gtemp = 1.d0 / 
     $                       (fpby3v * nwaty(iaty1) * nwaty(iaty2))
                        if (iaty1 .ne. iaty2) then
                           gtemp = gtemp / 2.d0
                        endif
                        rint = 0.d0
                        rint2 = 0.d0
c$$$  rmean = 0.d0
                     else
                        gorwfl = .false.
                     endif
                     if (gqrflg .and. 
     $                    isfqat(iaty1) .or. isfqat(iaty2)) then
                        gqrwfl = .true.
                        write(name, '(''gq'', 2i2.2, ''.'')') 
     $                       iaid1, iaid2
                        write(filnam, '(4a)')
     $                       datdir(1:index(datdir, ' ')-1), '/',
     $                       name(1:7), suffl(nsuff)
                        open(iutmp2, file = filnam, err = 9980,
     $                       iostat = ioval)
                     else
                        gqrwfl = .false.
                     endif
                     do 2540 ibin = 1, maxbin
                        gr = (ibin - 0.5) * rbw
                        stinef = strct / (itrgap * dtbig)
                        if (stinef .lt. 1.d0) then
                           stinef = 1.d0
                        endif
                        effbin = ceil(ngrbst / stinef)
                        stinef = ngrbst / effbin
                        if (gorwfl) then
                           if (cubflg .and. pbflag) then
                              if (gr .le. hbox) then
                                 gscale = 1.d0
                              else if (gr .le. hbox * dsqrt(2.d0)) then
                                 gscale = 1.d0 / (3 * hbox / gr - 2)
                              else
                                 gscale = 1.d0
                              endif
                           else
                              gscale = 1.d0
                           endif
c$$$  temp = gscale / (rideal * nmol *
c$$$  $                       ((ibin * rbw) ** 3 - 
c$$$  $                       ((ibin - 1) * rbw) ** 3))
c$$$  rmean = rg(ibin,iaty2,iaty1) * temp *
c$$$  $                       gscal2(iaty1,iaty2) / nbinst
c$$$  rmnsq = rg2(ibin,iaty2,iaty1) * temp ** 2 *
c$$$  $                       gscal2(iaty1,iaty2) ** 2 / nbinst
                           temp = gscale * gtemp / ((ibin * rbw) ** 3 -
     $                          ((ibin - 1) * rbw) ** 3)
c$$$  rmnold = rmean
                           rmean = rg(ibin,iaty2,iaty1) * temp / ngrbst
                           rmnsq = rg2(ibin,iaty2,iaty1) * temp ** 2 /
     $                         ngrbst
                           stdev = rmnsq - rmean * rmean
                           errbar = sqrt(stdev / (ngrbst / stinef))
                           if (iaty1 .eq. iaty2) then
                              rint = rint + rg(ibin,iaty2,iaty1) / 
     $                             (ngrbst * nwaty(iaty1))
                              rint2 = rint2 + rg(ibin,iaty2,iaty1) /
     $                             (ngrbst * nwaty(iaty2))
                           else
                              rint = rint + rg(ibin,iaty2,iaty1) /
     $                             (ngrbst * nwaty(iaty1) * 2.d0)
                              rint2 = rint2 + rg(ibin,iaty2,iaty1) /
     $                             (ngrbst * nwaty(iaty2) * 2.d0)
                           endif
                           write(iutmp, '(5f16.8)') 
     .                           gr, rmean, errbar, rint, rint2
                           if (iaid1 .eq. 8 .and. iaid2 .eq. 8) then
              endif
                        endif
                        if (gqrwfl) then
                           rmean1 = 0.d0
                           rmean2 = 0.d0
                           rmnsq1 = 0.d0
                           rmnsq2 = 0.d0
                           if (rg(ibin,iaty2,iaty1) .ne. 0) then
                              if (iaty1 .eq. iaty2) then
                                 rmean1 = (rgq(ibin,iaty2,iaty1) + 
     $                                rgq(ibin,iaty1,iaty2)) / 
     $                                rg(ibin,iaty2,iaty1)
                                 rmean2 = rmean1
                                 rmnsq1 = (rgq2(ibin,iaty2,iaty1) +
     $                                rgq2(ibin,iaty1,iaty2)) / 
     $                                rg(ibin,iaty2,iaty1)
                                 rmnsq2 = rmnsq1
                              else
                                 if (isfqat(iaty1)) then
                                    rmean1 = 2.d0 * 
     $                                   rgq(ibin,iaty1,iaty2) /
     $                                   rg(ibin,iaty2,iaty1)
                                    rmnsq1 = 2.d0 * 
     $                                   rgq2(ibin,iaty1,iaty2) /
     $                                   rg(ibin,iaty2,iaty1)
                                 endif
                                 if (isfqat(iaty2)) then
                                    rmean2 = 2.d0 * 
     $                                   rgq(ibin,iaty2,iaty1) /
     $                                   rg(ibin,iaty2,iaty1)
                                    rmnsq2 = 2.d0 * 
     $                                   rgq2(ibin,iaty2,iaty1) /
     $                                   rg(ibin,iaty2,iaty1)
                                 endif
                              endif
                           endif
                           stdev1 = rmnsq1 - rmean1 * rmean1
                           errb1 = sqrt(stdev1 / (ngrbst / stinef))
                           stdev2 = rmnsq2 - rmean2 * rmean2
                           errb2 = sqrt(stdev2 / (ngrbst / stinef))
                           write(iutmp2, *) gr, rmean1, errb1, 
     $                          rmean2, errb2
                        endif
 2540                continue
                     if (gorwfl) then
                        close(iutmp)
                     endif
                     if (gqrwfl) then
                        close(iutmp2)
                     endif
                  endif
 2550          continue
            endif
 2560    continue
      endif
      if (gcmflg) then
         do 2580 iaty = 1, naty
            if (.not. isrl(iaty)) then
               go to 2580
            endif
            write(name, '(''g'', i2.2, ''com.'')') ideaty(iaty)
            write(filnam, '(4a)') datdir(1:index(datdir, ' ')-1),
     $           '/', name(1:7), suffl(nsuff)
            open(iutmp, file = filnam, err = 9980, iostat = ioval)
            gtemp = 0.5d0 / (fpby3v * nwaty(iaty))
            rint = 0.d0
            rint2 = 0.d0
            do 2570 ibin = 1, maxbin
               gr = (ibin - 0.5) * rbw
               stinef = strct / (itrgap * dtbig)
               if (stinef .lt. 1.d0) then
                  stinef = 1.d0
               endif
               effbin = ceil(ngrbst / stinef)
               stinef = ngrbst / effbin
               temp = gtemp / ((ibin * rbw) ** 3 -
     $              ((ibin - 1) * rbw) ** 3)
               rmean = rcg(ibin,iaty) * temp / ngrbst
               rmnsq = rcg2(ibin,iaty) * temp ** 2 / ngrbst
               stdev = rmnsq - rmean ** 2
               errbar = sqrt(stdev / (ngrbst / stinef))
               rint = rint + rcg(ibin,iaty) / (ngrbst * 2.d0)
               write(iutmp, *) gr, rmean, errbar, rint
 2570       continue
            close(iutmp)
 2580    continue
      endif
      if (gzflg) then
         do 2600 iaty = 1, naty
            if (.not. isrl(iaty)) then
               go to 2600
            endif
            write(name, '(''gz'', i2.2, ''.'')') ideaty(iaty)
            write(filnam, '(4a)') datdir(1:index(datdir, ' ')-1),
     $           '/', name(1:5), suffl(nsuff)
            open(iutmp, file = filnam, err = 9980, iostat = ioval)
            temp = 1.d0 / (boxsiz(1) * boxsiz(2) * zbw)
            do 2590 ibin = 1, maxbin
               gz = (ibin - 0.5) * zbw + zbmn
               stinef = strct / (itrgap * dtbig)
               if (stinef .lt. 1.d0) then
                  stinef = 1.d0
               endif
               effbin = ceil(ngrbst / stinef)
               stinef = ngrbst / effbin
               rmean = rgz(ibin,iaty) * temp / ngrbst
               rmnsq = rgz2(ibin,iaty) * temp ** 2 / ngrbst
               stdev = rmnsq - rmean ** 2
               errbar = sqrt(stdev / (ngrbst / stinef))
               write(iutmp, *) gz, rmean, errbar
 2590       continue
            close(iutmp)
 2600    continue
      endif
      if (mzflg) then
         do 2620 imty = 1, nmty
            imid = midmty(imty)
            if (imid .eq. 2) then
               alpha = 'T '
            else if (imid .eq. 4) then
               alpha = 'qT '
            else if (imid .eq. 5) then
               alpha = 'LC '
            else if (imid .eq. 6) then
               alpha = 'DC '
            else
               alpha = 'un '
            endif
            write(name, '(''mz'', a, ''.'')') 
     $           alpha(1:index(alpha, ' ')-1)
            write(filnam, '(4a)') datdir(1:index(datdir, ' ')-1),
     $           '/', name(1:index(name, '.')), suffl(nsuff)
            open(iutmp, file = filnam, err = 9980, iostat = ioval)
            do 2610 ibin = 1, maxbin
               gz = (ibin - 0.5) * zbw + zbmn
               if (nmz(ibin,imty) .gt. 0) then
                  stinef = strct / (itrgap * dtbig)
                  if (stinef .lt. 1.d0) then
                     stinef = 1.d0
                  endif
                  effbin = ceil(ngrbst / stinef)
                  stinef = ngrbst / effbin
                  rmean1 = rmz(ibin,imty,1) / nmz(ibin,imty)
                  rmnsq1 = rmz2(ibin,imty,1) / nmz(ibin,imty)
                  stdev1 = rmnsq1 - rmean1 ** 2
                  errbr1 = sqrt(stdev1 / (ngrbst / stinef))
                  rmean2 = rmz(ibin,imty,2) / nmz(ibin,imty)
                  rmnsq2 = rmz2(ibin,imty,2) / nmz(ibin,imty)
                  stdev2 = rmnsq2 - rmean2 ** 2
                  errbr2 = sqrt(stdev2 / (ngrbst / stinef))
                  rmean3 = rmz(ibin,imty,3) / nmz(ibin,imty)
                  rmnsq3 = rmz2(ibin,imty,3) / nmz(ibin,imty)
                  stdev3 = rmnsq3 - rmean3 ** 2
                  errbr3 = sqrt(stdev3 / (ngrbst / stinef))
               else
                  rmean1 = 0.d0
                  errbr1 = 0.d0
                  rmean2 = 0.d0
                  errbr2 = 0.d0
                  rmean3 = 0.d0
                  errbr3 = 0.d0
               endif
               write(iutmp, *) gz, rmean1 * feA2D, errbr1 * feA2D, 
     $              rmean2 * feA2D, errbr2 * feA2D, 
     $              rmean3 * feA2D, errbr3 * feA2D
 2610       continue
            close(iutmp)
 2620    continue
      endif
c***********************************************************************
c     normalize and write out the residence time stuff:
c***********************************************************************
      if (trsflg) then
         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
     $        '/tres.', suffl(nsuff)
         open(iutmp, file = filnam, err = 9980, iostat = ioval)
         write(iuout, *) 'writing residence time curve to tres.',
     $        suffl(nsuff)
         nstep = nrec - i1rec
         if (nstep .ge. maxtbn) then
            nstep = maxtbn
         endif
         write(iutmp, *) 0.d0, ((1.d0, j = 1, i), i = 1, nmty)
         do 2700 istep = 1, nstep
            write(iutmp, *) istep * itrgap * dtbig, 
     $           ((real(nstnbr(istep,i,j)) / mxstnb(istep,i,j), 
     $           j = 1, i), i = 1, nmty)
 2700    continue
      endif
c***********************************************************************
c     normalize and write out the angle dist about Cl->com vector:
c***********************************************************************
      if (corflg) then
         stinef = strct / (itrgap * dtbig)
         if (stinef .lt. 1.d0) then
            stinef = 1.d0
         endif
         effbin = ceil(ngrbst / stinef)
         stinef = ngrbst / effbin
         rint = 0.d0
         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
     $        '/cortis.', suffl(nsuff)
         open(iutmp, file = filnam, err = 9980, iostat = ioval)
         do 2800 ibin = 1, nangbn
            theta = (ibin - 0.5) * abw
            dtheta = cos(theta - .5 * abw) - cos(theta + .5 * abw)
            scale = 1.d0 / (((dble(ntthet) / ngrbst) / 4.d0 / pi) * 
     $           2.d0 * pi * dtheta)
            rmean = scale * rthet(ibin) / ngrbst
            rmnsq = scale ** 2 * rthet2(ibin) / ngrbst
            stdev = rmnsq - rmean * rmean
            errbar = sqrt(stdev / (ngrbst / stinef))
            rint = rint + rmean * 2.d0 * pi * dtheta / (4 * pi) * 
     $           (dble(ntthet) / ngrbst)
            write(iutmp, *) theta * fr2deg, rmean, errbar, rint
 2800    continue
         close(iutmp)
      endif
      if (msdflg) then
         write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1),
     $      '/msd.', suffl(nsuff)
         open(iutmp, file = filnam, err = 9980, iostat = ioval)
         write(iuout, *) 'writing mean square displacement data to ',
     $      'msd.', suffl(nsuff)
c        calculate the diffusion coefficieint using weighted linear
c        least squares fit to the mean squared displacement
         smsd = 0.d0
         sxmsd = 0.d0
         symsd = 0.d0
         sxxmsd = 0.d0
         sxymsd = 0.d0
         syymsd = 0.d0
         nxmsd = 0
         write(iutmp, *) 0.d0, rmsd(0), 0.d0
         do 2900 imsdti = 1, maxtim
            if (nmsd(imsdti) .gt. 0) then
               xmsd = imsdti * dtbig
               rmsd(imsdti) = rmsd(imsdti) / nmsd(imsdti)
               rmssd(imsdti) = rmssd(imsdti) / nmsd(imsdti)
               varmsd = rmssd(imsdti) - rmsd(imsdti) * rmsd(imsdti)
               stdmsd = sqrt(varmsd)
               errmsd = stdmsd / sqrt(dble(nmsd(imsdti)))
               write(iutmp, *) xmsd, rmsd(imsdti), errmsd
c              don't use the first 1 ps or so to calculate the 
c              diffusion coefficient
               if (xmsd .ge. 10000) then
                  wtmsd = 1.d0 / errmsd / errmsd
                  smsd = smsd + wtmsd
                  sxmsd = sxmsd + xmsd * wtmsd
                  symsd = symsd + rmsd(imsdti) * wtmsd
                  sxxmsd = sxxmsd + xmsd * xmsd * wtmsd
                  sxymsd = sxymsd + xmsd * rmsd(imsdti) * wtmsd
                  syymsd = syymsd + rmsd(imsdti) * rmsd(imsdti) * wtmsd
                  nxmsd = nxmsd + 1
               endif
            endif
 2900    continue
      endif
c***********************************************************************
c     write out some miscellaneous properties:
c***********************************************************************
      write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
     $     '/props.', suffl(nsuff)
      open(iutmp, file = filnam, err = 9980, iostat = ioval)
      write(iuout, *) 'writing various properties to props.',
     $     suffl(nsuff)
      stinef = genct / (itrgap * dtbig)
      if (stinef .lt. 1.d0) then
         stinef = 1.d0
      endif
      effbin = ceil(nbinst / stinef)
      stinef = nbinst / effbin
      apeave = apesum / nbinst
      akeave = akesum / nbinst
      ateave = atesum / nbinst
      rktave = rktsum / nbinst
      rttave = rttsum / nbinst
      qtave = qtsum / nbinst
      dtave = dtsum / nbinst
      pave = psum / nbinst
      apestd = sqrt(apessq / nbinst - apeave ** 2)
      akestd = sqrt(akessq / nbinst - akeave ** 2)
      atestd = sqrt(atessq / nbinst - ateave ** 2)
      rktstd = sqrt(rktssq / nbinst - rktave ** 2)
      rttstd = sqrt(rttssq / nbinst - rttave ** 2)
      qtstd = sqrt(qtssq / nbinst - qtave ** 2)
      dtstd = sqrt(dtssq / nbinst - dtave ** 2)
      pstd = sqrt(pssq / nbinst - pave ** 2)
      rmean = apesum / nbinst
      rmnsq = apessq / nbinst
      stdev = rmnsq - rmean * rmean
      errbar = sqrt(stdev / (nbinst / stinef))
      write(iutmp, '(a,g10.4,a,g9.3,a)') 
     $     '<PE> = ', rmean * Escale, ' +/- ', errbar * Escale,
     $     ' kcal/mol'
      rmean = akesum / nbinst
      rmnsq = akessq / nbinst
      stdev1 = rmnsq - rmean * rmean
      errbar = sqrt(stdev1 / (nbinst / stinef))
      write(iutmp, '(a,g10.4,a,g9.3,a)') 
     $     '<KE> = ', rmean * Escale, ' +/- ', errbar * Escale,
     $     ' kcal/mol'
      rmean = atesum / nbinst
      rmnsq = atessq / nbinst
      stdev = rmnsq - rmean * rmean
      errbar = sqrt(stdev / (nbinst / stinef))
      write(iutmp, '(a,g10.4,a,g10.4,a)') 
     $     '<E>  = ', rmean * Escale, ' +/- ', errbar * Escale, 
     $     ' kcal/mol'
      write(iutmp, '(a,g10.4)')
     $     'sqrt(<dE^2>/<dKE^2>) = ', sqrt(stdev / stdev1)
      rmean = rktsum / nbinst
      rmnsq = rktssq / nbinst
      stdev = rmnsq - rmean * rmean
      errbar = sqrt(stdev / (nbinst / stinef))
      write(iutmp, '(a,g10.4,a,g9.3,a)')
     $     '<T>  = ', rmean, ' +/- ', errbar, ' K'
      rmean = rttsum / nbinst
      rmnsq = rttssq / nbinst
      stdev = rmnsq - rmean ** 2
      errbar = sqrt(stdev / (nbinst / stinef))
      write(iutmp, '(a,g10.4,a,g9.3,a)')
     $     '<tT> = ', rmean, ' +/- ', errbar, ' K'
      rmean = qtsum / nbinst
      rmnsq = qtssq / nbinst
      stdev = rmnsq - rmean * rmean
      errbar = sqrt(stdev / (nbinst / stinef))
      write(iutmp, '(a,g10.4,a,g9.3,a)')
     $     '<qT> = ', rmean, ' +/- ', errbar, ' K'
      if (hasdrd(molty(1))) then
         rmean = dtsum / nbinst
         rmnsq = dtssq / nbinst
         stdeve = rmnsq - rmean * rmean
         errbar = sqrt(stdev / (nbinst / stinef))
         write(iutmp, '(a,g10.4,a,g9.3,a)')
     $        '<dT> = ', rmean, ' +/- ', errbar, ' K'
      endif
      rmean = psum / nbinst
      rmnsq = pssq / nbinst
      stdev = rmnsq - rmean ** 2
      errbar = sqrt(stdev / (nbinst / stinef))
      write(iutmp, '(a,g10.4,a,g9.3,a)')
     $     '<P> = ', rmean * fpmd2a, ' +/- ', errbar * fpmd2a, ' atm'
      if (clsflg) then
         rmean = clssum / nbinst
         rmnsq = clsssq / nbinst
         stdev = rmnsq - rmean ** 2
         errbar = sqrt(stdev / (nbinst / stinef))
         write(iutmp, '(a,g10.4,a,g9.3)')
     $        'avg. cluster size = ', rmean, ' +/- ', errbar
      endif
      if (drflg) then
         rmean = drsum / nbinst
         rmnsq = drssq / nbinst
         stdev = rmnsq - rmean ** 2
         errbar = sqrt(stdev / (nbinst / stinef))
         write(iutmp, '(a,g10.4,a,g9.3,a)')
     $        '<r> = ', rmean, ' +/- ', errbar, ' A'
      endif
      if (rgyflg) then
         rmean = rgysum / nbinst
         rmnsq = rgyssq / nbinst
         stdev = rmnsq - rmean ** 2
         errbar = sqrt(stdev / (nbinst / stinef))
         write(iutmp, '(a,g10.4,a,g9.3,a,g9.3,a)')
     $        '<Rg> = ', rmean, ' +/- ', errbar, ', stdv = ', 
     $        sqrt(stdev), ' A'
      endif
      if (smuflg) then
         rmean = smusum(1) / nbinst
         rmnsq = smussq(1) / nbinst
         stdev = rmnsq - rmean ** 2
         errbar = sqrt(stdev / (nbinst / stinef))
         write(iutmp, '(a,g10.4,a,g9.3,a,g9.3,a)')
     $        '<Mx> = ', rmean * feA2D, ' +/- ', errbar * feA2D, ', 
     $        stdev = ', sqrt(stdev) * feA2D, ' D'
         rmean = smusum(2) / nbinst
         rmnsq = smussq(2) / nbinst
         stdev = rmnsq - rmean ** 2
         errbar = sqrt(stdev / (nbinst / stinef))
         write(iutmp, '(a,g10.4,a,g9.3,a,g9.3,a)')
     $        '<My> = ', rmean * feA2D, ' +/- ', errbar * feA2D, ', 
     $        stdev = ', sqrt(stdev) * feA2D, ' D'
         rmean = smusum(3) / nbinst
         rmnsq = smussq(3) / nbinst
         stdev = rmnsq - rmean ** 2
         errbar = sqrt(stdev / (nbinst / stinef))
         write(iutmp, '(a,g10.4,a,g9.3,a,g9.3,a)')
     $        '<Mz> = ', rmean * feA2D, ' +/- ', errbar * feA2D, ', 
     $        stdev = ', sqrt(stdev) * feA2D, ' D'
      endif
      if (msdflg) then
         dmsd = (smsd * sxymsd - sxmsd * symsd)
     $          / (smsd * sxxmsd - sxmsd * sxmsd)
         bmsd = symsd / smsd - dmsd * sxmsd / smsd
c        not sure about this one for weighted LS...
         sbmsd = sqrt((syymsd - 2.d0 * dmsd * sxymsd
     $                 - 2.d0 * bmsd * symsd + dmsd * dmsd * sxxmsd
     $                 + 2.d0 * dmsd * bmsd * sxmsd
     $                 + bmsd * bmsd * smsd) / smsd
     $                 / (dble(nxmsd) / (30000. / dtbig))
     $                / (sxxmsd / smsd - sxmsd * sxmsd / smsd / smsd))
         write(6, *) 'smsd = ', smsd
         write(iutmp, '(a,g10.4,a,g10.4,a,g10.4,a)')
     $        '<D> = ', dmsd / 6.d0, ' +/- ', sbmsd / 6.d0,
     $        ' cm^2/s, ballistic offset = ', bmsd, ' A^2'
      endif
c***********************************************************************
c     write out the status file (eop.SUF)
c***********************************************************************
      write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
     $     '/eop.', suffl(nsuff)
      write(iuout, *) 'writing binning status to ', filnam
      open(iutmp, file = filnam, form = 'unformatted', err = 9980,
     $     iostat = ioval)
      write(iutmp) ioeopv
      write(iutmp) trcflg, qbnflg, angflg, gorflg, gqrflg, gcmflg, 
     $     gzflg, mzflg, trsflg, clsflg, rgyflg, corflg, pesflg, smuflg
      write(iutmp) nbinst
      write(iutmp) maxbin
      if (qbnflg) then
         narrsz = maxbin * maxaty
         call iidump(iutmp, nq, narrsz)
      endif
      if (angflg) then
         narrsz = maxbin
         call iidump(iutmp, nang, narrsz)
         call iddump(iutmp, ang, narrsz)
         call iddump(iutmp, ang2, narrsz)
      endif
      if (gorflg) then
         narrsz = maxbin * maxaty * maxaty
         call iddump(iutmp, rg, narrsz)
         call iddump(iutmp, rg2, narrsz)
      endif
      if (gqrflg) then
         narrsz = maxbin * maxaty * maxaty
         call iddump(iutmp, rgq, narrsz)
         call iddump(iutmp, rgq2, narrsz)
      endif
      if (gcmflg) then
         narrsz = maxbin * maxaty
         call iddump(iutmp, rcg, narrsz)
         call iddump(iutmp, rcg2, narrsz)
      endif
      if (gzflg) then
         narrsz = maxbin * maxaty
         call iddump(iutmp, rgz, narrsz)
         call iddump(iutmp, rgz2, narrsz)
      endif
      if (mzflg) then
         narrsz = maxbin * maxmty * 3
         call iddump(iutmp, rmz, narrsz)
         call iddump(iutmp, rmz2, narrsz)
         narrsz = maxbin * maxmty
         call iidump(iutmp, nmz, narrsz)
      endif
      if (trsflg) then
         narrsz = maxtbn * maxmty * maxmty
         call iidump(iutmp, nstnbr, narrsz)
         call iidump(iutmp, mxstnb, narrsz)
      endif
      if (corflg) then
         narrsz = nangbn
         call iddump(iutmp, rthet, narrsz)
         call iddump(iutmp, rthet2, narrsz)
         write(iutmp) ntthet
      endif
      write(iutmp) apesum, apessq, akesum, akessq, atesum, atessq,
     $     rktsum, rktssq, rttsum, rttssq, qtsum, qtssq, dtsum, dtssq
      if (clsflg) then
         write(iutmp) clssum, clsssq
      endif
      if (drflg) then
         write(iutmp) drsum, drssq
      endif
      if (rgyflg) then
         write(iutmp) rgysum, rgyssq
      endif
      if (smuflg) then
         write(iutmp) (smusum(i), i = 1,3), (smussq(i), i = 1,3)
      endif
      close(iutmp)
c***********************************************************************
c     Close the output files and finish up:
c***********************************************************************
      if (trcflg) then
         close(11)
c$$$  close(12)
         close(20)
c$$$         close(30)
      endif
      stop
c***********************************************************************
c     error messages:
c***********************************************************************
 9980 write(iuout, *) 'props: error opening file ',filnam
      call ioerr(ioval)
      stop
 9990 write(iuout, *) 'props: error reading file ',filnam
      call ioerr(ioval)
      stop
 9991 write(iuout, *) 'props: unexpected end-of-file in ',filnam
      stop
      end
c
c***********************************************************************
c     brpair builds a bunch of molecule-size position arrays (r, r2,
c     1/r, dx, dy, dz) for a given pair of molecules.  It uses periodic
c     boundary conditions if applicable, and makes sure the same image
c     of each atom is used within a single molecule, using the first
c     atom as the reference....sjs 6/17/94   
c***********************************************************************
c     now uses dxinc, dyinc, dzinc passed as arguments.  call brhead
c     before calling this routine to set these.  This is so that
c     molecule separations can be obtained without filling the full
c     set of rmat arrays....sjs 4/12/95
c***********************************************************************
c     Try rewriting with BLAS routines
c***********************************************************************
c     pbrpr is the version for use in props, which calculates pair 
c     distances only for atoms with mass or fluctuating charge (and the
c     latter only if doing a gq calc)....sjs 6/12/95
c***********************************************************************
c
      subroutine pbrpr(imol1, imol2, dxinc, dyinc, dzinc, gqrflg)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c     
      logical gqrflg
c
      imty1 = molty(imol1)
      imty2 = molty(imol2)
      iatom1 = iatmol(imol1,1)
      iaty1 = iatype(imty1,1)
      do 200 imind2 = 2, numatm(imol2)
         iatom2 = iatmol(imol2,imind2)
         iaty2 = iatype(imty2,imind2)
         if ((isrl(iaty1) .and. isrl(iaty2)) .or. 
     $        (gqrflg .and. (isrl(iaty1) .and. isfqat(iaty2)) .or.
     $        (isfqat(iaty1) .and. isrl(iaty2)))) then
            dxmat(1,imind2) = pos(iatom2,1) - pos(iatom1,1) + dxinc
            dymat(1,imind2) = pos(iatom2,2) - pos(iatom1,2) + dyinc
            dzmat(1,imind2) = pos(iatom2,3) - pos(iatom1,3) + dzinc
            r2mat(1,imind2) = dxmat(1,imind2) * dxmat(1,imind2) +
     $           dymat(1,imind2) * dymat(1,imind2) +
     $           dzmat(1,imind2) * dzmat(1,imind2)
            rmat(1,imind2) = sqrt(r2mat(1,imind2))
            rinmat(1,imind2) = 1.0d0 / rmat(1,imind2)
         endif
  200 continue
      do 310 imind1 = 2, numatm(imol1)
         iatom1 = iatmol(imol1,imind1)
         iaty1 = iatype(imty1,imind1)
         do 300 imind2 = 1, numatm(imol2)
            iatom2 = iatmol(imol2,imind2)
            iaty2 = iatype(imty2,imind2)
            if ((isrl(iaty1) .and. isrl(iaty2)) .or. 
     $           (isrl(iaty1) .and. isfqat(iaty2)) .or.
     $           (isfqat(iaty1) .and. isrl(iaty2))) then
               dxmat(imind1,imind2) = pos(iatom2,1) - pos(iatom1,1) +
     $              dxinc
               dymat(imind1,imind2) = pos(iatom2,2) - pos(iatom1,2) +
     $              dyinc
               dzmat(imind1,imind2) = pos(iatom2,3) - pos(iatom1,3) +
     $              dzinc
               r2mat(imind1,imind2) = 
     $              dxmat(imind1,imind2) * dxmat(imind1,imind2) +
     $              dymat(imind1,imind2) * dymat(imind1,imind2) +
     $              dzmat(imind1,imind2) * dzmat(imind1,imind2)
               rmat(imind1,imind2) = sqrt(r2mat(imind1,imind2))
               rinmat(imind1,imind2) = 1.0d0 / rmat(imind1,imind2)
            endif
  300    continue
  310 continue
      return
      end
c
c***********************************************************************
c     foldr folds the pos() array of positions back into the main box,
c     storing them in fldpos().  All atoms in a given molecule must be 
c     in the same box, however....sjs 9/14/93
c***********************************************************************
c     pfoldr is like the foldr routine but it remembers how many box
c     lengths it had to fold each molecule back by, for use in the 
c     residence time calculations in props....sjs 4/27/95
c***********************************************************************
c
      subroutine pfoldr(imass, nfold)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      integer*4 nfold(maxmol,3)
c
      do 210 imol = 1, nmol
         imty = molty(imol)
         iaty = iatype(imty,1)
c$$$         if (mask(iaty)) then
         if (masmsk(iaty,imass)) then
            ioxy = iatmol(imol,1)
            nfold(imol,1) = anint((pos(ioxy,1) - boxby2(1)) / 
     $           boxsiz(1))
            nfold(imol,2) = anint((pos(ioxy,2) - boxby2(2)) / 
     $           boxsiz(2))
            nfold(imol,3) = anint((pos(ioxy,3) - boxby2(3)) / 
     $           boxsiz(3))
            dx = nfold(imol,1) * boxsiz(1)
            dy = nfold(imol,2) * boxsiz(2)
            dz = nfold(imol,3) * boxsiz(3)
c$$$            dx = anint((pos(ioxy,1) - boxby2(1)) / boxsiz(1)) * 
c$$$     $           boxsiz(1)
c$$$            dy = anint((pos(ioxy,2) - boxby2(2)) / boxsiz(2)) * 
c$$$     $           boxsiz(2)
c$$$            dz = anint((pos(ioxy,3) - boxby2(3)) / boxsiz(3)) * 
c$$$     $           boxsiz(3)
            fldpos(ioxy,1) = pos(ioxy,1) - dx
            fldpos(ioxy,2) = pos(ioxy,2) - dy
            fldpos(ioxy,3) = pos(ioxy,3) - dz
            do 200 iatinm = 2, numatm(imol)
               iatom = iatmol(imol,iatinm)
               fldpos(iatom,1) = pos(iatom,1) - dx
               fldpos(iatom,2) = pos(iatom,2) - dy
               fldpos(iatom,3) = pos(iatom,3) - dz
 200        continue
         endif
 210  continue
      return
      end
c         
c***********************************************************************
c     brhead figures out the separation between a pair of molecules,
c     using the first atom of each molecule as an approximation of its
c     center of mass (fine for water and ions).  The increments for how
c     many box lengths should be subtracted to wrap the pair into the
c     same image cell are passed back, for use in brpair....sjs 4/12/95
c***********************************************************************
c     pbrhd is the version for props, which sends back the number of
c     box lengths apart the two molecules are....sjs 4/27/95
c***********************************************************************
c
      subroutine pbrhd(imol1, imol2, dxinc, dyinc, dzinc, nfold)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      integer*4 nfold(3)
c
      iatom1 = iatmol(imol1,1)
      iatom2 = iatmol(imol2,1)
      dxmat(1,1) = pos(iatom2,1) - pos(iatom1,1)
      dymat(1,1) = pos(iatom2,2) - pos(iatom1,2)
      dzmat(1,1) = pos(iatom2,3) - pos(iatom1,3)
      if (pbflag) then
         nfold(1) = -1 * nint(dxmat(1,1) / boxsiz(1))
         nfold(2) = -1 * nint(dymat(1,1) / boxsiz(2))
         nfold(3) = -1 * nint(dzmat(1,1) / boxsiz(3))
         dxinc = dble(nfold(1)) * boxsiz(1)
         dyinc = dble(nfold(2)) * boxsiz(2)
         dzinc = dble(nfold(3)) * boxsiz(3)
c$$$         dxinc = -1 * dnint(dxmat(1,1) / boxsiz(1)) * boxsiz(1)
c$$$         dyinc = -1 * dnint(dymat(1,1) / boxsiz(2)) * boxsiz(2)
c$$$         dzinc = -1 * dnint(dzmat(1,1) / boxsiz(3)) * boxsiz(3)
      else
         nfold(1) = 0
         nfold(2) = 0
         nfold(3) = 0
         dxinc = 0.d0
         dyinc = 0.d0
         dzinc = 0.d0
      endif
      dxmat(1,1) = dxmat(1,1) + dxinc
      dymat(1,1) = dymat(1,1) + dyinc
      dzmat(1,1) = dzmat(1,1) + dzinc
      r2mat(1,1) = dxmat(1,1) * dxmat(1,1) +
     $     dymat(1,1) * dymat(1,1) +
     $     dzmat(1,1) * dzmat(1,1)
      rmat(1,1) = sqrt(r2mat(1,1))
      rinmat(1,1) = 1.0d0 / rmat(1,1)
      return
      end
c
c***********************************************************************
c     wrteor writes out an end-of-run file, either at the end of the 
c     run, or earlier, in case of a crash....sjs 10/5/93
c***********************************************************************
c     wrtsnp is a version for props that writes out at least some of
c     the eor info to a snapshot file....sjs 7/21/95
c***********************************************************************
c
      subroutine wrtsnp()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      write(filnam, '(3a)') datdir(1:index(datdir, ' ')-1), 
     $     '/snap.', suffix
      open(iutmp, file = filnam, err = 9990, iostat = ioval)
      write(iutmp, *, err = 9980, iostat = ioval)
     $     iosysv
      write(iutmp, *, err = 9980, iostat = ioval)
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
            write(iutmp, *, err = 9980, iostat = ioval)
     $           q(iatom), qvel(iatom)
 290     continue
         do 291 itind = 1, ntag(imty)
            imind = itind + nframe(imty)
            iatom = iatmol(imol,imind)
            write(iutmp, *, err = 9980, iostat = ioval)
     $           itind, ident(iatom)
            write(iutmp, *, err = 9980, iostat = ioval)
     $           q(iatom), qvel(iatom)
 291     continue
 300  continue
      write(iutmp, *, err = 9980, iostat = ioval) iocbnv
      write(iutmp, *, err = 9980, iostat = ioval) niters
      write(iutmp, *, err = 9980, iostat = ioval) eqflag
      write(iutmp, *, err = 9980, iostat = ioval) nscbas
      write(iutmp, *, err = 9980, iostat = ioval) wvirnx
      write(iutmp, *, err = 9980, iostat = ioval) maxbin
      write(iutmp, *, err = 9980, iostat = ioval) rktsum, rktssq
      close(iutmp)
      return
 9980 write(iuout, *) 'wrtsnp: error writing to ', filnam
      call ioerr(ioval)
      stop
 9990 write(iuout, *) 'wrtsnp: error opening ', filnam
      call ioerr(ioval)
      end
c
