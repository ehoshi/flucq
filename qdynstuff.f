c***********************************************************************
c     qdinit initializes whatever can be initialized in a subroutine
c     (the stuff in common blocks)....sjs 6/22/93
c***********************************************************************
c
      subroutine qdinit
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      movlt = .true.
      movhvy = .true.
      nerbin = 10000
      niters = 0
      nscbas = 0
      do 110 iatom = 1, maxatm
         hasmaa(iatom) = .false.
         coskx(iatom,1,1) = 1.
         coskx(iatom,1,2) = 1.
         coskx(iatom,1,3) = 1.
         sinkx(iatom,1,1) = 0.
         sinkx(iatom,1,2) = 0.
         sinkx(iatom,1,3) = 0.
  110 continue
      return
      end
c
c***********************************************************************
c     itinit initializes all the common-block stuff that needs to be
c     initialized (or incremented, evaluated, ...) every iteration.
c     ...sjs 9/7/93
c***********************************************************************
c
      subroutine itinit()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
c***********************************************************************
c     Resolve for the charges or Drude positions whenever it has been 
c     too long since the last time they were solved for, or at the 
c     beginning of the run:
c***********************************************************************
      if (qslvfl .and. (niters - lqslv) .ge. nqsolv) then
         qwrong = .true.
         lqslv = niters
      endif
      if (dslvfl .and. (niters - ldslv) .ge. ndsolv) then
         dwrong = .true.
         ldslv = niters
      endif
c***********************************************************************
c     If charges or Drude positions are screwed up, solve for them
c     (turning off atom motion during the resolving).  If they're both
c     screwed up, solve for each of them alternately until they're both
c     okay.
c***********************************************************************
      if (qwrong .and. dwrong) then
         fsiter = .false.
         if (qsiter) then
            if (dslvfl) then
               dsiter = .true.
            endif
            qsiter = .false.
         else 
            dsiter = .false.
            qsiter = .true.
         endif
         moving = .false.
      else if (dwrong) then
         fsiter = .false.
         dsiter = .true.
         qsiter = .false.
         moving = .false.
      else if (qwrong) then
         fsiter = .false.
         dsiter = .false.
         qsiter = .true.
         moving = .false.
      else
         dsiter = .false.
         qsiter = .false.
         if (fsiter .or. moving) then
            fsiter = .false.
            moving = .true.
         else
            fsiter = .true.
            moving = .false.
            write(iuout, *) 'doing a dry run for forces...'
         endif
      endif
c***********************************************************************
c     The time gets incremented only when dynamics is being done:
c***********************************************************************
      if (moving) then
         niters = niters + 1
      endif
      ake = 0.d0
      dke = 0.d0
      qke = 0.d0
      if (ioflg .and. nioint * (niters / nioint) .eq. niters) then
         ioiter = .true.
      else
         ioiter = .false.
      endif
      if (iplint * (niters / iplint) .eq. niters) then
         pliter = .true.
      else
         if (moving) then
            pliter = .false.
         else
            pliter = .true.
         endif
      endif
cmov      if (dbflag) then
cmov         do 160 ietype = 1, 9
cmov            ebit(ietype) = 0.0d0
cmov  160    continue
cmov      endif
cmov      if (dbflag) ebit(9) = ebit(9) - basepe
      return
      end
c
c***********************************************************************
c     Block data subroutine required for the picky IBMs.  Data needed
c     for the dynamics are initialized here...sjs 1/19/93
c***********************************************************************
c
      block data qdydat
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
c***********************************************************************
c     These are both the absolute and the expected upper and lower
c     limits for the partial charges on atoms of the specified atomic
c     number.
c***********************************************************************
      data qhigh(1), qhigh(6), qhigh(7), qhigh(8), qhigh(9), qhigh(17)
     $     /1.0, 4.0, 5.0, 6.0, 7.0, 7.0/
      data qlow(1), qlow(6), qlow(7), qlow(8), qlow(9), qlow(17)
     $     /-1.0, -4.0, -3.0, -2.0, -1.0, -1.0/
      data qlowg(1), qlowg(6), qlowg(7), qlowg(8), qlowg(9), qlowg(17)
     $     /0.2, -1.0, -3.0, -1.6, -1.0, -1.0/
      data qhighg(1), qhighg(6), qhighg(7), qhighg(8), qhighg(9),
     $     qhighg(17)
     $     /0.9, 1.0, 4.0, -0.5, -0.5,
     $     -0.5/
      end
c
c***********************************************************************
c     rmsdev calculates the rms deviation between the two arrays passed
c     to it.  The number of elements to be used must also be passed.
c     ...sjs 3/3/93
c***********************************************************************
c
      real*8 function rmsdev(a, b, n)
c
      include 'implic'
c
      real*8 a(n), b(n)
c
      rmsdev = 0.0d0
      do 200 i = 1, n
         rmsdev = rmsdev + (a(i) - b(i)) * (a(i) - b(i))
  200 continue
      rmsdev = rmsdev / n
      rmsdev = dsqrt(rmsdev)
      return
      end
c
c***********************************************************************
c     siding is a do-nothing subroutine that provides a convenient stop
c     for a debugging flag....sjs 6/9/93
c***********************************************************************
c
      subroutine siding
c
      return
      end
c
c***********************************************************************
c     qsolv2 is similar to qsolv, except that it rearranges the qmat
c     matrix to account for the charge constraints.  If this has already
c     been done, call qsolv instead.  qmat should already be symmetric, 
c     however....sjs 7/6/93
c***********************************************************************
c     If there are rigid charges around, the chi vector must be
c     calculated before the qmat array is destroyed....sjs 2/1/95
c***********************************************************************
c
      subroutine qsolv2()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      integer ipvt(maxatm)
      real*8 qmatsm(maxatm,maxatm),
     $     chivec(maxatm), chibig(maxatm)
c
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
c***********************************************************************
c     Fill the RHS vector with chi_N - chi_i + [(J_Nj - J_ij) * q_j]
c     (summed over rigid j):
c***********************************************************************
      if (.not. ctflag) then
         do 230 imol1 = 1, nmol
            imty1 = molty(imol1)
            if (hasfqt(imty1)) then
               iendi = ifqmol(imty1,nfqmol(imty1))
               iend = iatmol(imol1,iendi)
               tchi = trmq(iatype(imty1,iendi),1)
               do 220 iqind1 = 1, nfqmol(imty1) - 1
                  iatom1 = iatmol(imol1,ifqmol(imty1,iqind1))
                  fixqtm = 0.d0
                  do 210 imol2 = 1, nmol
                     imty2 = molty(imol2)
                     if (hasrqt(imty2)) then
                        do 200 iqind2 = 1, nrqmol(imty2)
                           iatom2 = iatmol(imol2,irqmol(imty2,iqind2))
                           fixqtm = fixqtm + (qmat(iend,iatom2) - 
     $                          qmat(iatom1,iatom2)) * q(iatom2)
 200                    continue
                     endif
 210              continue
                  chibig(iatom1) = tchi -
     $                 trmq(iatype(imty1,ifqmol(imty1,iqind1)),1) +
     $                 fixqtm
 220           continue
            endif
 230     continue
      else
         do 240 imol = nmol, 1, -1
            imty = molty(imol)
            if (hasfqt(imty)) then
               iendi = ifqmol(imty,nfqmol(imty))
               iend = iatmol(imol,iendi)
               tchi = trmq(iatype(imty,iendi),1)
               go to 250
            endif
 240     continue
 250     continue
         do 290 imol1 = 1, nmol
            imty1 = molty(imol1)
            do 280 iqind1 = 1, nfqmol(imty1)
               iatom1 = iatmol(imol1,ifqmol(imty1,iqind1))
               fixqtm = 0.d0
               do 270 imol2 = 1, nmol
                  imty2 = molty(imol2)
                  do 260 iqind2 = 1, nrqmol(imty2)
                     iatom2 = iatmol(imol2,irqmol(imty2,iqind2))
                     fixqtm = fixqtm + (qmat(iend,iatom2) -
     $                    qmat(iatom1,iatom2)) * q(iatom2)
 260              continue
 270           continue
               chibig(iatom1) = tchi - 
     $              trmq(iatype(imty1,ifqmol(imty1,iqind1)),1) + fixqtm
 280        continue
 290     continue
      endif
c***********************************************************************
c     Do the proper fiddling to take care of the charge constraints:
c***********************************************************************
      if (.not. ctflag) then
         do 330 imol1 = 1, nmol
            imty1 = molty(imol1)
            if (hasfqt(imty1)) then
               iend1 = iatmol(imol1,ifqmol(imty1,nfqmol(imty1)))
               do 320 imol2 = 1, nmol
                  imty2 = molty(imol2)
                  if (hasfqt(imty2)) then
                     iend2 = iatmol(imol2,ifqmol(imty2,nfqmol(imty2)))
                     do 310 iqind1 = 1, nfqmol(imty1) - 1
                        iatom1 = iatmol(imol1,ifqmol(imty1,iqind1))
                        do 300 iqind2 = 1, nfqmol(imty2) - 1
                           iatom2 = iatmol(imol2,ifqmol(imty2,iqind2))
                           qmat(iatom1,iatom2) = qmat(iatom1,iatom2) -
     $                          qmat(iatom1,iend2) - 
     $                          qmat(iend1,iatom2) + qmat(iend1,iend2)
 300                    continue
 310                 continue
                  endif
 320           continue
            endif
 330     continue
      else
         do 340 imol = nmol, 1, -1
            imty = molty(imol)
            if (hasfqt(imty)) then
               iend = iatmol(imol,ifqmol(imty,nfqmol(imty)))
               go to 350
            endif
 340     continue
 350     continue
         do 390 imol1 = 1, nmol
            imty1 = molty(imol1)
            do 380 imol2 = 1, nmol
               imty2 = molty(imol2)
               do 370 iqind1 = 1, nfqmol(imty1)
                  iatom1 = iatmol(imol1,ifqmol(imty1,iqind1))
                  do 360 iqind2 = 1, nfqmol(imty2)
                     iatom2 = iatmol(imol2,ifqmol(imty2,iqind2))
                     qmat(iatom1,iatom2) = qmat(iatom1,iatom2) -
     $                    qmat(iatom1,iend) - qmat(iend,iatom2) +
     $                    qmat(iend,iend)
 360              continue
 370           continue
 380        continue
 390     continue
      endif
c***********************************************************************
c     Collapse the big qmat array and chi vector down into an array and
c     vector for only independent charges:
c***********************************************************************
      do 410 iqind1 = 1, nqind
         chivec(iqind1) = chibig(indbig(iqind1))
         do 400 iqind2 = 1, nqind
            qmatsm(iqind1,iqind2) = qmat(indbig(iqind1),indbig(iqind2))
  400    continue
  410 continue
c$$$c***********************************************************************
c$$$c     Solve for the qs!  qmatsm is destroyed in dgef (Double-precision
c$$$c     Gaussian Elimination Factorizer), and chivec is
c$$$c     overwritten with the solution in dges (Double-precision Gaussian
c$$$c     Elimination Solver).  These are ESSL routines.
c$$$c***********************************************************************
c$$$      iopt = 0
c$$$      call dgef(qmatsm, maxatm, nqind, ipvt)
c$$$      call dges(qmatsm, maxatm, nqind, ipvt, chivec, iopt)
c***********************************************************************
c     Solve for the qs!  dgesv (Double-precision Gaussian Elimination
c     Solver) is a LAPACK routine.  qmatsm is destroyed, and chivec
c     is overwritten with the solution.
c***********************************************************************
      call dgesv(nqind, 1, qmatsm, maxatm, ipvt, chivec, maxatm, iexval)
      if (iexval .ne. 0) then
         write(iuout, *) 'qsolv2:  dgesv failed'
         stop
      endif
c***********************************************************************
c     Fill up the qans vector with the answers, solving for the 
c     dependent charges:
c***********************************************************************
      if (.not. ctflag) then
         do 610 imol = 1, nmol
            imty = molty(imol)
            if (hasfqt(imty)) then
               qsum = 0.0d0
               do 600 iqind = 1, nfqmol(imty) - 1
                  iatom = iatmol(imol,ifqmol(imty,iqind))
                  qans(iatom) = chivec(indsml(iatom))
                  qsum = qsum + qans(iatom)
 600           continue
               qans(iatmol(imol,ifqmol(imty,nfqmol(imty)))) = -qsum
            endif
 610     continue
      else
         qsum = 0.0d0
         do 620 iqind = 1, nqind
            iatom = indbig(iqind)
            qans(iatom) = chivec(iqind)
            qsum = qsum + chivec(iqind)
 620     continue
         qans(indbig(nqind + 1)) = -qsum
      endif
c***********************************************************************
c     The fluctuating charges (right now only on heavy atoms) have been
c     changed, therefore these atoms have "moved" (the forces between
c     them and others will have changed):
c***********************************************************************
      movhvy = .true.
      return
      end
c
c***********************************************************************
c     qcheck checks to see how many charges are outside of the specified
c     limits for their atom types, and complains about it if any are.
c     ...sjs 9/7/93
c***********************************************************************
c
      subroutine qcheck()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      nqout = 0
      do 205 imol = 1, nmol
         imty = molty(imol)
         do 200 iqind = 1, nfqmol(imty)
            iatom = iatmol(imol,ifqmol(imty,iqind))
            iatid = ident(iatom)
            if (q(iatom) .gt. qhigh(iatid) .or. 
     $           q(iatom) .lt. qlow(iatid)) then
               nqout = nqout + 1
               write(iuout, *) iatom, q(iatom)
            endif
 200     continue
 205  continue
      if (nqout .ne. 0)
     $     write(iuout, '(a, i6, a, i4, a)') 'Timestep ', niters, ':',
     $     nqout, ' atoms have illegal charges'
      return
      end
c
c***********************************************************************
c     dipols calculates the molecular and system dipoles....sjs 9/7/93
c***********************************************************************
c     updated to use an atom mask, putting partial molecular and system
c     dipoles into the arguments passed to it, but still putting the
c     full molecular and system dipoles into rmu and smu....sjs 3/15/95
c***********************************************************************
c
      subroutine dipols(getlt, gethvy)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      logical getlt, gethvy
c
      if (getlt .and. movlt) then
         smu(1,ilight) = 0.d0
         smu(2,ilight) = 0.d0
         smu(3,ilight) = 0.d0
      endif
      if ((getlt .or. gethvy) .and. (movlt .or. movhvy)) then
         smu(1,imixed) = 0.d0
         smu(2,imixed) = 0.d0
         smu(3,imixed) = 0.d0
      endif
      if (gethvy .and. movhvy) then
         smu(1,iheavy) = 0.d0
         smu(2,iheavy) = 0.d0
         smu(3,iheavy) = 0.d0
      endif
      do 210 imol = 1, nmol
         do 200 iatno = 1, numatm(imol)
            iaty = iatype(molty(imol),iatno)
            iatom = iatmol(imol,iatno)
            rmupcx = q(iatom) * fldpos(iatom,1)
            rmupcy = q(iatom) * fldpos(iatom,2)
            rmupcz = q(iatom) * fldpos(iatom,3)
            if (masmsk(iaty,ilight) .and. getlt .and. movlt) then
               smu(1,ilight) = smu(1,ilight) + rmupcx
               smu(2,ilight) = smu(2,ilight) + rmupcy
               smu(3,ilight) = smu(3,ilight) + rmupcz
            endif
            if ((getlt .or. gethvy) .and. (movlt .or. movhvy)) then
               smu(1,imixed) = smu(1,imixed) + rmupcx
               smu(2,imixed) = smu(2,imixed) + rmupcy
               smu(3,imixed) = smu(3,imixed) + rmupcz
            endif
            if (masmsk(iaty,iheavy) .and. gethvy .and. movhvy) then
               smu(1,iheavy) = smu(1,iheavy) + rmupcx
               smu(2,iheavy) = smu(2,iheavy) + rmupcy
               smu(3,iheavy) = smu(3,iheavy) + rmupcz
            endif
  200    continue
  210 continue
      return
      end
c
cmovc***********************************************************************
cmovc     Stuff left over from bnupdt:
cmovc***********************************************************************
cmovc     lower-diagonalize qbpre, qepre, and uqbpre, and remove the extra 
cmovc     factor of 2:
cmovc***********************************************************************
cmov      if (iobwfl) then
cmov         do 290 imol1 = 1, nmol
cmov            qbpre(imol1,imol1) = 0.5d0 * qbpre(imol1,imol1)
cmov            qepre(imol1,imol1) = 0.5d0 * qepre(imol1,imol1)
cmov            uqbpre(imol1,imol1) = 0.5d0 * uqbpre(imol1,imol1)
cmov            selfe(imol1) = 0.5d0 * selfe(imol1)
cmov            do 280 imol2 = 1, imol1 - 1
cmov               qbpre(imol1,imol2) = qbpre(imol1,imol2) +
cmov     $              qbpre(imol2,imol1)
cmov               qepre(imol1,imol2) = qepre(imol1,imol2) +
cmov     $              qepre(imol2,imol1)
cmov               uqbpre(imol1,imol2) = uqbpre(imol1,imol2) +
cmov     $              uqbpre(imol1,imol2)
cmov               qbpre(imol1,imol2) = qbpre(imol1,imol2) * 0.5d0
cmov               qepre(imol1,imol2) = qepre(imol1,imol2) * 0.5d0
cmov               uqbpre(imol1,imol2) = uqbpre(imol1,imol2) * 0.5d0
cmov 280        continue
cmov 290     continue
cmov         do 291 imol = 1, nmol
cmov            selfe(imol) = selfe(imol) + qbpre(imol,imol)
cmov            qepre(imol,imol) = qepre(imol,imol) - qbpre(imol,imol)
cmov 291     continue
cmovc***********************************************************************
cmovc     accumulate the denominators for weighting the self energies (these
cmovc     are the total charge-only interaction between one molecule and all
cmovc     other molecules, including all Ewald images but not the original
cmovc     molecule):
cmovc***********************************************************************
cmov         do 293 imol1 = 1, nmol
cmov            qbinde(imol1) = qbinde(imol1) + qepre(imol1,imol1)
cmov            do 292 imol2 = 1, imol1 - 1
cmov               qbinde(imol1) = qbinde(imol1) + qepre(imol1,imol2)
cmov               qbinde(imol2) = qbinde(imol2) + qepre(imol1,imol2)
cmov 292        continue
cmov 293     continue
cmovc***********************************************************************
cmovc     bin the normalized dimer energy, which is the total (q- plus 
cmovc     non-q) bare pair energy plus a piece of the self energy 
cmovc     proportional to the charge-only bare pair energy.  also count
cmovc     H-bonds:
cmovc***********************************************************************
cmov         do 310 imol1 = 1, nmol
cmov            imty1 = molty(imol1)
cmov            do 300 imol2 = 1, imol1 - 1
cmov               imty2 = molty(imol2)
cmov               dimenr = qbpre(imol1,imol2) + uqbpre(imol1,imol2) +
cmov     $              qbpre(imol1,imol2) * 
cmov     $              (selfe(imol1) / qbinde(imol1) +
cmov     $              selfe(imol2) / qbinde(imol2))
cmov               ibin = int((dimenr - prEbmn) * prEbwi) + 1
cmov               if (ibin .ge. 1 .and. ibin .le. maxbin)
cmov     $              npairE(ibin) = npairE(ibin) + 1
cmov               if (dimenr .lt. hbonde) then
cmov                  nhbond(imol1) = nhbond(imol1) + 1
cmov                  nhbond(imol2) = nhbond(imol2) + 1
cmov               endif
cmov 300        continue
cmov 310     continue
cmovc***********************************************************************
cmovc     bin the binding energies, which are a sum of full, Ewald pair
cmovc     energies, plus the self energy for the molecule in question, plus
cmovc     bits of the self energy (proportional to the Ewald charge-only
cmovc     interaction energy) for each "other half":
cmovc***********************************************************************
cmov         sbnde = 0.0
cmov         bindel = bigpos
cmov         bindeh = bigneg
cmov         do 318 imol1 = 1, nmol
cmov            binde = selfe(imol1)
cmov            do 313 imol2 = 1, imol1 - 1
cmov               binde = binde + qepre(imol1,imol2) + uqbpre(imol1,imol2)
cmov               binde = binde + qepre(imol1,imol2) / qbinde(imol2) *
cmov     $              selfe(imol2)
cmov 313        continue
cmov            binde = binde + qepre(imol1,imol1)
cmov            binde = binde + qepre(imol1,imol1) / qbinde(imol1) *
cmov     $           selfe(imol1)
cmov            do 315 imol2 = imol1 + 1, nmol
cmov               binde = binde + qepre(imol2,imol1) + uqbpre(imol2,imol1)
cmov               binde = binde + qepre(imol2,imol1) / qbinde(imol2) *
cmov     $              selfe(imol2)
cmov 315        continue
cmov            ibin = int((binde - bdEbmn) * bdEbwi) + 1
cmov            if (ibin .ge. 1 .and. ibin .le. maxbin) then
cmov               nbonde(ibin) = nbonde(ibin) + 1
cmov            endif
cmov            if (binde .lt. bindel) then
cmov               bindel = binde
cmov            endif
cmov            if (binde .gt. bindeh) then
cmov               bindeh = binde
cmov            endif
cmov            sbnde = sbnde + binde
cmov 318     continue
c
c***********************************************************************
c     conchk checks to make sure that all of the molecule frames meet
c     their distance constraints....sjs 6/14/94
c***********************************************************************
c
      subroutine conchk
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      do 210 imol = 1, nmol
         imty = molty(imol)
         do 200 iconi = 1, ncons(imty)
            iatom1 = iatmol(imol,icons(imty,iconi,1))
            iatom2 = iatmol(imol,icons(imty,iconi,2))
            dx = pos(iatom1,1) - pos(iatom2,1)
            dy = pos(iatom1,2) - pos(iatom2,2)
            dz = pos(iatom1,3) - pos(iatom2,3)
            tr2dif = dx * dx + dy * dy + dz * dz - r2cons(imty,iconi)
            if (abs(tr2dif) .gt. rattol) then
               write(iuout, '(a, i5, a, e10.3, a)') 
     $              'conchk: WARNING! molecule ',
     $              imol, ' fails a bond constraint by ', 
     $              sqrt(tr2dif + r2cons(imty,iconi)) - 
     $              sqrt(r2cons(imty,iconi)), ' A'
            endif
 200     continue
 210  continue
      return
      end
c
c***********************************************************************
c     sterfc sets up the erfca array with values of erfc in the
c     specified range, with the specified number of bins.  These values
c     are later read and interpolated by berfc....sjs 7/19/94
c***********************************************************************
c
      subroutine sterfc()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
c***********************************************************************
c     erfc is meaningless for args below zero:
c***********************************************************************
      if (erfclo .lt. 0.0d0) erfclo = 0.0d0
c***********************************************************************
c     on the IBMs, erfc is zero for args at or above 10.0
c***********************************************************************
      if (erfchi .gt. 10.0d0) erfchi = 10.0d0
c***********************************************************************
c     these will be used in berfc:
c***********************************************************************
      stpsiz = (erfchi - erfclo) / (nerbin - 1)
      stpb2 = stpsiz / 2.0d0
      stp2 = stpsiz * stpsiz
      stp2i = 1.0d0 / stp2
      stp22i = 0.5d0 * stp2i
c***********************************************************************
c     calculate the erfc values and stick them into an array:
c***********************************************************************
      do 100 ibin = 1, nerbin
         valarg = erfclo + stpsiz * (ibin - 1)
         erfca(ibin) = erfc(valarg)
 100  continue
      return
      end
c
c***********************************************************************
c     berfc gets an approximate value of erfc for arg, by
c     interpolating from the erfca array, which was set in sterfc.  The
c     quadratic interpolation is taken from Numerical Recipes.
c     ...sjs 7/19/94
c***********************************************************************
c
      real*8 function berfc(arg)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
c***********************************************************************
c     On the IBMs, erfc is zero for args at or above 10.0:
c***********************************************************************
      if (arg .gt. 10.0d0) then
         berfc = 0.0d0
         return
      endif
c***********************************************************************
c     Find the closest point for the quadratic fit:
c***********************************************************************
      ibinmd = int((arg + stpb2 - erfclo) / stpsiz)
c***********************************************************************
c     If it's at an endpoint, move it in:
c     (even if it's outside the endpoints, which probably results from
c     illegal input values (like NaN), set the bin to a value that
c     won't cause segfaults)
c***********************************************************************
      if (ibinmd .le. 0) then
         ibinmd = 1
      else if (ibinmd .ge. nerbin - 1) then
         ibinmd = nerbin - 2
      endif
c***********************************************************************
c     Calculate the interpolated value:
c***********************************************************************
      argmd = erfclo + ibinmd * stpsiz
      ibinmd = ibinmd + 1
      dmd = arg - argmd
      dmd2 = dmd * dmd
      stpdmd = stpsiz * dmd
      berfc = (dmd2 - stpdmd) * stp22i * erfca(ibinmd - 1) +
     $     (stp2 - dmd2) * stp2i * erfca(ibinmd) +
     $     (dmd2 + stpdmd) * stp22i * erfca(ibinmd + 1)
      return
      end
c
c***********************************************************************
c     setsys sets a bunch of variables related to the variables read in
c     from the system file (*.sys, *.eor)....sjs 7/28/94
c***********************************************************************
c
      subroutine setsys(anaflg)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      logical anaflg
c
      iabyte = (2 * natoms + 3) * nspbyt
      idbyte = (8 * natoms + 8) * nspbyt
      if (boxsiz(1) .eq. boxsiz(2) .and. boxsiz(2) .eq. boxsiz(3)) then
         cubflg = .true.
         boxby2(1) = boxsiz(1) / 2.0d0
         boxby2(2) = boxby2(1)
         boxby2(3) = boxby2(1)
         hbox = boxby2(1)
      else
         cubflg = .false.
         hbox = boxsiz(1)
         bbox = boxsiz(1)
         if (boxsiz(2) .lt. hbox) hbox = boxsiz(2)
         if (boxsiz(2) .gt. bbox) bbox = boxsiz(2)
         if (boxsiz(3) .lt. hbox) hbox = boxsiz(3)
         if (boxsiz(3) .gt. bbox) bbox = boxsiz(3)
         hbox = hbox / 2.0d0
         boxby2(1) = boxsiz(1) / 2.0d0
         boxby2(2) = boxsiz(2) / 2.0d0
         boxby2(3) = boxsiz(3) / 2.0d0
      endif
      erfclo = 0.
      if (cubflg) then
         erfchi = ewlkap * hbox * sqrt(3.)
      else
         erfchi = ewlkap * bbox * sqrt(3.) / 2.
      endif
      call sterfc()
      if (ewlflg) then
         if (cubflg) then
            ewlkap = ewlkap / boxsiz(1)
         else
            temp = 2 * hbox
            write(iuout, *) 'rdsim: WARNING! box is not cubic - using ',
     $           'the smallest side (', temp, ' A) to scale ewlkap.'
            ewlkap = ewlkap / temp
         endif
         r2kbrp = 2.0d0 * ewlkap / sqrt(pi)
         ewlkp2 = ewlkap * ewlkap
      endif
      vol = boxsiz(1) * boxsiz(2) * boxsiz(3)
      if (ewsrfl) srffac = epsinv * 2.0 * 2.0 * pi / 3.0 / vol
      ntstep = 0
      rktsum = 0.d0
      rktssq = 0.d0
      nmty = 0
      do 440 imid = 1, maxmid
         if (hasmid(imid)) then
            nmty = nmty + 1
            if (nmty .gt. maxmty) then
               write(iuout, *) 
     $              'rdsim: too many molecule types in system'
               stop
            endif
            mtymid(imid) = nmty
            midmty(nmty) = imid
         endif
 440     continue
      do 450 imol = 1, nmol
         molty(imol) = mtymid(molid(imol))
 450  continue
      return
      end
c
c***********************************************************************
c     chksys checks the variables set in reading the system file (*.sys,
c     *.eor) for illegal values and incompatibilities with other
c     variables....sjs 7/28/94
c***********************************************************************
c
      subroutine chksys(icjid)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      integer*4 icjid(mxcstj,3)
c
      if (iabyte / nspbyt .gt. maxdat) then
         write(iuout, *) 
     $        'chksys: Too many data items written to anneal.', suffix
         stop
      endif
      if (idbyte / nspbyt .gt. maxdat) then
         write(iuout, *) 'chksys: Too many data items will be ',
     $        'written to ', 'dyn.', suffix
         stop
      endif
      if (nmol .eq. 1 .and. vscflg) then
         write(iuout, *) 'chksys: WARNING! no v-scaling for 1 molecule'
         vscflg = .false.
      endif
      if (hbox .lt. 0) stop 'chksys: box size must be positive!'
      do 200 icint = 1, ncint
         imid1 = icjid(icint,1)
         if (.not. hasmid(imid1)) then
            write(iuout, *) 'chksys:  can''t use a custom J(r) for ',
     $           'molecule id ', imid1
            stop
         endif
         imid2 = icjid(icint,2)
         if (.not. hasmid(imid2)) then
            write(iuout, *) 'chksys:  can''t use a custom J(r) for ',
     $           'molecule id ', imid2
            stop
         endif
  200 continue
      return
      end
c
c***********************************************************************
c     trmjst fills trmj for each molecule type....sjs 8/1/94
c     and trmdj....sjs 12/9/94
c***********************************************************************
c
      subroutine trmjst()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      do 240 imty = 1, nmty
         do 190 imind1 = 1, natmty(imty)
            do 180 imind2 = 1, natmty(imty)
               trmj(imty,imind1,imind2) = 0.d0
               trmdj(imty,imind1,imind2) = 0.d0
 180        continue
 190     continue
         if (hasqt(imty)) then
            do 200 imol = 1, nmol
               if (molty(imol) .eq. imty) go to 210
 200        continue
            write(iuout, *) 'trmjst: can''t find molecule with id ', 
     $           midmty(imty)
            stop
 210        continue
            call brmono(imol)
            do 230 iqind1 = 1, nqmol(imty)
               imind1 = iqmol(imty,iqind1)
               iatyp1 = iatype(imty,imind1)
               zero = 0.0
               trmj(imty,imind1,imind1) = quadj(iatyp1, iatyp1, zero)
               do 220 iqind2 = 1, iqind1 - 1
                  imind2 = iqmol(imty,iqind2)
                  iatyp2 = iatype(imty,imind2)
                  trmj(imty,imind1,imind2) = 
     $                 quadj(iatyp1, iatyp2, rmat(imind1,imind2))
                  trmj(imty,imind2,imind1) = trmj(imty,imind1,imind2)
                  trmdj(imty,imind1,imind2) = 
     $                 quaddj(iatyp1, iatyp2, rmat(imind1,imind2))
                  trmdj(imty,imind2,imind1) = trmdj(imty,imind1,imind2)
 220           continue
 230        continue
         endif
 240  continue
      return
      end
c
c***********************************************************************
c     setmod sets a bunch of variables that need to be set after
c     reading in the model information....sjs 8/3/94
c***********************************************************************
c
      subroutine setmod(icjid)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      integer*4 icjid(mxcstj,3)
c
c***********************************************************************
c     initialize some stuff:
c***********************************************************************
      nfqaty = 0
      nmaaty = 0
      do 100 iaty = 1, maxaty
         isdrsp(iaty) = .false.
         masmsk(iaty,ilight) = .false.
         masmsk(iaty,iheavy) = .true.
         masmsk(iaty,imixed) = .true.
 100  continue
c***********************************************************************
c     count the number of fluc-q atom types (nfqaty) and the number of
c     atom types with mass (nmaaty).
c***********************************************************************
      do 200 iaty = 1, naty
         nwaty(iaty) = 0
         if (isfqat(iaty)) then
            nfqaty = nfqaty + 1
         endif
         if (ismasv(iaty)) then
            nmaaty = nmaaty + 1
         endif
 200  continue
c***********************************************************************
c     set up the atom->atom type array (itype) and count the number of 
c     atoms with each atom type (nwaty).
c***********************************************************************
      do 310 imol = 1, nmol
         imty = molty(imol)
         do 300 imind = 1, numatm(imol)
            iatom = iatmol(imol,imind)
            iaty = iatype(imty,imind)
            itype(iatom) = iaty
            nwaty(iaty) = nwaty(iaty) + 1
 300     continue
 310  continue
c***********************************************************************
c     set the flag for atoms with mass (hasmaa):
c***********************************************************************
      do 410 imol = 1, nmol
         imty = molty(imol)
         do 400 imaind = 1, nmamol(imty)
            iatom = iatmol(imol,imamol(imty,imaind))
            hasmaa(iatom) = .true.
 400     continue
 410  continue
c***********************************************************************
c     set the flags for Drude spring atoms (isdrsp, isdrat).  Also set
c     mass masks to label both base and satellite particles as light,
c     regardless of their physical mass:
c***********************************************************************
      do 520 imty = 1, nmty
         do 500 idrid = 1, ndrmol(imty)
            iaty = iatype(imty,idrmol(imty,idrid,2)) 
            isdrsp(iaty) = .true.
            isdrat(iaty) = .true.
            masmsk(iaty,ilight) = .true.
            masmsk(iaty,iheavy) = .false.
            iaty = iatype(imty,idrmol(imty,idrid,1))
            isdrat(iaty) = .true.
            masmsk(iaty,ilight) = .true.
            masmsk(iaty,iheavy) = .false.
 500     continue
 520  continue
c***********************************************************************
c     set the mass masks for all other atoms below the mass cutoff
c     which are not completely massless:
c     (only tag atoms are massless currently, and tag atoms can't be
c     light)
c***********************************************************************
      do 535 iaty = 1, naty
         if (atmass(iaty) .lt. hvywt .and. atmass(iaty) .gt. 0.d0) then
            masmsk(iaty,ilight) = .true.
            masmsk(iaty,iheavy) = .false.
         endif
  535 continue
c***********************************************************************
c     set the flags for whether molecule types have heavy and light 
c     atoms:
c***********************************************************************
      do 540 imty = 1, nmty
         hashvy(imty) = .false.
         haslt(imty) = .false.
         do 530 imind = 1, natmty(imty)
            iaty = iatype(imty,imind)
            if (masmsk(iaty,ilight)) then
               haslt(imty) = .true.
            else if (masmsk(iaty,iheavy)) then
               hashvy(imty) = .true.
            endif
  530    continue
  540 continue
c***********************************************************************
c     count the number of heavy and light atoms:
c***********************************************************************
      do 560 imol = 1, nmol
         imty = molty(imol)
         do 550 imind = 1, natmty(imty)
            iaty = iatype(imty,imind)
            if (masmsk(iaty,ilight)) then
               nlatom = nlatom + 1
            else
               nhatom = nhatom + 1
            endif
 550     continue
 560  continue
c***********************************************************************
c     count the number of atoms with fluc-q (nfqatm), mass (nmatom), and
c     Drude oscillators (ndatom); get the total charge of the system 
c     (qsys); and count up the number of charge degrees of freedom
c     (nqdof).
c***********************************************************************
      ndatom = 0
      nfqatm = 0
      nmatom = 0
      qsys = 0.0d0
      do 600 imol = 1, nmol
         imty = molty(imol)
         nfqatm = nfqatm + nfqmol(imty)
         nmatom = nmatom + nmamol(imty)
         ndatom = ndatom + ndrmol(imty)
         qsys = qsys + qmol(imty)
         if (hasfqt(imty)) then
            nqdof = nqdof + nfqmol(imty)
            if (.not. ctflag) then 
               nqdof = nqdof - 1
            endif
         endif
 600  continue
      if (ctflag) then
         nqdof = nqdof - 1
      endif
      write(iuout, *) 'There are ', nqdof, ' charge degrees of freedom.'
      if (qsys .ne. 0.0d0) then
         write(iuout, *) 'The system has a charge of ', qsys
      else
         write(iuout, *) 'The system is neutral.'
      endif
c$$$      rnkbv = dble(nmol) * k / vol
      if (.not. cbnflg .and. (qslvfl .or. dslvfl)) then
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
      endif
      if (qslvfl) then
         nqind = 0
         if (.not. ctflag) then
            do 670 imol = 1, nmol
               imty = molty(imol)
               do 660 iqind = 1, nfqmol(imty) - 1
                  iatom = iatmol(imol,ifqmol(imty,iqind))
                  nqind = nqind + 1
                  indsml(iatom) = nqind
                  indbig(nqind) = iatom
 660           continue
 670        continue
         else
            do 690 imol = 1, nmol
               imty = molty(imol)
               do 680 iqind = 1, nfqmol(imty)
                  iatom = iatmol(imol,ifqmol(imty,iqind))
                  nqind = nqind + 1
                  indsml(iatom) = nqind
                  indbig(nqind) = iatom
 680           continue
 690        continue
            nqind = nqind - 1
         endif
      endif
c***********************************************************************
c     sum up the system mass:
c***********************************************************************
      symass = 0.d0
      do 710 imol = 1, nmol
         imty = molty(imol)
         do 700 imaind = 1, nmamol(imty)
            iaty = iatype(imty,imaind)
            iatom = iatmol(imol,imamol(imty,imaind))
            rminv(iatom) = 1.0d0 / atmass(iaty)
            symass = symass + atmass(iaty)
 700     continue
 710  continue
c***********************************************************************
c     find the tag coeffs for each molecule's center of mass:
c***********************************************************************
      do 840 imty = 1, nmty
         rmolms = 0.d0
         do 800 ifind = 1, nframe(imty)
            comcof(imty,ifind) = 0.d0
 800     continue
         do 820 imaind = 1, nmamol(imty)
            imind = imamol(imty,imaind)
            iaty = iatype(imty,imind)
            if (isdrsp(iaty)) then
               go to 820
            endif
            if (isfrat(iaty)) then
               comcof(imty,imind) = comcof(imty,imind) + atmass(iaty)
            else
               do 810 ifind = 1, nframe(imty)
                  comcof(imty,ifind) = comcof(imty,ifind) + 
     $                 tagcof(imty,imind,ifind) * atmass(iaty)
 810           continue
            endif
            rmolms = rmolms + atmass(iaty)
 820     continue
         do 830 ifind = 1, nframe(imty)
            comcof(imty,ifind) = comcof(imty,ifind) / rmolms
 830     continue
 840  continue
c***********************************************************************
c     set the number of degrees of freedom and decide whether to use
c     RATTLE or not (charge degrees of freedom don't count b/c they're
c     cold; same for Drude oscillators; each constraint removes one dof;
c     keeping the system's c.o.m. motionless removes 3 dof; and if the
c     system as a whole is not rotating (not guaranteed, but generally
c     true for clusters and not for periodic systems) that's 3 more dof
c     (2 for a linear system, which I assume must be diatomic).
c***********************************************************************
      nrdof = 0
      ratflg = .false.
      do 900 imol = 1, nmol
         imty = molty(imol)
         nrdof = nrdof + 3 * nframe(imty)
         if (ncons(imty) .gt. 0) then
            nrdof = nrdof - ncons(imty)
            ratflg = .true.
         endif
         nrdof = nrdof - 3 * ndrmol(imty)
 900  continue
      nmdof = 3 * nmol
      nrdof = nrdof - 3
      nmdof = nmdof - 3
      if (.not. pbflag .and. nmol .gt. 1) then
         if (nmol .eq. 2) then
            nmdof = nmdof - 2
            if (natoms .eq. 2) then
               nrdof = nrdof - 2
            else
               nrdof = nrdof - 3
            endif
         else
            nmdof = nmdof - 3
            nrdof = nrdof - 3
         endif
c$$$         if (nmol .eq. 2 .and. natoms .eq. 2) then
c$$$            nrdof = nrdof - 2
c$$$         else
c$$$            nrdof = nrdof - 3
c$$$         endif
      endif
      if (nrdof .eq. 1) then
         write(iuout, '(a,i5,a)') 'There is ', nrdof,
     $        ' real (atomic) degree of freedom.'
      else
         write(iuout, '(a,i5,a)') 'There are ', nrdof, 
     $        ' real (atomic) degrees of freedom.'
      endif
      if (nmdof .eq. 1) then
         write(iuout, '(a,i5,a)') 'There is ', nmdof,
     $        ' molecular degree of freedom.'
      else
         write(iuout, '(a,i5,a)') 'There are ', nmdof, 
     $        ' molecular degrees of freedom.'
      endif
c***********************************************************************
c     add back in 3 dof when calculating the pressure to account for the
c     translations of the center of mass.
c***********************************************************************
      rnkbv = dble(nmdof + 3) / 3. * k / vol
      if (stfile .eq. 'none') then
         call filvel(T, iseed)
      endif
c***********************************************************************
c     set the mixed-LJ parameters, using the additive combining rule
c     for sigma:
c***********************************************************************
      do 1050 imty1 = 1, nmty
         do 1010 iljin1 = 1, nljmol(imty1)
            iaty1 = iatype(imty1,iljmol(imty1,iljin1))
            do 1000 iljin2 = 1, iljin1 - 1
               iaty2 = iatype(imty1,iljmol(imty1,iljin2))
               trmlj(iaty1,iaty2,1) = sqrt(trmlj(iaty1,iaty1,1) *
     $              trmlj(iaty2,iaty2,1))
               trmlj(iaty2,iaty1,1) = trmlj(iaty1,iaty2,1)
               trmlj(iaty1,iaty2,2) = 8. * trmlj(iaty1,iaty2,1)
               trmlj(iaty2,iaty1,2) = trmlj(iaty1,iaty2,2)
               trmlj(iaty1,iaty2,3) = 24. * trmlj(iaty1,iaty2,1)
               trmlj(iaty2,iaty1,3) = trmlj(iaty1,iaty2,3)
               trmlj(iaty1,iaty2,4) = 0.5 * (trmlj(iaty1,iaty1,4) +
     $              trmlj(iaty2,iaty2,4))
               trmlj(iaty2,iaty1,4) = trmlj(iaty1,iaty2,4)
               eps6 = trmlj(iaty1,iaty2,4) ** 3
               eps6 = eps6 * eps6
               trmlj(iaty1,iaty2,7) = eps6 * trmlj(iaty1,iaty2,2)
               trmlj(iaty2,iaty1,7) = trmlj(iaty1,iaty2,7)
               trmlj(iaty1,iaty2,8) = 3. * trmlj(iaty1,iaty2,7)
               trmlj(iaty2,iaty1,8) = trmlj(iaty1,iaty2,8)
               trmlj(iaty1,iaty2,9) = eps6 * trmlj(iaty1,iaty2,7)
               trmlj(iaty2,iaty1,9) = trmlj(iaty1,iaty2,9)
               trmlj(iaty1,iaty2,10) = 6. * trmlj(iaty1,iaty2,9)
               trmlj(iaty2,iaty1,10) = trmlj(iaty1,iaty2,10)
 1000        continue
 1010     continue
         do 1040 imty2 = 1, imty1 - 1
            do 1030 iljin1 = 1, nljmol(imty1)
               iaty1 = iatype(imty1,iljmol(imty1,iljin1))
               do 1020 iljin2 = 1, nljmol(imty2)
                  iaty2 = iatype(imty2,iljmol(imty2,iljin2))
                  trmlj(iaty1,iaty2,1) = sqrt(trmlj(iaty1,iaty1,1) *
     $                 trmlj(iaty2,iaty2,1))
                  trmlj(iaty2,iaty1,1) = trmlj(iaty1,iaty2,1)
                  trmlj(iaty1,iaty2,2) = 8. * trmlj(iaty1,iaty2,1)
                  trmlj(iaty2,iaty1,2) = trmlj(iaty1,iaty2,2)
                  trmlj(iaty1,iaty2,3) = 24. * trmlj(iaty1,iaty2,1)
                  trmlj(iaty2,iaty1,3) = trmlj(iaty1,iaty2,3)
                  trmlj(iaty1,iaty2,4) = 0.5 * (trmlj(iaty1,iaty1,4) +
     $                 trmlj(iaty2,iaty2,4))
                  trmlj(iaty2,iaty1,4) = trmlj(iaty1,iaty2,4)
                  eps6 = trmlj(iaty1,iaty2,4) ** 3
                  eps6 = eps6 * eps6
                  trmlj(iaty1,iaty2,7) = eps6 * trmlj(iaty1,iaty2,2)
                  trmlj(iaty2,iaty1,7) = trmlj(iaty1,iaty2,7)
                  trmlj(iaty1,iaty2,8) = 3. * trmlj(iaty1,iaty2,7)
                  trmlj(iaty2,iaty1,8) = trmlj(iaty1,iaty2,8)
                  trmlj(iaty1,iaty2,9) = eps6 * trmlj(iaty1,iaty2,7)
                  trmlj(iaty2,iaty1,9) = trmlj(iaty1,iaty2,9)
                  trmlj(iaty1,iaty2,10) = 6. * trmlj(iaty1,iaty2,9)
                  trmlj(iaty2,iaty1,10) = trmlj(iaty1,iaty2,10)
 1020           continue
 1030        continue
 1040     continue
 1050  continue
c***********************************************************************
c     define, count, and index the "real" atom types (those with mass, 
c     excluding Drude oscillator spring sites), 
c***********************************************************************
      nrlaty = 0
      do 1100 iaty = 1, naty
         isrl(iaty) = ismasv(iaty) .and. .not. isdrsp(iaty)
         if (isrl(iaty)) then
            nrlaty = nrlaty + 1
         endif
 1100 continue
      do 1120 imty = 1, nmty
         nrlmol(imty) = 0
         do 1110 imind = 1, natmty(imty)
            iaty = iatype(imty,imind)
            if (isrl(iaty)) then
               nrlmol(imty) = nrlmol(imty) + 1
               irlmol(imty,nrlmol(imty)) = imind
            endif
 1110    continue
 1120 continue
c***********************************************************************
c     set the info for custom (non-default) J(r) interactions:
c***********************************************************************
      if (ncint .gt. 0) then
         hscstj = .true.
         do 1210 imty1 = 1, nmty
            do 1200 imty2 = 1, nmty
               icustj(imty1,imty2) = interJ
 1200       continue
 1210    continue
         do 1220 icint = 1, ncint
            imty1 = mtymid(icjid(icint,1))
            imty2 = mtymid(icjid(icint,2))
            icustj(imty1,imty2) = icjid(icint,3)
            if (imty1 .ne. imty2) then
               icustj(imty2,imty1) = icustj(imty1,imty2)
            endif
 1220    continue
      endif
c***********************************************************************
c     fill the J() array (only if needed):
c***********************************************************************
      if ((intraJ .eq. Jreg .or. interJ .eq. Jreg .or. hscstJ) .and.
     $     neede) then
         call Jinit()
      endif
c***********************************************************************
c     make sure everything is okay:
c***********************************************************************
      call chkmod()
      return
      end
c
c***********************************************************************
c     chkmod checks a bunch of variables that are set while or after
c     reading in the model information....sjs 8/3/94
c***********************************************************************
c
      subroutine chkmod()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      logical qdlj1,
     $     hsmaid(maxaid)
c
      if (qslvfl .and. nfqaty .eq. 0) then
         write(iuout, *)
     $        'chkmod: WARNING! can''t solve for fixed charges'
         qslvfl = .false.
         qwrong = .false.
      else if (nfqaty .ne. 0 .and. .not. qslvfl) then
         write(iuout, *) 'chkmod: WARNING! doing fluc-q dynamics ',
     $        'without solving for charges'
      endif
      if (nfqaty .eq. 0 .and. interJ .ne. i1byr) then
         write(iuout, *) 'chkmod: must use 1/r intermolecular ',
     $        'interaction for fixed-q models'
         stop
      endif
      if (ndatom .eq. 0 .and. dslvfl) then
         write(iuout, *) 'chkmod: WARNING! can''t solve for Drude ',
     $        'positions if there are no Drude atoms'
         dslvfl = .false.
         dwrong = .false.
      else if (ndatom .ne. 0 .and. .not. dslvfl) then
         write(iuout, *) 'chkmod: WARNING! using Drude oscillators ',
     $        'without solving for their positions'
      endif
      do 200 imty = 1, nmty
         write(iuout, *) 'There is ', 
     $        modnam(imty)(1:index(modnam(imty),'  ')), 
     $        'in the system'
 200  continue
      do 407 imty = 1, nmty
         if (hasbdt(imty)) then
            write(iuout, *) 'chkmod: can''t do subflg stuff for ',
     $           'flexible molecules yet'
            stop
         endif
         if (hasfqt(imty) .and. hasrqt(imty)) then
            write(iuout, *) 'chkmod: mpeset can''t do mixed fluc ',
     $           'and fixed charges yet'
            stop
         endif
         do 405 imind = nframe(imty) + 1, nframe(imty) + ntag(imty)
            if (ismasv(iatype(imty,imind))) then
               write(iuout, *) 'chkmod:  not ready for massed tag ',
     $              'atoms yet (KE wrong, among others?)'
            endif
 405     continue
         do 406 iconi = 1, ncons(imty)
            iaty1 = iatype(imty,icons(imty,iconi,1))
            iaty2 = iatype(imty,icons(imty,iconi,2))
            if ((masmsk(iaty1,ilight) .and. masmsk(iaty2,iheavy)) .or.
     $           (masmsk(iaty1,iheavy) .and. masmsk(iaty2,ilight))) then
               write(iuout, *) 'chkmod:  RATTLE/RESPA combo can''t ',
     $              'handle constraints b/w light/heavy atoms'
            endif
 406     continue
         inat = nframe(imty)
         if (hasdrd(imty)) then
            inat = inat - ndrmol(imty)
         endif
         if (inat .eq. 1) then
            incon = 0
         else if (inat .eq. 2) then
            incon = 1
         else
            incon = 3 * inat - 6
         endif
         if (ncons(imty) .ne. incon) then
            write(iuout, *) 'chkmod: wrong number of constraints ',
     $           'for ', modnam(imty)(1:index(modnam(imty),'  '))
            write(iuout, *) '        (need to fix Jtype for mixed ',
     $           'flexible/fixed molecs with q-q interactions'
            write(iuout, *) '        (also need to fix comcof stuff)'
            write(iuout, *) '        (also need to fix pressure)'
            stop
         endif
 407  continue
      do 409 iatom = 1, natoms
         if (ident(iatom) .ne. ideaty(itype(iatom))) then
            write(iuout, *) 'chkmod: wrong atom type for atom ', iatom
            stop
         endif
 409  continue
      do 420 imol = 1, nmol
         imty = molty(imol)
         do 410 itind = nframe(imty) + 1, nframe(imty) + ntag(imty)
            iaty = iatype(imty,itind)
            iatom = iatmol(imol,itind)
            if (hasmaa(iatom)) then
               write(iuout, *) 
     $              'chkmod: not ready to use massive tag atoms ',
     $              '(KE, force transfer, ...)'
               stop
            endif
            if (masmsk(iaty,ilight)) then
               write(iuout, *) 'chkmod: tag atoms can''t be light'
               stop
            endif
 410     continue
         do 415 ifind = 1, nframe(imty)
            iatom = iatmol(imol,ifind)
            if (.not. hasmaa(iatom)) then
               write(iuout, *)
     $              'chkmod: not ready for massless frame atoms (KE)'
               stop
            endif
 415     continue
 420  continue
      do 425 iaid = 1, maxaid
         hsmaid(iaid) = .false.
 425  continue
      do 430 iaty = 1, naty
         if (interJ .eq. Jreg .and. isqat(iaty) .and. 
     $        .not. isztat(iaty)) then
            write(iuout, *) 'chmod: all atoms must have zeta values in',
     $           ' order to do J(r) intermolecularly'
            stop
         endif
         if (isrl(iaty)) then
            iaid = ideaty(iaty)
            if (hsmaid(iaid)) then
               write(iuout, *) 'chkmod: need to fix g(r) stuff for ',
     $              'multiple atom types with the same id, or else ',
     $              'some gxxyy files will be written over.'
               stop
            else
               hsmaid(iaid) = .true.
            endif
         endif
         qdlj1 = .false.
 430  continue
      do 530 imty = 1, nmty
         do 525 iconi = 1, ncons(imty)
            iaty1 = iatype(imty,icons(imty,iconi,1))
            iaty2 = iatype(imty,icons(imty,iconi,2))
            if (.not. ismasv(iaty1) .or. .not. ismasv(iaty2)) then
               write(iuout, *) 'chkmod:  shouldn''t have massless ',
     $              'atoms in constraints'
               stop
            endif
 525     continue
 530  continue
      return
      end
c
c***********************************************************************
c     gtcomv evaluates the system's center of mass linear momentum.
c     ...sjs 5/22/95
c***********************************************************************
c
      subroutine gtcomv(cmvel)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
      real*8 cmvel(3)
c     
      cmvel(1) = 0.d0
      cmvel(2) = 0.d0
      cmvel(3) = 0.d0
      do 170 imol = 1, nmol
         imty = molty(imol)
         do 160 imaind = 1, nmamol(imty)
            imind = imamol(imty,imaind)
            iaty = iatype(imty,imind)
            iatom = iatmol(imol,imind)
            cmvel(1) = cmvel(1) + atmass(iaty) * vel(iatom,1)
            cmvel(2) = cmvel(2) + atmass(iaty) * vel(iatom,2)
            cmvel(3) = cmvel(3) + atmass(iaty) * vel(iatom,3)
  160    continue
  170 continue
      cmvel(1) = cmvel(1) / symass
      cmvel(2) = cmvel(2) / symass
      cmvel(3) = cmvel(3) / symass
      return
      end
c
c***********************************************************************
c     setdet sets a bunch of variables that need to be (or can be) set
c     after reading in the details.dat file....sjs 8/9/94
c***********************************************************************
c
      subroutine setdet()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      qminv = 1. / qmass
      dths = dtbig
      dtl = dtbig / nl
      dtlf = dtl / nlf
      dthf = dtbig / nhf
      eqflag = .not. vscflg
      Efield(1) = Efield(1) * fE2Emd
      Efield(2) = Efield(2) * fE2Emd
      Efield(3) = Efield(3) * fE2Emd
      return
      end
c
c***********************************************************************
c     chkdet checks to see that the various variables read in from
c     details.dat are valid and compatible with each other....sjs 8/9/94
c***********************************************************************
c
      subroutine chkdet(icjid)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      integer*4 icjid(mxcstj,3)
c
      if (2 * (nl / 2) .ne. nl) then
         write(iuout, *) 'chkdet: nl must be even!'
         stop
      endif
      if (fcutr .lt. 2.0d0) then
         write(iuout, *) 'chkdet: WARNING! leading edge of cutoff ',
     $        '(fcutr) is dangerously small'
      endif
      if (fcutd .lt. 0.d0) then
         write(iuout, *) 'chkdet: cutoff range (fcutd) must be ',
     $        'non-negative'
      endif
      if (interJ .eq. Jisol) then
         write(iuout, *) 'chkdet: Jisol invalid for ',
     $        'intermolecular interactions'
         stop
      endif
      if (ncint .gt. mxcstj) then
         write(iuout, *) 'chkdet: more than ', mxcstj, 
     $        ' custom J(r) interactions found'
         stop
      endif
      do 100 icint = 1, ncint
         if (icjid(icint,3) .eq. Jisol) then
            write(iuout, *) 'chkdet: Jisol invalid for ',
     $           'intermolecular interactions'
         endif
  100 continue
      if (intraJ .eq. i1byr) then
         write(iuout, *) 'chkdet: WARNING! 1/r interactions should ',
     $        'not be used intramolecularly'
      endif
      if (ewlflg) then
         if (.not. pbflag) then
            write(iuout, *) 'chkdet: can''t do Ewald without periodic ',
     $           'boundary conditions'
            stop
         endif
      else
         if (ewlkfl) then
            write(iuout, *) 'chkdet: WARNING! can''t do Ewald k-space ',
     $           'stuff without doing Ewald stuff'
            ewlkfl = .false.
         endif
         if (ewsrfl) then
            write(iuout, *) 'chkdet: WARNING! can''t do Ewald surface ',
     $           'term without doing Ewald'
         endif
      endif
      if (stint .lt. 0.d0) then
         stint = 3000.d0
      endif
      if (vscflg .and. nscint * dtbig .lt. stint) then
         write(iuout, *) 'chkdet: WARNING! can''t rescale more often ',
     $        'than the corr. time of the temperature'
      endif
      if (vscflg .and. nsceql .lt. 15 * nscint) then
         write(iuout, *) 'chkdet: WARNING! very few T scaling ',
     $        'intervals in the T scaling equilibration time'
      endif
      if (cbnflg .and. stfile .eq. 'none') then
         write(iuout, *) 'chkdet: WARNING! cannot continue binning ',
     $        'without a startup file'
         cbnflg = .false.
      endif
      return
      end
c
cmovc***********************************************************************
cmovc     testpr does some debugging stuff on the pair energies.
cmovc     ...sjs 11/14/94
cmovc***********************************************************************
cmovc
cmov      subroutine testpr()
cmovc
cmov      include 'implic'
cmov      include 'qpar'
cmov      include 'genpar'
cmov      include 'commons'
cmovc
cmov      real*8 qbinde(maxmol), tselfe(maxmol),
cmov     $     tqbpre(maxmol,maxmol), tqepre(maxmol,maxmol),
cmov     $     tuqbpe(maxmol,maxmol)
cmovc
cmov      do 100 imol = 1, nmol
cmov         qbinde(imol) = 0.0d0
cmov 100  continue
cmov      do 290 imol1 = 1, nmol
cmov         tqbpre(imol1,imol1) = 0.5d0 * qbpre(imol1,imol1)
cmov         tqepre(imol1,imol1) = 0.5d0 * qepre(imol1,imol1) 
cmov         tuqbpe(imol1,imol1) = 0.5d0 * uqbpre(imol1,imol1)
cmov         tselfe(imol1) = 0.5d0 * selfe(imol1)
cmov         do 280 imol2 = 1, imol1 - 1
cmov            tqbpre(imol1,imol2) = qbpre(imol1,imol2) +
cmov     $           qbpre(imol2,imol1)
cmov            tqepre(imol1,imol2) = qepre(imol1,imol2) +
cmov     $           qepre(imol2,imol1)
cmov            tuqbpe(imol1,imol2) = uqbpre(imol1,imol2) +
cmov     $           uqbpre(imol2,imol1)
cmov            tqbpre(imol1,imol2) = tqbpre(imol1,imol2) * 0.5d0
cmov            tqepre(imol1,imol2) = tqepre(imol1,imol2) * 0.5d0
cmov            tuqbpe(imol1,imol2) = tuqbpe(imol1,imol2) * 0.5d0
cmov 280     continue
cmov 290  continue
cmov      do 291 imol = 1, nmol
cmov         tselfe(imol) = tselfe(imol) + tqbpre(imol,imol)
cmov         tqepre(imol,imol) = tqepre(imol,imol) - tqbpre(imol,imol)
cmov 291  continue
cmov      do 293 imol1 = 1, nmol
cmov         qbinde(imol1) = qbinde(imol1) + tqepre(imol1,imol1)
cmov         do 292 imol2 = 1, imol1 - 1
cmov            qbinde(imol1) = qbinde(imol1) + tqepre(imol1,imol2)
cmov            qbinde(imol2) = qbinde(imol2) + tqepre(imol1,imol2)
cmov 292     continue
cmov 293  continue
cmov      sdimer = 0.0d0
cmov      dimerl = bigpos
cmov      dimerh = bigneg
cmov      do 310 imol1 = 1, nmol
cmov         imty1 = molty(imol1)
cmov         do 300 imol2 = 1, imol1 - 1
cmov            imty2 = molty(imol2)
cmov            dimenr = tqbpre(imol1,imol2) + tuqbpe(imol1,imol2) +
cmov     $           tqbpre(imol1,imol2) * 
cmov     $           (tselfe(imol1) / qbinde(imol1) +
cmov     $           tselfe(imol2) / qbinde(imol2))
cmov            if (dimenr .lt. dimerl) then
cmov               dimerl = dimenr
cmov            endif
cmov            if (dimenr .gt. dimerh) then
cmov               dimerh = dimenr
cmov            endif
cmov            sdimer = sdimer + dimenr
cmov 300     continue
cmov 310  continue
cmov      sbnde = 0.0
cmov      bindel = bigpos
cmov      bindeh = bigneg
cmov      do 318 imol1 = 1, nmol
cmov         binde = tselfe(imol1)
cmov         do 313 imol2 = 1, imol1 - 1
cmov            binde = binde + tqepre(imol1,imol2) + tuqbpe(imol1,imol2)
cmov            binde = binde + tqepre(imol1,imol2) / qbinde(imol2) *
cmov     $           tselfe(imol2)
cmov 313     continue
cmov         binde = binde + tqepre(imol1,imol1) + tuqbpe(imol1,imol1)
cmov         binde = binde + tqepre(imol1,imol1) / qbinde(imol1) *
cmov     $        tselfe(imol1)
cmov         do 315 imol2 = imol1 + 1, nmol
cmov            binde = binde + tqepre(imol2,imol1) + tuqbpe(imol2,imol1)
cmov            binde = binde + tqepre(imol2,imol1) / qbinde(imol2) *
cmov     $           tselfe(imol2)
cmov 315     continue
cmov         if (binde .lt. bindel) then
cmov            bindel = binde
cmov         endif
cmov         if (binde .gt. bindeh) then
cmov            bindeh = binde
cmov         endif
cmov         sbnde = sbnde + binde
cmov 318  continue
cmov      sqbpre = 0.0d0
cmov      dsqbpr = 0.0d0
cmov      dsqepr = 0.0d0
cmov      sqepre = 0.0d0
cmov      suqbpe = 0.0d0
cmov      sselfe = 0.0d0
cmov      sqbnde = 0.0d0
cmov      selfel = bigpos
cmov      selfeh = bigneg
cmov      qbprel = bigpos
cmov      qbpreh = bigneg
cmov      uqbpel = bigpos
cmov      uqbpeh = bigneg
cmov      bpel = bigpos
cmov      bpeh = bigneg
cmov      qeprel = bigpos
cmov      qepreh = bigneg
cmov      qbndel = bigpos
cmov      qbndeh = bigneg
cmov      do 400 imol1 = 1, nmol
cmov         dsqbpr = dsqbpr + tqbpre(imol1,imol1)
cmov         suqbpe = suqbpe + tuqbpe(imol1,imol1)
cmov         sselfe = sselfe + tselfe(imol1)
cmov         sqbnde = sqbnde + qbinde(imol1)
cmov         sqepre = sqepre + tqepre(imol1,imol1)
cmov         dsqepr = dsqepr + tqepre(imol1,imol1)
cmov         if (tselfe(imol1) .lt. selfel) then
cmov            selfel = tselfe(imol1)
cmov         endif
cmov         if (tselfe(imol1) .gt. selfeh) then
cmov            selfeh = tselfe(imol1)
cmov         endif
cmov         if (qbinde(imol1) .lt. qbndel) then
cmov            qbndel = qbinde(imol1)
cmov         endif
cmov         if (qbinde(imol1) .gt. qbndeh) then
cmov            qbndeh = qbinde(imol1)
cmov         endif
cmov         if (tqepre(imol1,imol1) .lt. qeprel) then
cmov            qeprel = tqepre(imol1,imol1)
cmov         endif
cmov         if (tqepre(imol1,imol1) .gt. qepreh) then
cmov            qepreh = tqepre(imol1,imol1)
cmov         endif
cmov         if (tuqbpe(imol1,imol1) .gt. uqbpeh) then
cmov            uqbpeh = tuqbpe(imol1,imol1)
cmov         endif
cmov         if (tuqbpe(imol1,imol1) .lt. uqbpel) then
cmov            uqbpel = tuqbpe(imol1,imol1)
cmov         endif
cmov         do 390 imol2 = 1, imol1 - 1
cmov            sqbpre = sqbpre + tqbpre(imol1,imol2)
cmov            dsqepr = dsqepr + 2.0d0 * tqepre(imol1,imol2)
cmov            sqepre = sqepre + tqepre(imol1,imol2)
cmov            suqbpe = suqbpe + tuqbpe(imol1,imol2)
cmov            if (tqbpre(imol1,imol2) .lt. qbprel) then
cmov               qbprel = tqbpre(imol1,imol2)
cmov            endif
cmov            if (tqbpre(imol1,imol2) .gt. qbpreh) then
cmov               qbpreh = tqbpre(imol1,imol2)
cmov            endif
cmov            if (tuqbpe(imol1,imol2) .lt. uqbpel) then
cmov               uqbpel = tuqbpe(imol1,imol2)
cmov            endif
cmov            if (tuqbpe(imol1,imol2) .gt. uqbpeh) then
cmov               uqbpeh = tuqbpe(imol1,imol2)
cmov            endif
cmov            if (tqbpre(imol1,imol2) + tuqbpe(imol1,imol2) .lt. bpel) 
cmov     $           then
cmov               bpel = tqbpre(imol1,imol2) + tuqbpe(imol1,imol2)
cmov            endif
cmov            if (tqbpre(imol1,imol2) + tuqbpe(imol1,imol2) .gt. bpeh) 
cmov     $           then
cmov               bpeh = tqbpre(imol1,imol2) + tuqbpe(imol1,imol2)
cmov            endif
cmov            if (tqepre(imol1,imol2) .lt. qeprel) then
cmov               qeprel = tqepre(imol1,imol2)
cmov            endif
cmov            if (tqepre(imol1,imol2) .gt. qepreh) then
cmov               qepreh = tqepre(imol1,imol2)
cmov            endif
cmov 390     continue
cmov 400  continue
cmov      write(iuout, *) 'sum of binding energies =  ', 
cmov     $     sbnde / nmol * fJm2kc
cmov      write(iuout, *) 'sum of binding E check =   ',
cmov     $     (dsqepr + 2.0d0 * (suqbpe + sselfe)) 
cmov     $     / nmol * fJm2kc
cmov      write(iuout, *) 'ape =                      ', 
cmov     $     ape / nmol * fJm2kc
cmov      write(iuout, *) 'qpe =                      ',
cmov     $     qpe / nmol * fJm2kc
cmov      write(iuout, *) 'ape + qpe =                ',
cmov     $     (ape + qpe) / nmol * fJm2kc
cmov      write(iuout, *) 'system E check =           ',
cmov     $     (0.5d0 * dsqepr + suqbpe + sselfe) / nmol * fJm2kc
cmov      write(iuout, *) 'doubled system E check =   ',
cmov     $     (0.5d0 * dsqepr + suqbpe + sselfe) / nmol *
cmov     $     fJm2kc * 2
cmov      write(iuout, *) '-----------------------------------------------'
cmov      write(iuout, *) 'dsqepr                   = ',
cmov     $     dsqepr / nmol * fJm2kc
cmov      write(iuout, *) 'dsqbpr                   = ',
cmov     $     dsqbpr / nmol * fJm2kc
cmov      write(iuout, *) 'suqbpe                   = ',
cmov     $     suqbpe / nmol * fJm2kc
cmov      write(iuout, *) 'sselfe                   = ',
cmov     $     sselfe / nmol * fJm2kc
cmov      write(iuout, *) 'total charge bare pair E = ', 
cmov     $     sqbpre / nmol * fJm2kc
cmov      write(iuout, *) 'total non-q bare pair E =  ', 
cmov     $     suqbpe / nmol * fJm2kc
cmov      write(iuout, *) 'total bare pair E =        ',
cmov     $     (sqbpre + suqbpe) / nmol * fJm2kc
cmov      write(iuout, *) 'total q-only Ewald pair E =',
cmov     $     sqepre / nmol * fJm2kc
cmov      write(iuout, *) 'total q-only Ewald bind E =',
cmov     $     sqbnde / nmol * fJm2kc
cmov      write(iuout, *) 'total Ewald pair E =       ',
cmov     $     (sqepre + suqbpe) / nmol * fJm2kc
cmov      write(iuout, *) '-----------------------------------------------'
cmov      write(iuout, *) 'lowest norm. bare pair E = ',
cmov     $     dimerl * fJm2kc
cmov      write(iuout, *) 'highest norm. bare pair E =',
cmov     $     dimerh * fJm2kc
cmov      write(iuout, *) '-----------------------------------------------'
cmov      write(iuout, *) 'lowest binding energy =    ',
cmov     $     bindel * fJm2kc
cmov      write(iuout, *) 'highest binding energy =   ',
cmov     $     bindeh * fJm2kc
cmov      write(iuout, *) '-----------------------------------------------'
cmov      write(iuout, *) 'lowest q-only bare pair E =',
cmov     $     qbprel * fJm2kc
cmov      write(iuout, *) 'highest q-only bare pr E = ',
cmov     $     qbpreh * fJm2kc
cmov      write(iuout, *) '-----------------------------------------------'
cmov      write(iuout, *) 'lowest non-q bare pair E = ',
cmov     $     uqbpel * fJm2kc
cmov      write(iuout, *) 'highest non-q bare pair E =',
cmov     $     uqbpeh * fJm2kc
cmov      write(iuout, *) '-----------------------------------------------'
cmov      write(iuout, *) 'lowest total bare pair E = ',
cmov     $     bpel * fJm2kc
cmov      write(iuout, *) 'highest total bare pair E =',
cmov     $     bpeh * fJm2kc
cmov      write(iuout, *) '-----------------------------------------------'
cmov      write(iuout, *) 'lowest q-only Ewald pr E = ',
cmov     $     qeprel * fJm2kc
cmov      write(iuout, *) 'highest q-only Ewald pr E =',
cmov     $     qepreh * fJm2kc
cmov      write(iuout, *) '-----------------------------------------------'
cmov      write(iuout, *) 'lowest q-only Ewld bind E =',
cmov     $     qbndel * fJm2kc
cmov      write(iuout, *) 'highest q-only Ewd bind E =',
cmov     $     qbndeh * fJm2kc
cmov      write(iuout, *) '-----------------------------------------------'
cmov      write(iuout, *) 'lowest self energy =       ',
cmov     $     selfel * fJm2kc
cmov      write(iuout, *) 'highest self energy =      ',
cmov     $     selfeh * fJm2kc
cmov      write(iuout, *) 'average self energy =      ',
cmov     $     sselfe / nmol * fJm2kc
cmov      write(iuout, *) '-----------------------------------------------'
cmov      write(iuout, *) 'total norm. bare pair E =  ',
cmov     $     sdimer / nmol * fJm2kc
cmov      write(iuout, *) ' '
cmov      return
cmov      stop
cmov      end
cmovc
c***********************************************************************
c     dsolv solves for the position of the spring atoms on all of the
c     Drude oscillator atoms, using the electric field at the spring
c     atom of the Drude pair, and assuming the Drude atom pair to be a point
c     dipole....sjs 12/14/94
c***********************************************************************
c
      subroutine dsolv()
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
c***********************************************************************
c     initialize some stuff:
c***********************************************************************
      rms = 0.d0
c***********************************************************************
c     looping over all Drude atoms, move the spring atom to create a
c     dipole that will minimize the energy, using the field at the base
c     atom.  (r = q E / mw^2).  also keep track of a rms average:
c***********************************************************************
      do 210 imol = 1, nmol
         imty = molty(imol)
         do 200 idrind = 1, ndrmol(imty)
            ibsat = iatmol(imol,idrmol(imty,idrind,1))
            ispat = iatmol(imol,idrmol(imty,idrind,2))
            qbmw2 = q(ispat) / trmdrd(imty,idrind)
            xnew = pos(ibsat,1) + qbmw2 * field(ispat,1)
            ynew = pos(ibsat,2) + qbmw2 * field(ispat,2)
            znew = pos(ibsat,3) + qbmw2 * field(ispat,3)
            dx = xnew - pos(ispat,1)
            dy = ynew - pos(ispat,2)
            dz = znew - pos(ispat,3)
            pos(ispat,1) = xnew
            pos(ispat,2) = ynew
            pos(ispat,3) = znew
            vel(ispat,1) = vel(ibsat,1)
            vel(ispat,2) = vel(ibsat,2)
            vel(ispat,3) = vel(ibsat,3)
            rms = rms + dx * dx + dy * dy + dz * dz
 200     continue
 210  continue
      rms = sqrt(rms / ndatom)
      if (ioflg) then
         write(iuout, '(a, e10.3, a)') 
     $        'rms deviation in spring atoms = ', rms, ' A'
      endif
      if (rms .gt. drmsct) then
         dwrong = .true.
         if (qslvfl) then
            qwrong = .true.
         endif
         if (ioflg) then
            write(iuout, *) '  Drude sites not converged...'
         endif
      else
         dwrong = .false.
         if (ioflg) then
            write(iuout, *) '  Drude sites converged okay'
         endif
      endif
      movlt = .true.
      return
      end
c
c***********************************************************************
c     getdke gets the KE and temperature of the Drude oscillators.
c     ...sjs 2/7/95
c***********************************************************************
c
      subroutine getdke()
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
      do 210 imol = 1, nmol
         imty = molty(imol)
         do 200 idrind = 1, ndrmol(imty)
            ibsat = iatmol(imol,idrmol(imty,idrind,1))
            ispat = iatmol(imol,idrmol(imty,idrind,2))
            ibsaty = iatype(imty,idrmol(imty,idrind,1))
            ispaty = iatype(imty,idrmol(imty,idrind,2))
            vx = vel(ibsat,1) - vel(ispat,1)
            vy = vel(ibsat,2) - vel(ispat,2)
            vz = vel(ispat,3) - vel(ispat,3)
c$$$            redmas = rmass(ibsat) * rmass(ispat) / 
c$$$     $           (rmass(ibsat) + rmass(ispat))
            redmas = atmass(ibsaty) * atmass(ispaty) /
     $           (atmass(ibsaty) + atmass(ispaty))
            dke = dke + 0.5d0 * redmas * (vx * vx + vy * vy + vz * vz)
 200     continue
 210  continue
      drudet = 2. * dke / (3 * ndatom * k)
      return
      end
c




