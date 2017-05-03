c***********************************************************************
c     dohmbd calculates the PE and forces (actually 2PE and -dV/dx) for
c     the harmonic bonds in a given molecule....sjs 6/16/94
c***********************************************************************
c     This is copied from wattra, and I think it expects dx-type
c     variables to be dx(1,2) = x(1) - x(2), which is the opposite of
c     the way the d[xyz]mat arrays are calcd, so the signs are switched.
c     Fix this later, after it's debugged.
c***********************************************************************
c     
      subroutine dohmbd(imol, getlt, gethvy, gete)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      logical getlt, gethvy, gete
c
      logical dopr
c
cmovc***********************************************************************
cmovc     need to fix for dimer energy stuff
cmovc***********************************************************************
cmov      stop
      imty = molty(imol)
      do 200 iibd = 1, nbdmol(imty)
         imind1 = ibdmol(imty,iibd,1)
         imind2 = ibdmol(imty,iibd,2)
         iaty1 = iatype(imty,imind1)
         iaty2 = iatype(imty,imind2)
         dr = rmat(imind1,imind2) - trmbd(imty,iibd,1)
         pepc = trmbd(imty,iibd,2) * dr
         if (masmsk(iaty1,ilight) .and. masmsk(iaty2,ilight)) then
            imass = ilight
         else if (masmsk(iaty1,iheavy) .and. masmsk(iaty2,iheavy)) then
            imass = iheavy
         else
            imass = imixed
         endif
         dopr = ((imass .eq. ilight .and. getlt .and. movlt) .or.
     $        (imass .eq. imixed .and. (getlt .or. gethvy) .and. 
     $        (movlt .or. movhvy)) .or.
     $        (imass .eq. iheavy .and. gethvy .and. movhvy))
         if (.not. dopr) then
            go to 200
         endif
         if (gete) then
            ape(imass) = ape(imass) + pepc * dr
            if (dbflag) then
               ebit(3,imass) = ebit(3,imass) + pepc * dr
            endif
            if (bdeflg) then
               selfe(imol,imass) = selfe(imol,imass) + pepc * dr
            endif
         endif
cmov            if (dbflag) ebit(3) = ebit(3) + pepc * dr
         fcpc = pepc * rinmat(imind1,imind2)
         dx = dxmat(imind1,imind2)
         dy = dymat(imind1,imind2)
         dz = dzmat(imind1,imind2)
         fcx = fcpc * dx
         fcy = fcpc * dy
         fcz = fcpc * dz
         force(iatom1,1,inear,imass) = force(iatom1,1,inear,imass) + fcx
         force(iatom2,1,inear,imass) = force(iatom2,1,inear,imass) - fcx
         force(iatom1,2,inear,imass) = force(iatom1,2,inear,imass) + fcy
         force(iatom2,2,inear,imass) = force(iatom2,2,inear,imass) - fcy
         force(iatom1,3,inear,imass) = force(iatom1,3,inear,imass) + fcz
         force(iatom2,3,inear,imass) = force(iatom2,3,inear,imass) - fcz
c***********************************************************************
c     There is no pressure resulting from intramolecular forces
c***********************************************************************
 200  continue
      return
      end
c
c***********************************************************************
c     dolj calculates the PE and forces (actually 2PE and -dV/dx) for
c     the Lennard-Jones interactions between a pair of molecules.
c***********************************************************************
c     This is copied from watter, which expects dx-type variables to be
c     dx(1,2) = x(2,1), which is the way the d[xyz]mat variables are
c     calcd, so no sign change is necessary.
c***********************************************************************
c     Removed the charge-dependent LJ stuff, changed the spline stuff
c     to be a split between near and far forces, added atom mask.
c     ...sjs 3/16/95
c***********************************************************************
c     Switching function-induced forces are now applied to the atoms on 
c     which the switch is based....sjs 5/1/96
c***********************************************************************
c     The LJ potential and forces are cut off with a splined cutoff.
c***********************************************************************
c     The LJ forces are split between near and far.
c***********************************************************************
c
      subroutine dolj(imol1, imol2, getlt, gethvy, getnr, getfar,
     $     gete)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      logical getlt, gethvy, getnr, getfar, gete
c
      logical dopr
c
      imty1 = molty(imol1)
      imty2 = molty(imol2)
      do 210 iljin1 = 1, nljmol(imty1)
         imind1 = iljmol(imty1,iljin1)
         iatom1 = iatmol(imol1,imind1)
         iaty1 = iatype(imty1,imind1)
         ihead1 = iatmol(imol1,1)
         iatyh1 = iatype(imty1,1)
         do 200 iljin2 = 1, nljmol(imty2)
            imind2 = iljmol(imty2,iljin2)
            iatom2 = iatmol(imol2,imind2)
            iaty2 = iatype(imty2,imind2)
            ihead2 = iatmol(imol2,1)
            iatyh2 = iatype(imty2,1)
            if (masmsk(iaty1,ilight) .and. masmsk(iaty2,ilight) .and.
     $           masmsk(iatyh1,ilight) .and. masmsk(iatyh2,ilight)) then
               imass = ilight
            else if (masmsk(iaty1,iheavy) .and.
     $              masmsk(iaty2,iheavy) .and. 
     $              masmsk(iatyh1,iheavy) .and.
     $              masmsk(iatyh2,iheavy)) then
               imass = iheavy
            else
               imass = imixed
            endif
            dopr = ((imass .eq. ilight .and. getlt .and. movlt) .or.
     $           (imass .eq. imixed .and.
     $           (getlt .or. gethvy) .and. (movlt .or. movhvy)) .or.
     $           (imass .eq. iheavy .and. gethvy .and. movhvy))
            if (.not. dopr) then
               go to 200
            endif
            cr6 = trmlj(iaty1,iaty2,7)
            cr8 = trmlj(iaty1,iaty2,8)
            cr12 = trmlj(iaty1,iaty2,9)
            cr14 = trmlj(iaty1,iaty2,10)
            rinv = rinmat(imind1,imind2)
            rinv2 = rinv ** 2
            rinv6 = rinv2 ** 3
            rhdinv = rinmat(1,1)
            rhd2 = r2mat(1,1)
            peterm = (cr12 * rinv6 - cr6) * rinv6
            fcvpc = (cr8 - cr14 * rinv6) * rinv6 * rinv2
            fcvpc = fcvpc * spline
            fcsppc = 0.5d0 * peterm * rhdinv
            fcsppc = fcsppc * spldrv
            peterm = peterm * spline
            if (gete) then
               ape(imass) = ape(imass) + peterm
               if (dbflag) then
                  ebit(1,imass) = ebit(1,imass) + peterm
               endif
               if (bdeflg) then
                  rlje(imol1,imass) = rlje(imol1,imass) + peterm
                  rlje(imol2,imass) = rlje(imol2,imass) + peterm
               endif
            endif
cmov            if (dbflag) ebit(1) = ebit(1) + peterm
cmov            if (bniter) then
cmov               uqbpre(imol1,imol2) = uqbpre(imol1,imol2) + peterm
cmov            endif
c$$$            if (bniter) then
c$$$               dimere(imol1,imol2) = dimere(imol1,imol2) + peterm
c$$$            endif
            dx = dxmat(imind1,imind2)
            dy = dymat(imind1,imind2)
            dz = dzmat(imind1,imind2)
            dxhead = dxmat(1,1)
            dyhead = dymat(1,1)
            dzhead = dzmat(1,1)
c$$$            fcpc = fcvpc * spline + fcsppc * spldrv
            fcfx = fcvpc * dx
            fcfy = fcvpc * dy
            fcfz = fcvpc * dz
            fcspfx = fcsppc * dxhead
            fcspfy = fcsppc * dyhead
            fcspfz = fcsppc * dzhead
            if (getnr) then
c$$$               fcx = fcpc * dx * rspln
c$$$               fcy = fcpc * dy * rspln
c$$$               fcz = fcpc * dz * rspln
c$$$               fcx = fcvpc * dx * rspln
c$$$               fcy = fcvpc * dy * rspln
c$$$               fcz = fcvpc * dz * rspln
               fcx = fcfx * rspln
               fcy = fcfy * rspln
               fcz = fcfz * rspln
c$$$               fcspx = fcsppc * dxhead * rspln
c$$$               fcspy = fcsppc * dyhead * rspln
c$$$               fcspz = fcsppc * dzhead * rspln            
               fcspx = fcspfx * rspln
               fcspy = fcspfy * rspln
               fcspz = fcspfz * rspln
               force(iatom2,1,inear,imass) = 
     $              force(iatom2,1,inear,imass) - fcx
               force(iatom1,1,inear,imass) = 
     $              force(iatom1,1,inear,imass) + fcx
               force(iatom2,2,inear,imass) = 
     $              force(iatom2,2,inear,imass) - fcy
               force(iatom1,2,inear,imass) = 
     $              force(iatom1,2,inear,imass) + fcy
               force(iatom2,3,inear,imass) = 
     $              force(iatom2,3,inear,imass) - fcz
               force(iatom1,3,inear,imass) = 
     $              force(iatom1,3,inear,imass) + fcz
               force(ihead2,1,inear,imass) =
     $              force(ihead2,1,inear,imass) - fcspx
               force(ihead1,1,inear,imass) =
     $              force(ihead1,1,inear,imass) + fcspx
               force(ihead2,2,inear,imass) =
     $              force(ihead2,2,inear,imass) - fcspy
               force(ihead1,2,inear,imass) =
     $              force(ihead1,2,inear,imass) + fcspy
               force(ihead2,3,inear,imass) =
     $              force(ihead2,3,inear,imass) - fcspz
               force(ihead1,3,inear,imass) =
     $              force(ihead1,3,inear,imass) + fcspz
            endif
            if (getfar) then
c$$$               fcx = fcpc * dx * (1.d0 - rspln)
c$$$               fcy = fcpc * dy * (1.d0 - rspln)
c$$$               fcz = fcpc * dz * (1.d0 - rspln)
c$$$               fcx = fcvpc * dx * (1.d0 - rspln)
c$$$               fcy = fcvpc * dy * (1.d0 - rspln)
c$$$               fcz = fcvpc * dz * (1.d0 - rspln)
               fcx = fcfx * (1.d0 - rspln)
               fcy = fcfy * (1.d0 - rspln)
               fcz = fcfz * (1.d0 - rspln)
c$$$               fcspx = fcsppc * dxhead * (1.d0 - rspln)
c$$$               fcspy = fcsppc * dyhead * (1.d0 - rspln)
c$$$               fcspz = fcsppc * dzhead * (1.d0 - rspln)
               fcspx = fcspfx * (1.d0 - rspln)
               fcspy = fcspfy * (1.d0 - rspln)
               fcspz = fcspfz * (1.d0 - rspln)
               force(iatom2,1,ifar,imass) = 
     $              force(iatom2,1,ifar,imass) - fcx
               force(iatom1,1,ifar,imass) = 
     $              force(iatom1,1,ifar,imass) + fcx
               force(iatom2,2,ifar,imass) = 
     $              force(iatom2,2,ifar,imass) - fcy
               force(iatom1,2,ifar,imass) = 
     $              force(iatom1,2,ifar,imass) + fcy
               force(iatom2,3,ifar,imass) = 
     $              force(iatom2,3,ifar,imass) - fcz
               force(iatom1,3,ifar,imass) = 
     $              force(iatom1,3,ifar,imass) + fcz
               force(ihead2,1,ifar,imass) =
     $              force(ihead2,1,ifar,imass) - fcspx
               force(ihead1,1,ifar,imass) =
     $              force(ihead1,1,ifar,imass) + fcspx
               force(ihead2,2,ifar,imass) =
     $              force(ihead2,2,ifar,imass) - fcspy
               force(ihead1,2,ifar,imass) =
     $              force(ihead1,2,ifar,imass) + fcspy
               force(ihead2,3,ifar,imass) =
     $              force(ihead2,3,ifar,imass) - fcspz
               force(ihead1,3,ifar,imass) =
     $              force(ihead1,3,ifar,imass) + fcspz
            endif
c***********************************************************************
c     The virial as calculated here contributes positively to the 
c     pressure (i.e., no more sign changes).  Attractive forces =>
c     fcvpc > 0 => -fcfx * dxhead ~= -fcfx * dx = - fcvpc * dx**2 < 0
c     => lowered pressure.
c***********************************************************************
c     dot the force into the distance between:
c       atoms, for flexible molecules, or
c       head atoms, for rigid molecules, or
c       atom and head atom for a hybrid pair.
c     This accounts for constraint forces.  See Allen & Tildesley or 
c     paper by Smith.  There is also a contribution from the spline-
c     induced force, which always uses distance between head atoms.
c***********************************************************************
            if (ioiter .and. gete) then
               if (ncons(imty1) .gt. 0) then
                  ipind1 = 1
               else
                  ipind1 = imind1
               endif
               if (ncons(imty2) .gt. 0) then
                  ipind2 = 1
               else
                  ipind2 = imind2
               endif
               dxp = dxmat(ipind1,ipind2)
               dyp = dymat(ipind1,ipind2)
               dzp = dzmat(ipind1,ipind2)
               wvirlj(1) = wvirlj(1) - fcfx * dxp
               wvirlj(2) = wvirlj(2) - fcfy * dyp
               wvirlj(3) = wvirlj(3) - fcfz * dzp
c$$$               wvirlj = wvirlj - fcfx * dxp - fcfy * dyp - fcfz * dzp
               wvirlj(1) = wvirlj(1) - fcspfx * dxhead
               wvirlj(2) = wvirlj(2) - fcspfy * dyhead
               wvirlj(3) = wvirlj(3) - fcspfz * dzhead
c$$$               wvirlj = wvirlj - fcsppc * rhd2
            endif
 200     continue
 210  continue
      return
      end
c
c***********************************************************************
c     doqtra calculates the charge-charge energy and forces (2PE and
c     -dV/dx) inside a given molecule....sjs 6/17/94
c***********************************************************************
c
      subroutine doqtra(imol, getlt, gethvy, getnr,
     $     getfar, gete)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      logical getlt, gethvy, getnr, getfar, gete
c
      logical dopr, intfl
c
c$$$      logical subflg
c$$$c
      qijunk = 0.0d0
c$$$      subflg = .false.
      imty = molty(imol)
      do 210 iqind1 = 1, nqmol(imty)
         imind1 = iqmol(imty,iqind1)
         iatom1 = iatmol(imol,imind1)
         iaty1 = iatype(imty,imind1)
c***********************************************************************
c     if this atom is either the base or spring atom of a Drude pair,
c     set things up so that its partner will be easy to recognize:
c***********************************************************************
         if (isdrat(iaty1)) then
c$$$         if (hasdrd(imty)) then
            do 190 idrind = 1, ndrmol(imty)
               if (idrmol(imty,idrind,1) .eq. imind1) then
c$$$                  isdrat = .true.
                  idrpid = idrind
                  idroth = 2
                  go to 191
               else if (idrmol(imty,idrind,2) .eq. imind1) then
c$$$                  isdrat = .true.
                  idrpid = idrind
                  idroth = 1
                  go to 191
               endif
 190        continue
 191        continue
         endif
         do 200 iqind2 = 1, iqind1 - 1
            imind2 = iqmol(imty,iqind2)
            iatom2 = iatmol(imol,imind2)
            iaty2 = iatype(imty,imind2)
            if (masmsk(iaty1,ilight) .and. masmsk(iaty2,ilight)) then
               imass = ilight
            else if (masmsk(iaty1,iheavy) .and. 
     $              masmsk(iaty2,iheavy)) then
               imass = iheavy
            else
               imass = imixed
            endif
c***********************************************************************
c     don't bother calculating anything if we don't care about these
c     atoms right now.
c***********************************************************************
            dopr = ((imass .eq. ilight .and. getlt .and. movlt) .or.
     $           (imass .eq. imixed .and. 
     $           (getlt .or. gethvy) .and. (movlt .or. movhvy)) .or.
     $           (imass .eq. iheavy .and. gethvy .and. movhvy))
c***********************************************************************
c     The atoms don't interact if they're from a single Drude pair, or
c     neither is fluc-q and there is no flexible bond between them.  But
c     even if they don't interact, call qqjtra to take out the central
c     cell piece of the k-space Ewald sum.
c***********************************************************************
            intfl = (isfqat(iaty1) .or. isfqat(iaty2) .or. 
     $           isdrat(iaty1) .or. isdrat(iaty2) .or. hasbdt(imty))
            if (isdrat(iaty1) .and. isdrat(iaty2)) then
               intfl = intfl .and. 
     $              .not. (idrmol(imty,idrpid,idroth) .eq. imind2)
            endif
c$$$            dopr = dopr .and. ((isfqat(iaty1) .or. isfqat(iaty2)) .or.
c$$$     $           isdrat(iaty1) .or. isdrat(iaty2) .or. hasbdt(imty))
c$$$            if (isdrat(iaty1) .and. isdrat(iaty2)) then
c$$$               dopr = dopr .and. .not. 
c$$$     $              (idrmol(imty,idrpid,idroth) .eq. imind2)
c$$$            endif
c$$$            if ((.not. mask(iaty1) .and. .not. mask(iaty2)) .or.
c$$$     $           (.not. isfqat(iaty1) .and. .not. isfqat(iaty2) .and.
c$$$     $           .not. isdrat(iaty1) .and. .not. isdrat(iaty2) .and. 
c$$$     $           .not. hasbdt(imty)) .or. (isdrat(iaty1) .and. 
c$$$     $           isdrat(iaty2) .and. 
c$$$     $           idrmol(imty,idrpid,idroth) .eq. imind2)) then
c$$$               go to 199
c$$$            endif
            if (.not. dopr) then
               go to 199
            endif
c$$$            if (isdrat) then
c$$$               if (idrmol(imty,idrpid,idroth) .eq. imind2) then
c$$$                  go to 199
c$$$               endif
c$$$            endif
c$$$            qij = q(iatom1) * q(iatom2)
c$$$            dx = -dxmat(imind1,imind2)
c$$$            dy = -dymat(imind1,imind2)
c$$$            dz = -dzmat(imind1,imind2)
            call qqjtra(imol, 
c$$$     $           imind1, imind2, iatom1, iatom2, qij, qijunk,
c$$$     $           subflg)
     $           imind1, imind2, imass, getnr,
     $           getfar, intfl, gete)
 199        continue
 200     continue
 210  continue
      return
      end
c
c***********************************************************************
c     doqter calculates the charge-charge energy and forces (2PE and
c     -dV/dx) between a given pair of molecules....sjs 6/21/94
c***********************************************************************
c     
      subroutine doqter(imol1, imol2, getlt, gethvy, getnr,
     $     getfar, gete)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      logical getlt, gethvy, getnr, getfar, gete, intfl
c
      logical dopr
c
      imty1 = molty(imol1)
      imty2 = molty(imol2)
      do 210 iqind1 = 1, nqmol(imty1)
         imind1 = iqmol(imty1,iqind1)
         iaty1 = iatype(imty1,imind1)
         iatyh1 = iatype(imty1,1)
         do 200 iqind2 = 1, nqmol(imty2)
            imind2 = iqmol(imty2,iqind2)
            iaty2 = iatype(imty2,imind2)
            iatyh2 = iatype(imty2,1)
            if (masmsk(iaty1,ilight) .and. masmsk(iaty2,ilight) .and.
     $           masmsk(iatyh1,ilight) .and. masmsk(iatyh2,ilight)) then
               imass = ilight
            else if (masmsk(iaty1,iheavy) .and. 
     $              masmsk(iaty2,iheavy) .and.
     $              masmsk(iatyh1,iheavy) .and.
     $              masmsk(iatyh2,iheavy)) then
               imass = iheavy
            else
               imass = imixed
            endif
            dopr = ((imass .eq. ilight .and. getlt .and. movlt) .or.
     $           (imass .eq. imixed .and.
     $           (getlt .or. gethvy) .and. (movlt .or. movhvy)) .or.
     $           (imass .eq. iheavy .and. gethvy .and. movhvy))
            if (.not. dopr) then
               go to 200
            endif
            intfl = .true.
            call qqjter(imol1, imol2, imind1, imind2, imass, getnr,
     $           getfar, gete, intfl)
            imyun = 50 + imind1
 200     continue
 210  continue
      return
      end
c
c$$$c***********************************************************************
c$$$c     dor4 calculates the energy and forces (2PE and -dV/dx) associated
c$$$c     with the r4 and Gaussian parts of the RER potential, for a given
c$$$c     pair of molecules...sjs 6/21/94
c$$$c***********************************************************************
c$$$c     This is copied from watter, and I think it expects dx-type
c$$$c     variables to be dx(1,2) = x(2) - x(1), which is how they are
c$$$c     calculated in the d[xyz]mat arrays, so no sign change is needed.
c$$$c***********************************************************************
c$$$c
c$$$      subroutine dor4(imol1, imol2, spline, spldrv)
c$$$c
c$$$      include 'implic'
c$$$      include 'qpar'
c$$$      include 'genpar'
c$$$      include 'commons'
c$$$c
c$$$      imty = molty(imol1)
c$$$      if (imty .ne. molty(imol2)) then
c$$$         write(iuout, *) 'don''t know how to do RER potentials on ',
c$$$     $        'heterogeneous molecules'
c$$$         stop
c$$$      endif
c$$$      if (nr4mol(imty) .ne. 1) then
c$$$         write(iuout, *) 'don''t know how to do RER potentials on ',
c$$$     $        'more than one site per molecule'
c$$$         stop
c$$$      endif
c$$$      irind = ir4mol(imty,1)
c$$$      iaty = iatype(imty,irind)
c$$$      iatom1 = iatmol(imol1,irind)
c$$$      iatom2 = iatmol(imol2,irind)
c$$$      rinv2 = rinmat(irind,irind) * rinmat(irind,irind)
c$$$      rinv4 = rinv2 * rinv2
c$$$      rinv6 = rinv4 * rinv2
c$$$      rOOdev = rmat(irind,irind) - trmr4(iaty,1)
c$$$      Ogemdr = trmr4(iaty,2) * rOOdev
c$$$      fixpc = trmr4(iaty,3) * exp(-Ogemdr * rOOdev)
c$$$      fixpc2 = 2.0d0 * fixpc
c$$$      peterm = trmr4(iaty,4) * rinv4 - fixpc2
c$$$      fcpc = fixpc2 * rinmat(irind,irind) * Ogemdr - 
c$$$     $     trmr4(iaty,5) * rinv6
c$$$      if (splflg) then
c$$$         fcsppc = (0.5 * trmr4(iaty,4) * rinv4 - fixpc) *
c$$$     $        rinmat(irind,irind)
c$$$         peterm = peterm * spline
c$$$         fcpc = fcpc * spline + fcsppc * spldrv
c$$$      endif
c$$$      ape = ape + peterm
c$$$c$$$      if (bniter) then
c$$$c$$$         dimere(imol1,imol2) = dimere(imol1,imol2) + peterm
c$$$c$$$      endif
c$$$      fcx = fcpc * dxmat(irind,irind)
c$$$      fcy = fcpc * dymat(irind,irind)
c$$$      fcz = fcpc * dzmat(irind,irind)
c$$$      force(iatom2,1) = force(iatom2,1) - fcx
c$$$      force(iatom1,1) = force(iatom1,1) + fcx
c$$$      force(iatom2,2) = force(iatom2,2) - fcy
c$$$      force(iatom1,2) = force(iatom1,2) + fcy
c$$$      force(iatom2,3) = force(iatom2,3) - fcz
c$$$      force(iatom1,3) = force(iatom1,3) + fcz
c$$$cmov      if (bniter) then
c$$$cmov         uqbpre(imol1,imol2) = uqbpre(imol1,imol2) + peterm
c$$$cmov      endif
c$$$      if (ioiter) then
c$$$         wvir = wvir - fcx * dx - fcy * dy - fcz * dz
c$$$      endif
c$$$      return
c$$$      end
c$$$c
c***********************************************************************
c     mpeset calculates the PE of a monomer of the given molecule id,
c     and stores the result in the the rmonpe array.  caution:  this
c     is the correct energy, not twice the energy as is used in some
c     other energy calcs....sjs 6/16/93
c***********************************************************************
c
      subroutine mpeset(imty)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      real*8 J
c
      integer ipvt(mxatml)
      real*8  cvec(mxatml),
     $        qarr(mxatml,mxatml)
c     
      rmonpe(imty) = 0.
c***********************************************************************
c     TIP4P-MQ has a self-polarization energy.  Other fixed-charge
c     models have no self energy.
c***********************************************************************
      if (maxmid .ne. 18) then
         write(iuout, *) 'mpeset: subroutine not up to date'
         stop
      endif
      if (nfqmol(imty) .eq. 0) then
         if (midmty(imty) .eq. 7) then
            rmonpe(imty) = -5.7 * fkc2Jm
         endif
         return
      endif
      nqmlm1 = nfqmol(imty) - 1
      if (nfqmol(imty) .eq. 1) then
         imind = ifqmol(imty,1)
         iaty = iatype(imty,imind)
         qval = -trmq(iaty,1) / trmj(imty,imind,imind)
         rmonpe(imty) = rmonpe(imty) + trmq(iaty,1) * qval +
     $        0.5 * trmj(imty,imind,imind) * qval * qval
      else
         imend = ifqmol(imty,nfqmol(imty))
         iatynd = iatype(imty,ifqmol(imty,nfqmol(imty)))
         do 110 iqind1 = 1, nqmlm1
            imind1 = ifqmol(imty,iqind1)
            qarr(iqind1,iqind1) = trmj(imty,imind1,imind1) -
     $           trmj(imty,imind1,imend) - trmj(imty,imend,imind1) +
     $           trmj(imty,imend,imend)
            cvec(iqind1) = trmq(iatynd,1) - 
c$$$     $           trmq(iatype(imty,ifqmol(imty,iqind1)),1)
     $           trmq(iatype(imty,imind1),1)
            do 100 iqind2 = 1, iqind1 - 1
               imind2 = ifqmol(imty,iqind2)
               qarr(iqind1,iqind2) = trmj(imty,imind1,imind2) - 
     $              trmj(imty,imind1,imend) - trmj(imty,imend,imind2) +
     $              trmj(imty,imend,imend)
               qarr(iqind2,iqind1) = qarr(iqind1,iqind2)
 100        continue
 110     continue
         iopt = 0
c$$$         call dgef(qarr, mxatml, nqmlm1, ipvt)
c$$$         call dges(qarr, mxatml, nqmlm1, ipvt, cvec, iopt)
         call dgesv(nqmlm1, 1, qarr, mxatml, ipvt, cvec, mxatml, iexval)
         if (iexval .ne. 0) then
            write(iuout, *) 'mpeset:  dgesv failed'
            stop
         endif
         cvec(nfqmol(imty)) = 0.
         do 200 iqind = 1, nqmlm1
            cvec(nfqmol(imty)) = cvec(nfqmol(imty)) - cvec(iqind)
 200     continue
         do 300 iqind1 = 1, nfqmol(imty)
            imind1 = ifqmol(imty,iqind1)
            rmonpe(imty) = rmonpe(imty) +
     $           trmq(iatype(imty,ifqmol(imty,iqind1)),1) *
     $           cvec(iqind1)
            rmonpe(imty) = rmonpe(imty) +
     $           0.5 * trmj(imty,imind1,imind1) * 
     $           cvec(iqind1) * cvec(iqind1)
            do 290 iqind2 = 1, iqind1 - 1
               imind2 = ifqmol(imty,iqind2)
               rmonpe(imty) = rmonpe(imty) + trmj(imty,imind1,imind2) *
     $              cvec(iqind1) * cvec(iqind2)
 290        continue
 300     continue
      endif
      return
      end
c
c***********************************************************************
c     qqjter calculates the intermolecular Coulomb PE and forces from
c     one atom-atom pair.  When Ewald summing is being used, only the
c     real-space, central box PE and forces are calculated.
c     ...sjs 9/17/93
c***********************************************************************
c     added the stuff to do RESPA force split.  Left out the part
c     optimized for a 1/r intermolecular J(r), this would speed things
c     up a little bit (eliminate some J calls)....sjs 3/17/95
c***********************************************************************
c     switching function-induced forces are now applied to the atoms on
c     which the switch is based....sjs 5/1/96
c***********************************************************************
c     The RESPA split between near and far forces (both real and charge
c     forces) is accomplished with a splined switching function on the
c     forces themselves, not a switch on the potential.
c***********************************************************************
c     Overview of this routine:  take the full i-j interaction/forces 
c     and split them between near and far.  This double-counts the
c     contribution from the Ewald k-space term, so correct for this by
c     subtracting k-space forces from the far terms.  And the near
c     (real-space) contribution is smoothed to zero at the real-space
c     cutoff (different from the near/far cutoff), so only include the
c     splined contribution in the short forces.  All these fixes make
c     the "plain" 1/r contribution look not very much like 1/r.
c***********************************************************************
c
      subroutine qqjter(imol1, imol2, imind1, imind2, imass, getnr,
     $     getfar, gete, intfl)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      logical getnr, getfar, gete, intfl
c
      real*8 J
c
      
      imty1 = molty(imol1)
      imty2 = molty(imol2)
      iaty1 = iatype(imty1,imind1)
      iaty2 = iatype(imty2,imind2)
      iatom1 = iatmol(imol1,imind1)
      iatom2 = iatmol(imol2,imind2)
      ihead1 = iatmol(imol1,1)
      ihead2 = iatmol(imol2,1)
      rij = rmat(imind1,imind2)
      rij2 = r2mat(imind1,imind2)
      rijinv = rinmat(imind1,imind2)
      rhdinv = rinmat(1,1)
      rhd2 = r2mat(1,1)
c***********************************************************************
c     pkpc is (except for charges) the piece of the central cell 
c     Coulombic interaction that is accounted for in the k-space part
c     of the Ewald sum (in the limit of infinitely many k-vectors).
c     This term is only needed in calculating the energy and
c     recalculating the charges.
c***********************************************************************
      if ((ewlkfl .or. (ewlflg .and. intfl)) .and. 
     $     (gete .or. qsiter .or. getfar .or. dsiter))
     $     then
         pkpc = epsinv * rijinv * (1.d0 - berfc(ewlkap * rij))
      else
         pkpc = 0.d0
      endif
      fckspc = pkpc
c***********************************************************************
c     pepc is (except for charges) the net Coulombic energy that gets
c     added in this subroutine.  if Ewald sums are being performed, this
c     is the real-space term and is close to erfc(Kr)/r.  if necessary,
c     this piece gets smoothed down to zero at the real-space cutoff.
c***********************************************************************
      if (intfl) then
         if (hscstj) then
            intjty = icustj(imty1,imty2)
         else
            intjty = interJ
         endif
         pepc = J(imty1, imty2, imind1, imind2, rij, intjty)
      else
         pepc = 0.d0
      endif
      pepc = pepc - pkpc
      fcespc = pepc * spldrv
      pepc = pepc * spline
c***********************************************************************
c     prpc is (except for charges) the full Coulombic interaction
c     between the atom pair in question.  This would be exactly J(r)
c     except for the splined cutoff on the real-space term.
c***********************************************************************
      prpc = pepc + pkpc
      fcrspc = fcespc + fckspc
      qiqj = q(iatom1) * q(iatom2)
      qiqjx2 = 2.0d0 * qiqj
      qpetrm = qiqjx2 * pepc
c***********************************************************************
c     Put the net interaction into the PE, since it doesn't care about
c     near vs. far:
c***********************************************************************
      if (gete) then
         qpe(imass) = qpe(imass) + qpetrm
c***********************************************************************
c     The bond-energy hack includes the full 1/r pair energy, since
c     Ewald energies are excluded from the bond energy sum
c***********************************************************************
         if (bdeflg) then
            coule(imol1,imass) = coule(imol1,imass) + qiqjx2 * prpc
            coule(imol2,imass) = coule(imol2,imass) + qiqjx2 * prpc
            dimere(imol1,imol2,imass) = dimere(imol1,imol2,imass) 
     .           + qiqjx2 * prpc
            dimere(imol2,imol1,imass) = dimere(imol2,imol1,imass) 
     .           + qiqjx2 * prpc
         endif
         if (dbflag) then
            ebit(7,imass) = ebit(7,imass) + qpetrm
         endif
      endif
cmovc$$$      if (bniter) dimere(imol1,imol2) = dimere(imol1,imol2) +
cmovc$$$     $     qiqjx2 * prpc
cmovc$$$      if (bniter) dimere(imol1,imol2) = dimere(imol1,imol2) + qpetrm
cmov      if (bniter) then
cmov         bpetrm = qiqjx2 * prpc
cmov         epetrm = bpetrm - qpetrm
cmov         qbpre(imol1,imol2) = qbpre(imol1,imol2) + bpetrm
cmov         qepre(imol1,imol2) = qepre(imol1,imol2) + qpetrm
cmov      endif
c***********************************************************************
c     add dV/dq to the electronegativity:
c***********************************************************************
c***********************************************************************
c     Partition the electronegativity contribution from the full
c     (splined) interaction between the near and far parts of the
c     electronegativity based on the rRESPA switch.  A correction for 
c     the Ewald piece currently in the far part will be taken care of 
c     below, if needed.
c***********************************************************************
      if (isfqat(iaty1)) then
         if (getnr) then
            chi(iatom1,inear,imass) = chi(iatom1,inear,imass) +
     $           q(iatom2) * prpc * rspln
         endif
         if (getfar) then
            chi(iatom1,ifar,imass) = chi(iatom1,ifar,imass) +
     $           q(iatom2) * prpc * (1.d0 - rspln)
         endif
      endif
      if (isfqat(iaty2)) then
         if (getnr) then
            chi(iatom2,inear,imass) = chi(iatom2,inear,imass) +
     $           q(iatom1) * prpc * rspln
         endif
         if (getfar) then
            chi(iatom2,ifar,imass) = chi(iatom2,ifar,imass) +
     $           q(iatom1) * prpc * (1.d0 - rspln)
         endif
      endif
c***********************************************************************
c     But the Q matrix gets only the contribution from the net
c     interaction, since it doesn't care about near vs. far:
c***********************************************************************
      if (qsiter) then
         qmat(iatom1,iatom2) = qmat(iatom1,iatom2) + pepc
         qmat(iatom2,iatom1) = qmat(iatom2,iatom1) + pepc
      endif
c***********************************************************************
c     fcrpc is part of the Coulomb force between the bare pair, without
c     considering Ewald sums, accounting for cutoffs, but not the
c     cutoff-induced force.
c***********************************************************************
      if (intfl) then
         if (hscstj) then
            intjty = icustj(imty1,imty2)
         else
            intjty = interJ
         endif
         fcrpc = dJdr(imty1, imty2, imind1, imind2, rij, intjty)
      else
         fcrpc = 0.d0
      endif
      fcrpc = fcrpc * spline
c***********************************************************************
c     fckpc part of the central cell Coulomb force that is 
c     accounted for in the k-space part of the Ewald sum (in the limit 
c     of infinitely many k-vectors), corrected by the real-space cutoff,
c     but not including any cutoff-induced force.
c     We will only need this if we care about the slow (far) forces or
c     if we are recalculating the Drude positions.
c***********************************************************************
      if ((ewlkfl .or. (ewlflg .and. intfl)) .and. (getfar .or. dsiter))
     $     then
         fckpc = (r2kbrp * exp(-ewlkp2 * rij2) * epsinv - pkpc) * rijinv
      else
         fckpc = 0.d0
      endif
      fckpc = fckpc * spline
c***********************************************************************
c     fcepc is part of the net Coulombic force that is introduced in
c     this subroutine, and includes the effect of the cutoff, but not
c     any cutoff-induced force:
c***********************************************************************
      fcepc = fcrpc - fckpc
      dx = dxmat(imind1,imind2)
      dy = dymat(imind1,imind2)
      dz = dzmat(imind1,imind2)
      dxhead = dxmat(1,1)
      dyhead = dymat(1,1)
      dzhead = dzmat(1,1)
      fcrpc = fcrpc * rijinv
      fckpc = fckpc * rijinv
      fcepc = fcepc * rijinv
      fcrspc = fcrspc * rhdinv
      fckspc = fckspc * rhdinv
      fcespc = fcespc * rhdinv
c***********************************************************************
c     E = F / q.  Include the full net force when calculating the field,
c     since there is no split between near and far:
c***********************************************************************
      if (dsiter) then
         fcpc = fcepc
         fcx = fcpc * dx
         fcy = fcpc * dy
         fcz = fcpc * dz
         fcsppc = fcespc
         fcspx = fcsppc * dxhead
         fcspy = fcsppc * dyhead
         fcspz = fcsppc * dzhead
         field(iatom2,1) = field(iatom2,1) - fcx * q(iatom1)
         field(iatom1,1) = field(iatom1,1) + fcx * q(iatom2)
         field(iatom2,2) = field(iatom2,2) - fcy * q(iatom1)
         field(iatom1,2) = field(iatom1,2) + fcy * q(iatom2)
         field(iatom2,3) = field(iatom2,3) - fcz * q(iatom1)
         field(iatom1,3) = field(iatom1,3) + fcz * q(iatom2)
         field(ihead2,1) = field(ihead2,1) - fcspx * q(iatom1)
         field(ihead1,1) = field(ihead1,1) + fcspx * q(iatom2)
         field(ihead2,2) = field(ihead2,2) - fcspy * q(iatom1)
         field(ihead1,2) = field(ihead1,2) + fcspy * q(iatom2)
         field(ihead2,3) = field(ihead2,3) - fcspz * q(iatom1)
         field(ihead1,3) = field(ihead1,3) + fcspz * q(iatom2)
      endif
c***********************************************************************
c     Partition the full (splined) bare pair force between near and far.
c     A correction for the k-space far forces will be taken care of 
c     below, if needed.
c***********************************************************************
      if (getnr) then
         fcpc = qiqj * fcrpc * rspln
         fcx = fcpc * dx
         fcy = fcpc * dy
         fcz = fcpc * dz
         fcsppc = qiqj * fcrspc * rspln
         fcspx = fcsppc * dxhead
         fcspy = fcsppc * dyhead
         fcspz = fcsppc * dzhead
         force(iatom2,1,inear,imass) = force(iatom2,1,inear,imass) - fcx
         force(iatom1,1,inear,imass) = force(iatom1,1,inear,imass) + fcx
         force(iatom2,2,inear,imass) = force(iatom2,2,inear,imass) - fcy
         force(iatom1,2,inear,imass) = force(iatom1,2,inear,imass) + fcy
         force(iatom2,3,inear,imass) = force(iatom2,3,inear,imass) - fcz
         force(iatom1,3,inear,imass) = force(iatom1,3,inear,imass) + fcz
         force(ihead2,1,inear,imass) = 
     $        force(ihead2,1,inear,imass) - fcspx
         force(ihead1,1,inear,imass) =
     $        force(ihead1,1,inear,imass) + fcspx
         force(ihead2,2,inear,imass) = 
     $        force(ihead2,2,inear,imass) - fcspy
         force(ihead1,2,inear,imass) =
     $        force(ihead1,2,inear,imass) + fcspy
         force(ihead2,3,inear,imass) = 
     $        force(ihead2,3,inear,imass) - fcspz
         force(ihead1,3,inear,imass) =
     $        force(ihead1,3,inear,imass) + fcspz
      endif
      if (getfar) then
         fcpc = qiqj * fcrpc * (1.d0 - rspln)
         fcx = fcpc * dx
         fcy = fcpc * dy
         fcz = fcpc * dz
         fcsppc = qiqj * fcrspc * (1.d0 - rspln)
         fcspx = fcsppc * dxhead
         fcspy = fcsppc * dyhead
         fcspz = fcsppc * dzhead
         force(iatom2,1,ifar,imass) = force(iatom2,1,ifar,imass) - fcx
         force(iatom1,1,ifar,imass) = force(iatom1,1,ifar,imass) + fcx
         force(iatom2,2,ifar,imass) = force(iatom2,2,ifar,imass) - fcy
         force(iatom1,2,ifar,imass) = force(iatom1,2,ifar,imass) + fcy
         force(iatom2,3,ifar,imass) = force(iatom2,3,ifar,imass) - fcz
         force(iatom1,3,ifar,imass) = force(iatom1,3,ifar,imass) + fcz
         force(ihead2,1,ifar,imass) = force(ihead2,1,ifar,imass) - fcspx
         force(ihead1,1,ifar,imass) = force(ihead1,1,ifar,imass) + fcspx
         force(ihead2,2,ifar,imass) = force(ihead2,2,ifar,imass) - fcspy
         force(ihead1,2,ifar,imass) = force(ihead1,2,ifar,imass) + fcspy
         force(ihead2,3,ifar,imass) = force(ihead2,3,ifar,imass) - fcspz
         force(ihead1,3,ifar,imass) = force(ihead1,3,ifar,imass) + fcspz
      endif
c***********************************************************************
c     if necessary, subtract from the far electronegativity and far
c     force the (splined) contribution from the Ewald k-space 
c     calculation, since the full (splined) bare pair interaction was 
c     split between near and far in this subroutine:
c***********************************************************************
      if (ewlflg .and. getfar) then
         fcpc = -qiqj * fckpc 
         fcx = fcpc * dx
         fcy = fcpc * dy
         fcz = fcpc * dz
         fcsppc = -qiqj * fckspc
         fcspx = fcsppc * dxhead
         fcspy = fcsppc * dyhead
         fcspz = fcsppc * dzhead
         force(iatom2,1,ifar,imass) = force(iatom2,1,ifar,imass) - fcx
         force(iatom1,1,ifar,imass) = force(iatom1,1,ifar,imass) + fcx
         force(iatom2,2,ifar,imass) = force(iatom2,2,ifar,imass) - fcy
         force(iatom1,2,ifar,imass) = force(iatom1,2,ifar,imass) + fcy
         force(iatom2,3,ifar,imass) = force(iatom2,3,ifar,imass) - fcz
         force(iatom1,3,ifar,imass) = force(iatom1,3,ifar,imass) + fcz
         force(ihead2,1,ifar,imass) = force(ihead2,1,ifar,imass) - fcspx
         force(ihead1,1,ifar,imass) = force(ihead1,1,ifar,imass) + fcspx
         force(ihead2,2,ifar,imass) = force(ihead2,2,ifar,imass) - fcspy
         force(ihead1,2,ifar,imass) = force(ihead1,2,ifar,imass) + fcspy
         force(ihead2,3,ifar,imass) = force(ihead2,3,ifar,imass) - fcspz
         force(ihead1,3,ifar,imass) = force(ihead1,3,ifar,imass) + fcspz
         chi(iatom1,ifar,imass) = chi(iatom1,ifar,imass) -
     $        q(iatom2) * pkpc
         chi(iatom2,ifar,imass) = chi(iatom2,ifar,imass) -
     $        q(iatom1) * pkpc
      endif
c***********************************************************************
c     The virial as calculated here contributes positively to the
c     pressure (i.e., no more sign changes).  Attraction => opposite
c     charges => qiqj < 0, fcepc < 0 => fcpc > 0 => -fcx * dxhead = 
c     -fcpc * dx * dxhead ~= -fcpc * dx**2 < 0 => lowered pressure
c***********************************************************************
c     dot the force into the distance between:
c       atoms, for flexible molecules, or
c       head atoms, for rigid molecules, or
c       atom and head atom for a hybrid pair.
c     This accounts for constraint forces.  See Allen & Tildesley or 
c     paper by Smith.  There is also a contribution from the spline-
c     induced force, which always uses distance between head atoms.
c***********************************************************************
      if (ioiter .and. gete) then
         fcpc = qiqj * fcepc
         fcx = fcpc * dx
         fcy = fcpc * dy
         fcz = fcpc * dz
         if (ncons(imty1) .gt. 0) then
            ipind1 = 1
         else
            ipind1 = imind1
         endif
         if (ncons(imty2) .gt. 0) then
            ipind2 = 1
         else
            ipind2 = imind2
         endif
         dxp = dxmat(ipind1,ipind2)
         dyp = dymat(ipind1,ipind2)
         dzp = dzmat(ipind1,ipind2)
c$$$         wvirqt = wvirqt - fcx * dxp - fcy * dyp - fcz * dzp
         wvirqt(1) = wvirqt(1) - fcx * dxp
         wvirqt(2) = wvirqt(2) - fcy * dyp
         wvirqt(3) = wvirqt(3) - fcz * dzp
c$$$         wvirqt = wvirqt - qiqj * fcespc * rhd2
         wvirqt(1) = wvirqt(1) - qiqj * fcespc * dxhead ** 2
         wvirqt(2) = wvirqt(2) - qiqj * fcespc * dyhead ** 2
         wvirqt(3) = wvirqt(3) - qiqj * fcespc * dzhead ** 2
      endif
      return
      end
c
c***********************************************************************
c     qqjtra calculates the intramolecular Coulomb PE and forces for
c     one atom-atom pair....sjs 9/18/93
c***********************************************************************
c     is the sign wrong on the force and wvir calculation?
c***********************************************************************
c     dx, dy, dz aren't passed in correctly...
c***********************************************************************
c     I've axed the subflg stuff - if flexible molecules are ever used
c     again, the interaction/forces b/w bonded atoms should use dq,
c     (difference from gas-phase charges), or else use a bond potential
c     designed to work okay with the full charges.
c***********************************************************************
c
      subroutine qqjtra(imol, 
     $     imind1, imind2, imass, getnr,
     $     getfar, intfl, gete)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      logical getnr, getfar, gete, intfl
c
      real*8 J
c
      imty = molty(imol)
      iaty1 = iatype(imty,imind1)
      iaty2 = iatype(imty,imind2)
      iatom = iatmol(imol,imind1)
      iatom2 = iatmol(imol,imind2)
c$$$      id1 = ident(iatom)
c$$$      id2 = ident(iatom2)
      rij = rmat(imind1,imind2)
      rij2 = r2mat(imind1,imind2)
      rijinv = rinmat(imind1,imind2)
      qij = q(iatom) * q(iatom2)
c***********************************************************************
c     prpc is the Coulombic interaction that would be counted
c     without Ewald sums or cutoffs:
c***********************************************************************
      if (intfl) then
         prpc = J(imty, imty, imind1, imind2, rij, intraJ)
      else
         prpc = 0.d0
      endif
c***********************************************************************
c     pkpc is the Coulombic interaction between these two atoms from the
c     minimum-image interaction that is accounted for in the k-space 
c     of the Ewald sum (in the
c     limit of infinitely many k-vectors).  We will need to know it if
c     we want the energy or if we want to recalculate the charges or if
c     we will want the far forces or if we need fields to recalculate
c     the Drude displacements:
c***********************************************************************
      if ((ewlkfl .or. (ewlflg .and. intfl)) .and. 
     $     (gete .or. qsiter .or. getfar .or. dsiter)) then
         if (rij .eq. 0.d0) then
            pkpc = r2kbrp * epsinv
         else
            pkpc = (1.d0 - berfc(ewlkap * rij)) * epsinv * rijinv
         endif
      else
         pkpc = 0.d0
      endif
c***********************************************************************
c     pepc is the net Coulombic interaction between these two
c     atoms that is counted in this subroutine:
c***********************************************************************
      pepc = prpc - pkpc
      qpetrm = 2.0d0 * qij * pepc
c***********************************************************************
c     Include the net interaction in the PE:
c***********************************************************************
      if (gete) then
         qpe(imass) = qpe(imass) + qpetrm
c***********************************************************************
c     the bond energy hack includes the full J(r) pair energy, without
c     the Ewald correction, since Ewald terms are excluded from the
c     bonding energy (solvation energy) calculation
c***********************************************************************
         if (bdeflg) then
            coule(imol,imass) = coule(imol,imass) + 4.d0 * qij * prpc
            dimere(imol,imol,imass) = dimere(imol,imol,imass) 
     .           + 2.d0 * qij * prpc
         endif
c$$$         qpetra = qpetra + qpetrm
c$$$         qpetra = qpetra + qpetrm
         if (dbflag) then
            ebit(8,imass) = ebit(8,imass) + qpetrm
         endif
      endif
cmov      if (bniter) then
cmovc***********************************************************************
cmovc     In the k-space part of the Ewald sum, the full molecule i -
cmovc     molecule i interaction (twice what is included in the PE) is put
cmovc     into qepre.  Consequently, a piece of the central cell q-q 
cmovc     interaction (the part that is not included here) is also doubled,
cmovc     when it shouldn't be.  This piece is subtracted out here.  Other
cmovc     than this correction, what goes into qepre is the total molecule
cmovc     i - molecule i interaction, exactly what goes into qpe.
cmovc***********************************************************************
cmov         bpetrm = 2.0d0 * qij * pepcj
cmov         epetrm = bpetrm - qpetrm
cmov         qbpre(imol,imol) = qbpre(imol,imol) + bpetrm
cmov         qepre(imol,imol) = qepre(imol,imol) + qpetrm - epetrm
cmov      endif
c***********************************************************************
c     Put the electronegativity contribution from the bare pair
c     interaction into the near part of the electronegativity.  The
c     correction for the Ewald part will be taken care of below, if 
c     needed.
c***********************************************************************
      if (getnr) then
         if (isfqat(iaty1)) then
            chi(iatom,inear,imass) = chi(iatom,inear,imass) +
     $           q(iatom2) * prpc
         endif
         if (isfqat(iaty2)) then
            chi(iatom2,inear,imass) = chi(iatom2,inear,imass) +
     $           q(iatom) * prpc
         endif
      endif
c***********************************************************************
c     The Q matrix, on the other hand, gets its contribution from the
c     net interaction, since it isn't split between near and far:
c***********************************************************************
      if (qsiter) then
         qmat(iatom,iatom2) = qmat(iatom,iatom2) + pepc
         qmat(iatom2,iatom) = qmat(iatom2,iatom) + pepc
      endif
      dx = dxmat(imind1,imind2)
      dy = dymat(imind1,imind2)
      dz = dzmat(imind1,imind2)
c$$$      if ((dsiter .or. hasbdt(imty) .or. isdrat(iaty1) .or.
c$$$     $     isdrat(iaty2)) .and. getnr) then
c***********************************************************************
c     fcrpc is the Coulombic force from the pair interaction without
c     Ewald sums or cutoffs:
c***********************************************************************
      if (intfl) then
         fcrpc = dJdr(imty, imty, imind1, imind2, rij, intraJ)
      else
         fcrpc = 0.d0
      endif
c***********************************************************************
c     fckpc is the Coulombic force between the minimum-image atom pair
c     that is accounted for in the Ewald part of the k-space sum (in the
c     limit of infinitely many k-vectors).  We will need to know what it
c     is if we are going to get the slow (far) forces or if we are going
c     to recalculate the Drude positions:
c***********************************************************************
      if ((ewlkfl .or. (ewlflg .and. intfl)) .and. 
     $     (getfar .or. dsiter)) then
         if (rij .eq. 0.d0) then
            fckpc = 0.d0
         else
c$$$               fckpc = (pkpc - r2kbrp * exp(-ewlkp2 * rij2) * epsinv) *
            fckpc = (r2kbrp * exp(-ewlkp2 * rij2) * epsinv - pkpc) *
     $           rijinv
         endif
      else
         fckpc = 0.d0
      endif
c***********************************************************************
c     fcepc is the net Coulombic force that is counted in this
c     subroutine
c***********************************************************************
      fcepc = fcrpc - fckpc
      if (rij .eq. 0.d0) then
         if (fcepc .ne. 0.d0) then
            write(iuout, *) 'qqjtra: fix for overlapping Drudes'
            stop
         endif
         if (fcrpc .ne. 0.d0) then
            write(iuout, *) 'qqjtra: fix for overlapping Drudes'
            stop
         endif
         if (fckpc .ne. 0.d0) then
            write(iuout, *) 'qqjtra: fix for overlapping Drudes'
            stop
         endif
      else
         fcepc = fcepc * rijinv
         fcrpc = fcrpc * rijinv
         fckpc = fckpc * rijinv
      endif
c$$$      endif
c***********************************************************************
c     E = F / q.  Use the net force, since the field array isn't broken
c     up into near and far:
c***********************************************************************
      if (dsiter) then
         fcpc = fcepc
         fcx = fcpc * dx
         fcy = fcpc * dy
         fcz = fcpc * dz
         field(iatom2,1) = field(iatom2,1) - fcx * q(iatom)
         field(iatom,1) = field(iatom,1) + fcx * q(iatom2)
         field(iatom2,2) = field(iatom2,2) - fcy * q(iatom)
         field(iatom,2) = field(iatom,2) + fcy * q(iatom2)
         field(iatom2,3) = field(iatom2,3) - fcz * q(iatom)
         field(iatom,3) = field(iatom,3) + fcz * q(iatom2)
      endif
c***********************************************************************
c     Only flexible or Drude atoms need forces applied intramolecularly.
c     Apply the full bare pair force to the near part of the force 
c     array.  The correction for the Ewald part in the far force array
c     will be done below, if needed.
c***********************************************************************
      if ((hasbdt(imty) .or. isdrat(iaty1) .or. isdrat(iaty2)) .and.
     $     getnr) then
         fcpc = qij * fcrpc
         fcx = fcpc * dx
         fcy = fcpc * dy
         fcz = fcpc * dz
         force(iatom2,1,inear,imass) = force(iatom2,1,inear,imass) - fcx
         force(iatom,1,inear,imass) = force(iatom,1,inear,imass) + fcx
         force(iatom2,2,inear,imass) = force(iatom2,2,inear,imass) - fcy
         force(iatom,2,inear,imass) = force(iatom,2,inear,imass) + fcy
         force(iatom2,3,inear,imass) = force(iatom2,3,inear,imass) - fcz
         force(iatom,3,inear,imass) = force(iatom,3,inear,imass) + fcz
c***********************************************************************
c     There is no pressure resulting from this component of the force,
c     since it is purely intramolecular.
c***********************************************************************
      endif
c***********************************************************************
c     if necessary, subtract out the part of the Ewald k-space 
c     interaction that represents the central cell interaction, since
c     it was mistakenly attributed to the far force array.  Note that
c     whether this correction gets applied depends on whether there was
c     an Ewald contribution b/w these two atoms, not on whether there is
c     an intramolecular force between them.
c***********************************************************************
      if (ewlflg .and. getfar) then
         fcpc = -qij * fckpc
         fcx = fcpc * dx
         fcy = fcpc * dy
         fcz = fcpc * dz
         force(iatom2,1,ifar,imass) = force(iatom2,1,ifar,imass) - 
     $        fcx
         force(iatom,1,ifar,imass) = force(iatom,1,ifar,imass) + fcx
         force(iatom2,2,ifar,imass) = force(iatom2,2,ifar,imass) - 
     $        fcy
         force(iatom,2,ifar,imass) = force(iatom,2,ifar,imass) + fcy
         force(iatom2,3,ifar,imass) = force(iatom2,3,ifar,imass) - 
     $        fcz
         force(iatom,3,ifar,imass) = force(iatom,3,ifar,imass) + fcz
         chi(iatom,ifar,imass) = chi(iatom,ifar,imass) -
     $        q(iatom2) * pkpc
         chi(iatom2,ifar,imass) = chi(iatom2,ifar,imass) -
     $        q(iatom) * pkpc
c***********************************************************************
c     This is a correction for a force that was erroneously included in
c     the Ewald k-space calculation.  Thus a pressure component was
c     erroneously included as well, and must be removed here.
c***********************************************************************
c     The virial as calculated here contributes positively to the
c     pressure (i.e., no more sign changes).  Attraction => opposite
c     charges => qij < 0, fckpc < 0 => fcpc < 0 => -fcx * dx = 
c     -fcpc * dx ** 2 > 0 => increased pressure to compensate for 
c     lowering of pressure in Ewald sum.
c***********************************************************************
c     dot the force into the distance between:
c       atoms, for a flexible molecule, or
c       zero, for a rigid molecule.
c     There is no contribution from a spline-induced force, since the
c     force is intramolecular (and the distance between head atoms is
c     zero anyway).
c***********************************************************************
         if (ioiter .and. gete .and. ncons(imty) .eq. 0) then
            fcpc = -qij * fcrpc
            fcx = fcpc * dx
            fcy = fcpc * dy
            fcz = fcpc * dz
c$$$            wvirqa = wvirqa - fcx * dx - fcy * dy - fcz * dz
            wvirqa(1) = wvirqa(1) - fcx * dx
            wvirqa(2) = wvirqa(2) - fcy * dy
            wvirqa(3) = wvirqa(3) - fcz * dz
         endif
      endif
      return
      end
c
c***********************************************************************
c     dodrud calculates the PE (actually 2PE) and forces for the
c     harmonic spring term in a molecule with Drude oscillator atoms.
c     ...sjs 10/21/94
c***********************************************************************
      subroutine dodrud(imol, getlt, gethvy, gete)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c     
      logical getlt, gethvy, gete
c
      logical dopr
c
      imty = molty(imol)
      do 200 idrind = 1, ndrmol(imty)
         ibsind = idrmol(imty,idrind,1)
         iatyb = iatype(imty,ibsind)
         ispind = idrmol(imty,idrind,2)
         iatys = iatype(imty,ispind)
c$$$         if (mask(iatyb) .or. mask(iatys)) then
         if (masmsk(iatyb,ilight) .and. masmsk(iatys,ilight)) then
            imass = ilight
         else if (masmsk(iatyb,iheavy) .and. masmsk(iatys,iheavy)) then
            imass = iheavy
         else
            imass = imixed
         endif
         dopr = ((imass .eq. ilight .and. getlt .and. movlt) .or.
     $        (imass .eq. imixed .and. 
     $        (getlt .or. gethvy) .and. (movlt .or. movhvy)) .or.
     $        (imass .eq. iheavy .and. gethvy .and. movhvy))
         if (.not. dopr) then
            go to 200
         endif
         ibsat = iatmol(imol,ibsind)
         ispat = iatmol(imol,ispind)
         peterm = trmdrd(imty,idrind) * r2mat(ispind,ibsind)
         if (gete) then
            ape(imass) = ape(imass) + peterm
            if (dbflag) then
               ebit(2,imass) = ebit(2,imass) + peterm
            endif
            if (bdeflg) then
               selfe(imol,imass) = selfe(imol,imass) + peterm
            endif
         endif
cmov            if (dbflag) ebit(2) = ebit(2) + peterm
         dx = dxmat(ispind,ibsind)
         dy = dymat(ispind,ibsind)
         dz = dzmat(ispind,ibsind)
         fcx = trmdrd(imty,idrind) * dx
         fcy = trmdrd(imty,idrind) * dy
         fcz = trmdrd(imty,idrind) * dz
         force(ispat,1,inear,imass) = force(ispat,1,inear,imass) + fcx
         force(ibsat,1,inear,imass) = force(ibsat,1,inear,imass) - fcx
         force(ispat,2,inear,imass) = force(ispat,2,inear,imass) + fcy
         force(ibsat,2,inear,imass) = force(ibsat,2,inear,imass) - fcy
         force(ispat,3,inear,imass) = force(ispat,3,inear,imass) + fcz
         force(ibsat,3,inear,imass) = force(ibsat,3,inear,imass) - fcz
c***********************************************************************
c     There is no pressure resulting from intramolecular forces.
c***********************************************************************
 200  continue
      return
      end
c
c***********************************************************************
c     getf calls the appropriate energy/force routines, depending on
c     the atom-type mask passed to it and whether long- and/or short-
c     range forces are needed.  (Currently the atom-type mask is
c     ignored.)...sjs 3/10/95
c***********************************************************************
c
      subroutine getf(getnr, getfar, getlt, gethvy, gete)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
      logical getnr, getfar, getlt, gethvy, incut, gete, 
     $     doone, dotwo, dopr, domass, dodist
c     
c***********************************************************************
c     Zero the various forces and energies that will be accumulated 
c     here:
c***********************************************************************
      do 120 imass = 1, nmass
c***********************************************************************
c     We won't bother to calculate or sum a component of the force
c     array if (a) we don't care about it for this step, or (b) none
c     of the atoms have moved since the last time we summed it:
c***********************************************************************
         domass = ((imass .eq. ilight .and. getlt .and. movlt) .or.
     $        (imass .eq. imixed .and.
     $        (getlt .or. gethvy) .and. (movlt .or. movhvy)) .or.
     $        (imass .eq. iheavy .and. gethvy .and. movhvy))
         if (.not. domass) then
c$$$         if ((imass .eq. ilight .and. 
c$$$     $        (.not. getlt .or. .not. movlt)) .or.
c$$$     $        (imass .eq. imixed .and. 
c$$$     $        ((.not. getlt .and. .not. gethvy) .or. 
c$$$     $        (.not. movlt .and. .not. movhvy))) .or.
c$$$     $        (imass .eq. iheavy .and. 
c$$$     $        (.not. gethvy .or. .not. movhvy))) then
            go to 120
         endif
         do 110 idist = 1, ndist
            dodist = ((idist .eq. inear .and. getnr) .or.
     $           (idist .eq. ifar .and. getfar))
            if (.not. dodist) then
c$$$            if ((idist .eq. inear .and. .not. getnr) .or.
c$$$     $           (idist .eq. ifar .and. .not. getfar)) then
               go to 110
            endif
            do 100 iatom = 1, natoms
               force(iatom,1,idist,imass) = 0.d0
               force(iatom,2,idist,imass) = 0.d0
               force(iatom,3,idist,imass) = 0.d0
               qforce(iatom,idist,imass) = 0.d0
               chi(iatom,idist,imass) = 0.d0
 100        continue
 110     continue
 120  continue
      if (dsiter) then
         do 130 iatom = 1, natoms
            field(iatom,1) = 0.d0
            field(iatom,2) = 0.d0
            field(iatom,3) = 0.d0
 130     continue
      endif
      if (gete) then
c$$$         qpe = -basepe
         if (bdeflg) then
            do 9090 imol = 1, nmol
               solve(imol) = 0.d0
 9090       continue
         endif
         if (getlt .and. movlt) then
            qpe(ilight) = 0.d0
            ape(ilight) = 0.d0
            if (bdeflg) then
               do 131 imol = 1, nmol
                  selfe(imol,ilight) = 0.d0
                  coule(imol,ilight) = 0.d0
                  rlje(imol,ilight) = 0.d0
                  do 1311 jmol = 1, nmol
                     dimere(imol,jmol,ilight) = 0.d0
 1311             continue
 131           continue
            endif
         endif
         if ((getlt .or. gethvy) .and. (movlt .or. movhvy)) then
            qpe(imixed) = 0.d0
            ape(imixed) = 0.d0
            if (bdeflg) then
               do 132 imol = 1, nmol
                  selfe(imol,imixed) = 0.d0
                  coule(imol,imixed) = 0.d0
                  rlje(imol,imixed) = 0.d0
                  do 1321 jmol = 1, nmol
                     dimere(imol,jmol,imixed) = 0.d0
 1321             continue
 132           continue
            endif
         endif
         if (gethvy .and. movhvy) then
            qpe(iheavy) = -basepe
            ape(iheavy) = 0.d0
            if (bdeflg) then
               do 133 imol = 1, nmol
                  selfe(imol,iheavy) = -2.d0 * rmonpe(molty(imol))
                  coule(imol,iheavy) = 0.d0
                  rlje(imol,iheavy) = 0.d0
                  do 1331 jmol = 1, nmol
                     dimere(imol,jmol,iheavy) = 0.d0
 1331             continue
 133           continue
            endif
         endif
c$$$         qpeter = 0.d0
c$$$         qpetra = 0.d0
c$$$         ape = 0.d0
         if (dbflag) then
            do 135 ii = 1, 10
               if (getlt .and. movlt) then
                  ebit(ii,ilight) = 0.d0
               endif
               if ((getlt .or. gethvy) .and. (movlt .or. movhvy)) then
                  ebit(ii,imixed) = 0.d0
               endif
               if (gethvy .and. movhvy) then
                  ebit(ii,iheavy) = 0.d0
                  ebit(9,iheavy) = -basepe
               endif
  135       continue
         endif
      endif
      if (ioiter .and. gete) then
c$$$         wvirlj = 0.d0
         wvirlj(1) = 0.d0
         wvirlj(2) = 0.d0
         wvirlj(3) = 0.d0
c$$$         wvirqt = 0.d0
         wvirqt(1) = 0.d0
         wvirqt(2) = 0.d0
         wvirqt(3) = 0.d0
c$$$         wvirqa = 0.d0
         wvirqa(1) = 0.d0
         wvirqa(2) = 0.d0
         wvirqa(3) = 0.d0
         wvirk = 0.d0
         wvirkd = 0.d0
         wvirdr = 0.d0
      endif
      if (qsiter) then
         do 170 imol1 = 1, nmol
            imty1 = molty(imol1)
            do 160 iqind1 = 1, nqmol(imty1)
               iatom1 = iatmol(imol1,iqmol(imty1,iqind1))
               qvec(iatom1) = 0.0d0
               do 150 imol2 = 1, nmol
                  imty2 = molty(imol2)
                  do 140 iqind2 = 1, nqmol(imty2)
                     iatom2 = iatmol(imol2,iqmol(imty2,iqind2))
                     qmat(iatom2,iatom1) = 0.0d0
 140              continue
 150           continue
 160        continue
 170     continue
      endif
c***********************************************************************
c     These are the Ewald forces.  They are all long-ranged:
c***********************************************************************
      if (ewlflg .and. getfar) then
         call ewlrcp(getlt, gethvy, gete)
         call ewlslf(getlt, gethvy, gete)
         if (ewsrfl) then
            call ewlsrf(getlt, gethvy, gete)
         endif
      endif
c***********************************************************************
c     loop over molecules, calculating intramolecular interactions:
c***********************************************************************
      do 510 imol1 = 1, nmol
         imty1 = molty(imol1)
c***********************************************************************
c     don't bother calculating any intramolecular interactions if the
c     molecule has no atoms of the mass type we care about at the
c     moment, or if they haven't moved:
c***********************************************************************
         doone = ((haslt(imty1) .and. getlt .and. movlt) .or. 
     $            ((haslt(imty1) .and. hashvy(imty1)) .and. 
     $             (getlt .or. gethvy) .and. (movlt .or. movhvy)) .or.
     $            (hashvy(imty1) .and. gethvy .and. movhvy))
c***********************************************************************
c     These are the intramolecular interactions, which are
c     all short-ranged:
c***********************************************************************
         if (getnr .and. doone) then
            call doqslf(imol1, getlt, gethvy, gete)
            call brmono(imol1)
            if (hasbdt(imty1)) then
               call dohmbd(imol1, getlt, gethvy, gete)
            endif
            if (hasdrd(imty1)) then
               call dodrud(imol1, getlt, gethvy, gete)
            endif
         endif
c***********************************************************************
c     intramolecular Coulomb interactions are always short-ranged, but 
c     when doing Ewald, a correction has to be made for what has already
c     been included in the far part by the k-space Ewald sum:
c***********************************************************************
         if (doone) then
            call doqtra(imol1, getlt, gethvy, getnr, getfar, gete)
         endif
c***********************************************************************
c     Now do a second loop over molecules, calculating intermolecular
c     interactions:
c***********************************************************************
         do 500 imol2 = 1, imol1 - 1
            imty2 = molty(imol2)
c***********************************************************************
c     don't bother if the interactions between these molecules either
c     are ones we don't care about, or if the molecules haven't moved:
c***********************************************************************
            dopr = (((haslt(imty1) .and. haslt(imty2)) .and. 
     $               getlt .and. movlt) .or.
     $              (((haslt(imty1) .and. hashvy(imty2)) .or. 
     $                (hashvy(imty1) .and. haslt(imty2))) .and. 
     $               (getlt .or. gethvy) .and. (movlt .or. movhvy)) .or.
     $              ((hashvy(imty1) .and. hashvy(imty2)) .and.
     $               gethvy .and. movhvy))
c$$$            dotwo = ((haslt(imty2) .and. getlt) .or.
c$$$     $           (hashvy(imty2) .and. gethvy))
c$$$            dopr = (doone .or. dotwo)
            if (.not. dopr) then
               go to 500
            endif
c***********************************************************************
c     occasionally, recheck to see whether the pair belongs on the 
c     pair list or not.  fcutr is the near edge of the cutoff, fcutd
c     is the width of the splined cutoff region, and fcutb is the width
c     of the buffer region.
c***********************************************************************
            if (pliter) then
               call brhead(imol1, imol2, dxinc, dyinc, dzinc)
               if (rmat(1,1) .lt. fcutr + fcutd + fcutb) then
                  lclose(imol1,imol2) = .true.
                  call brpair(imol1, imol2, dxinc, dyinc, dzinc)
               else
                  lclose(imol1,imol2) = .false.
                  go to 500
               endif
c***********************************************************************
c     If it's not time for a pair list update, only calculate distances
c     if the pair is on the pairlist:
c***********************************************************************
            else
               if (lclose(imol1,imol2)) then
                  call brhead(imol1, imol2, dxinc, dyinc, dzinc)
                  call brpair(imol1, imol2, dxinc, dyinc, dzinc)
               else
                  go to 500
               endif
            endif
c***********************************************************************
c     keep two sets of splines.  spline and spldrv determine how much of
c     the interaction gets calculated at all; rspln determines
c     how the force gets partitioned between near and far in the rRESPA
c     split.  Both get applied at a molecular level.
c***********************************************************************
            if (rmat(1,1) .lt. fcutr) then
c$$$               tgetnr = .true. .and. getnr
c$$$               tgetfr = .false. .and. getfar
               incut = .true.
               spline = 1.d0
               spldrv = 0.d0
            else if (rmat(1,1) .gt. fcutr + fcutd) then
c$$$               tgetnr = .false. .and. getnr
c$$$               tgetfr = .true. .and. getfar
c$$$               go to 500
               incut = .false.
               spline = 0.d0
               spldrv = 0.d0
            else
               incut = .true.
               splr = (rmat(1,1) - fcutr) / fcutd
               spline = 1.d0 + splr * splr * (2.d0 * splr - 3.d0)
               spldrv = 6.d0 / fcutd * splr * (splr - 1.d0)
            endif
            if (rmat(1,1) .lt. fnearr) then
               rspln = 1.d0
            else if (rmat(1,1) .gt. fnearr + fneard) then
               rspln = 0.d0
            else
               rsplr = (rmat(1,1) - fnearr) / fneard
               rspln = 1.d0 + rsplr * rsplr * (2.d0 * rsplr - 3.d0)
            endif
c***********************************************************************
c     Lennard-Jones and Coulomb interactions are switched off outside
c     the cutoff so are only needed if the pair is inside fcutr+fcutd:
c***********************************************************************
            if (incut) then
               call dolj(imol1, imol2, getlt, gethvy, getnr, getfar,
     $              gete)
               call doqter(imol1, imol2, getlt, gethvy, getnr,
     $              getfar, gete)
            endif
 500     continue
 510  continue
c***********************************************************************
c     find the charge forces that will satisfy the proper charge
c     constraints:
c***********************************************************************
      call doqcon(getnr, getfar)
c***********************************************************************
c     distribute the forces onto frame atoms:
c***********************************************************************
      call fdist(getlt, gethvy, getnr, getfar)
c***********************************************************************
c     we actually calculated 2*PE, so halve it:
c***********************************************************************
      if (gete) then
         if (getlt .and. movlt) then
            qpe(ilight) = 0.5d0 * qpe(ilight)
            ape(ilight) = 0.5d0 * ape(ilight) + qpe(ilight)
            if (bdeflg) then
               do 600 imol = 1, nmol
                  coule(imol,ilight) = 0.5d0 * coule(imol,ilight)
                  rlje(imol,ilight) = 0.5d0 * rlje(imol,ilight)
                  selfe(imol,ilight) = 0.5d0 * selfe(imol,ilight)
                  do 599 jmol = 1, nmol 
                     dimere(imol,jmol,ilight) = 0.5d0 
     .                    * dimere(imol,jmol,ilight)
 599              continue
 600           continue
            endif
         endif
         if ((getlt .or. gethvy) .and. (movlt .or. movhvy)) then
            qpe(imixed) = 0.5d0 * qpe(imixed)
            ape(imixed) = 0.5d0 * ape(imixed) + qpe(imixed)
            if (bdeflg) then
               do 610 imol = 1, nmol
                  coule(imol,imixed) = 0.5d0 * coule(imol,imixed)
                  rlje(imol,imixed) = 0.5d0 * rlje(imol,imixed)
                  selfe(imol,imixed) = 0.5d0 * selfe(imol,imixed)
                  do 609 jmol = 1, nmol
                     dimere(imol,jmol,imixed) = 0.5d0
     .                    * dimere(imol,jmol,imixed)
 609              continue
 610           continue
            endif
         endif
         if (gethvy .and. movhvy) then
            qpe(iheavy) = 0.5d0 * qpe(iheavy)
            ape(iheavy) = 0.5d0 * ape(iheavy) + qpe(iheavy)
            if (bdeflg) then
               do 620 imol = 1, nmol
                  coule(imol,iheavy) = 0.5d0 * coule(imol,iheavy)
                  rlje(imol,iheavy) = 0.5d0 * rlje(imol,iheavy)
                  selfe(imol,iheavy) = 0.5d0 * selfe(imol,iheavy)
                  do 619 jmol = 1, nmol
                     dimere(imol,jmol,iheavy) = 0.5d0
     .                    * dimere(imol,jmol,iheavy)
 619              continue
 620           continue
            endif
         endif
         if (dbflag) then
            do 1510 ii = 1, 10
               if (getlt .and. movlt) then
                  ebit(ii,ilight) = 0.5d0 * ebit(ii,ilight)
               endif
               if ((getlt .or. gethvy) .and. (movlt .or. movhvy)) then
                  ebit(ii,imixed) = 0.5d0 * ebit(ii,imixed)
               endif
               if (gethvy .and. movhvy) then
                  ebit(ii,iheavy) = 0.5d0 * ebit(ii,iheavy)
               endif
 1510       continue
         endif
      endif
c***********************************************************************
c     reset the flags for any forces we've just calculated:
c***********************************************************************
      if (getlt .and. movlt) then
         movlt = .false.
      endif
      if (gethvy .and. movhvy) then
         movhvy = .false.
      endif
      return
      end
c
c***********************************************************************
c     doqslf calculates the electronegativity equalization self terms.
c     ...sjs 3/17/95
c***********************************************************************
c
      subroutine doqslf(imol, getlt, gethvy, gete)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
      logical getlt, gethvy, gete
c
      imty = molty(imol)
      do 200 iqind = 1, nfqmol(imty)
         imind = ifqmol(imty,iqind)
         iaty = iatype(imty,imind)
c$$$         if (mask(iaty)) then
         iatom = iatmol(imol,imind)
         pepc = trmq(iaty,1) + q(iatom) * trmj(imty,imind,imind)
c$$$            chinr(iatom) = chinr(iatom) + pepc
         if (masmsk(iaty,ilight) .and. getlt .and. movlt) then
            chi(iatom,inear,ilight) = chi(iatom,inear,ilight) + pepc
         else if (masmsk(iaty,iheavy) .and. gethvy .and. movhvy) then
            chi(iatom,inear,iheavy) = chi(iatom,inear,iheavy) + pepc
         endif
c$$$            qpetrm = 2.d0 * q(iatom) * (trmq(iaty,1) +
c$$$     $           q(iatom) * trmj(imty,imind,imind))
         if (gete) then
            if (masmsk(iaty,ilight) .and. getlt .and. movlt) then
               qpe(ilight) = qpe(ilight) + 
     $              q(iatom) * (trmq(iaty,1) + pepc)
               if (bdeflg) then
                  selfe(imol,ilight) = selfe(imol,ilight) + 
     $                 q(iatom) * (trmq(iaty,1) + pepc)
               endif
            else if (masmsk(iaty,iheavy) .and. gethvy .and. movhvy) then
               qpe(iheavy) = qpe(iheavy) +
     $              q(iatom) * (trmq(iaty,1) + pepc)
c$$$               qpe = qpe + q(iatom) * (trmq(iaty,1) + pepc)
               if (bdeflg) then
                  selfe(imol,iheavy) = selfe(imol,iheavy) +
     $                 q(iatom) * (trmq(iaty,1) + pepc)
               endif
            endif
c$$$            qpetra = qpetra + q(iatom) * (trmq(iaty,1) + pepc)
            if (dbflag) then
               if (masmsk(iaty,ilight) .and. getlt .and. movlt) then
                  ebit(10,ilight) = ebit(10,ilight) +
     $                 q(iatom) * (trmq(iaty,1) + pepc)
               else if (masmsk(iaty,iheavy) .and. 
     $                 gethvy .and. movhvy) then
                  ebit(10,iheavy) = ebit(10,iheavy) + 
     $                 q(iatom) * (trmq(iaty,1) + pepc)
               endif
            endif
         endif
cmov            if (bniter) then
cmov               selfe(imol) = selfe(imol) + qpetrm
cmov            endif
         if (qsiter) then
            qmat(iatom,iatom) = qmat(iatom,iatom) +
     $           trmj(imty,imind,imind)
         endif
c$$$         endif
 200  continue
      return
      end
c
c***********************************************************************
c     doqcon readjusts the charge forces (electronegativities) to ensure
c     that charge is conserved for either the system or for each 
c     molecule....sjs 3/17/95
c***********************************************************************
c
      subroutine doqcon(getnr, getfar)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
      logical getnr, getfar
c
      real*8 chitot(ndist,nmass), chiave(ndist,nmass)
c
c***********************************************************************
c     the charge force on a single charge that will keep the charge 
c     constraint(s) satisfied is the average electronegativity (chi)
c     over the set of constrained charges minus the electronegativity
c     on the charge of interest.  To keep the total (near, far; light,
c     mixed, heavy) charge force equal to what it should be, the average
c     chi should be averaged only over the appropriate contributions.
c***********************************************************************
c     But only update the stuff that has changed.
c***********************************************************************
      if (ctflag) then
         if (movlt) then
            chitot(inear,ilight) = 0.d0
            chitot(ifar,ilight) = 0.d0
         endif
         if (movlt .or. movhvy) then
            chitot(inear,imixed) = 0.d0
            chitot(ifar,imixed) = 0.d0
         endif
         if (movhvy) then
            chitot(inear,iheavy) = 0.d0
            chitot(ifar,iheavy) = 0.d0
         endif
         do 210 imol = 1, nmol
            imty = molty(imol)
            do 200 iqind = 1, nfqmol(imty)
               iatom = iatmol(imol,ifqmol(imty,iqind))
               if (getnr) then
                  if (movlt) then
                     chitot(inear,ilight) = chitot(inear,ilight) +
     $                    chi(iatom,inear,ilight)
                  endif
                  if (movlt .or. movhvy) then
                     chitot(inear,imixed) = chitot(inear,imixed) +
     $                    chi(iatom,inear,imixed)
                  endif
                  if (movhvy) then
                     chitot(inear,iheavy) = chitot(inear,iheavy) +
     $                    chi(iatom,inear,iheavy)
                  endif
               endif
               if (getfar) then
                  if (movlt) then
                     chitot(ifar,ilight) = chitot(ifar,ilight) +
     $                    chi(iatom,ifar,ilight)
                  endif
                  if (movlt .or. movhvy) then
                     chitot(ifar,imixed) = chitot(ifar,imixed) +
     $                    chi(iatom,ifar,imixed)
                  endif
                  if (movhvy) then
                     chitot(ifar,iheavy) = chitot(ifar,iheavy) +
     $                    chi(iatom,ifar,iheavy)
                  endif
               endif
 200        continue
 210     continue
         if (movlt) then
            chiave(inear,ilight) = chitot(inear,ilight) / nfqatm
            chiave(ifar,ilight) = chitot(ifar,ilight) / nfqatm
         endif
         if (movlt .or. movhvy) then
            chiave(inear,imixed) = chitot(inear,imixed) / nfqatm
            chiave(ifar,imixed) = chitot(ifar,imixed) / nfqatm
         endif
         if (movhvy) then
            chiave(inear,iheavy) = chitot(inear,iheavy) / nfqatm
            chiave(ifar,iheavy) = chitot(ifar,iheavy) / nfqatm
         endif
         do 230 imol = 1, nmol
            imty = molty(imol)
            do 220 iqind = 1, nfqmol(imty)
               iatom = iatmol(imol,ifqmol(imty,iqind))
               if (getnr) then
                  if (movlt) then
                     qforce(iatom,inear,ilight) = chiave(inear,ilight) -
     $                    chi(iatom,inear,ilight)
                  endif
                  if (movlt .or. movhvy) then
                     qforce(iatom,inear,imixed) = chiave(inear,imixed) -
     $                    chi(iatom,inear,imixed)
                  endif
                  if (movhvy) then
                     qforce(iatom,inear,iheavy) = chiave(inear,iheavy) -
     $                    chi(iatom,inear,iheavy)
                  endif
               endif
               if (getfar) then
                  if (movlt) then
                     qforce(iatom,ifar,ilight) = chiave(ifar,ilight) -
     $                    chi(iatom,ifar,ilight)
                  endif
                  if (movlt .or. movhvy) then
                     qforce(iatom,ifar,imixed) = chiave(ifar,imixed) -
     $                    chi(iatom,ifar,imixed)
                  endif
                  if (movhvy) then
                     qforce(iatom,ifar,iheavy) = chiave(ifar,iheavy) -
     $                    chi(iatom,ifar,iheavy)
                  endif
               endif
 220        continue
 230     continue
      else
         do 320 imol = 1, nmol
            imty = molty(imol)
            if (hasfqt(imty)) then
               if (movlt) then
                  chitot(inear,ilight) = 0.d0
                  chitot(ifar,ilight) = 0.d0
               endif
               if (movlt .or. movhvy) then
                  chitot(inear,imixed) = 0.d0
                  chitot(ifar,imixed) = 0.d0
               endif
               if (movhvy) then
                  chitot(inear,iheavy) = 0.d0
                  chitot(ifar,iheavy) = 0.d0
               endif
               do 300 iqind = 1, nfqmol(imty)
                  iatom = iatmol(imol,ifqmol(imty,iqind))
                  if (getnr) then
                     if (movlt) then
                        chitot(inear,ilight) = chitot(inear,ilight) +
     $                       chi(iatom,inear,ilight)
                     endif
                     if (movlt .or. movhvy) then
                        chitot(inear,imixed) = chitot(inear,imixed) +
     $                       chi(iatom,inear,imixed)
                     endif
                     if (movhvy) then
                        chitot(inear,iheavy) = chitot(inear,iheavy) +
     $                       chi(iatom,inear,iheavy)
                     endif
                  endif
                  if (getfar) then
                     if (movlt) then
                        chitot(ifar,ilight) = chitot(ifar,ilight) +
     $                       chi(iatom,ifar,ilight)
                     endif
                     if (movlt .or. movhvy) then
                        chitot(ifar,imixed) = chitot(ifar,imixed) +
     $                       chi(iatom,ifar,imixed)
                     endif
                     if (movhvy) then
                        chitot(ifar,iheavy) = chitot(ifar,iheavy) +
     $                       chi(iatom,ifar,iheavy)
                     endif
                  endif
 300           continue
               if (movlt) then
                  chiave(inear,ilight) = chitot(inear,ilight) / 
     $                 nfqmol(imty)
                  chiave(ifar,ilight) = chitot(ifar,ilight) / 
     $                 nfqmol(imty)
               endif
               if (movlt .or. movhvy) then
                  chiave(inear,imixed) = chitot(inear,imixed) / 
     $                 nfqmol(imty)
                  chiave(ifar,imixed) = chitot(ifar,imixed) /
     $                 nfqmol(imty)
               endif
               if (movhvy) then
                  chiave(inear,iheavy) = chitot(inear,iheavy) /
     $                 nfqmol(imty)
                  chiave(ifar,iheavy) = chitot(ifar,iheavy) / 
     $                 nfqmol(imty)
               endif
               do 310 iqind = 1, nfqmol(imty)
                  iatom = iatmol(imol,ifqmol(imty,iqind))
                  if (getnr) then
                     if (movlt) then
                        qforce(iatom,inear,ilight) = 
     $                       chiave(inear,ilight) -
     $                       chi(iatom,inear,ilight)
                     endif
                     if (movlt .or. movhvy) then
                        qforce(iatom,inear,imixed) = 
     $                       chiave(inear,imixed) -
     $                       chi(iatom,inear,imixed)
                     endif
                     if (movhvy) then
                        qforce(iatom,inear,iheavy) = 
     $                       chiave(inear,iheavy) -
     $                       chi(iatom,inear,iheavy)
                     endif
                  endif
                  if (getfar) then
                     if (movlt) then
                        qforce(iatom,ifar,ilight) = 
     $                       chiave(ifar,ilight) -
     $                       chi(iatom,ifar,ilight)
                     endif
                     if (movlt .or. movhvy) then
                        qforce(iatom,ifar,imixed) = 
     $                       chiave(ifar,imixed) -
     $                       chi(iatom,ifar,imixed)
                     endif
                     if (movhvy) then
                        qforce(iatom,ifar,iheavy) = 
     $                       chiave(ifar,iheavy) -
     $                       chi(iatom,ifar,iheavy)
                     endif
                  endif
 310           continue
            endif
 320     continue
      endif
c$$$      if (getnr) then
c$$$         if (ctflag) then
c$$$            chitot = 0.d0
c$$$            do 210 imol = 1, nmol
c$$$               imty = molty(imol)
c$$$               do 200 iqind = 1, nfqmol(imty)
c$$$                  iatom = iatmol(imol,ifqmol(imty,iqind))
c$$$                  chitot = chitot + chinr(iatom)
c$$$ 200           continue
c$$$ 210        continue
c$$$            chiave = chitot / nfqatm
c$$$            do 230 imol = 1, nmol
c$$$               imty = molty(imol)
c$$$               do 220 iqind = 1, nfqmol(imty)
c$$$                  iatom = iatmol(imol,ifqmol(imty,iqind))
c$$$                  qforcn(iatom) = chiave - chinr(iatom)
c$$$ 220           continue
c$$$ 230        continue
c$$$         else
c$$$            do 320 imol = 1, nmol
c$$$               imty = molty(imol)
c$$$               if (hasfqt(imty)) then
c$$$                  chitot = 0.d0
c$$$                  do 300 iqind = 1, nfqmol(imty)
c$$$                     iatom = iatmol(imol,ifqmol(imty,iqind))
c$$$                     chitot = chitot + chinr(iatom)
c$$$ 300              continue
c$$$                  chiave = chitot / nfqmol(imty)
c$$$                  do 310 iqind = 1, nfqmol(imty)
c$$$                     iatom = iatmol(imol,ifqmol(imty,iqind))
c$$$                     qforcn(iatom) = chiave - chinr(iatom)
c$$$ 310              continue
c$$$               endif
c$$$ 320        continue
c$$$         endif
c$$$      endif
c$$$      if (getfar) then
c$$$         if (ctflag) then
c$$$            chitot = 0.d0
c$$$            do 410 imol = 1, nmol
c$$$               imty = molty(imol)
c$$$               do 400 iqind = 1, nfqmol(imty)
c$$$                  iatom = iatmol(imol,ifqmol(imty,iqind))
c$$$                  chitot = chitot + chifar(iatom)
c$$$ 400           continue
c$$$ 410        continue
c$$$            chiave = chitot / nfqatm
c$$$            do 430 imol = 1, nmol
c$$$               imty = molty(imol)
c$$$               do 420 iqind = 1, nfqmol(imty)
c$$$                  iatom = iatmol(imol,ifqmol(imty,iqind))
c$$$                  qforcf(iatom) = chiave - chifar(iatom)
c$$$ 420           continue
c$$$ 430        continue
c$$$         else
c$$$            do 520 imol = 1, nmol
c$$$               imty = molty(imol)
c$$$               if (hasfqt(imty)) then
c$$$                  chitot = 0.d0
c$$$                  do 500 iqind = 1, nfqmol(imty)
c$$$                     iatom = iatmol(imol,ifqmol(imty,iqind))
c$$$                     chitot = chitot + chifar(iatom)
c$$$ 500              continue
c$$$                  chiave = chitot / nfqmol(imty)
c$$$                  do 510 iqind = 1, nfqmol(imty)
c$$$                     iatom = iatmol(imol,ifqmol(imty,iqind))
c$$$                     qforcf(iatom) = chiave - chifar(iatom)
c$$$ 510              continue
c$$$               endif
c$$$ 520        continue
c$$$         endif
c$$$      endif
      return
      end
c
