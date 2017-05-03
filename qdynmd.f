c***********************************************************************
c     filvel fills the velocity vector with Maxwell-distributed
c     velocities.  Each dimension is shifted to insure a zero total
c     momentum...sjs 6/3/92
c***********************************************************************
c     Now done by molecules....sjs 6/9/94
c***********************************************************************
      subroutine filvel(tmp, jseed)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c     
      logical lspare
      real*8  tmom(3)
c
      return
      tmom(1) = 0.0d0
      tmom(2) = 0.0d0
      tmom(3) = 0.0d0
c***********************************************************************
c     3 * nmol degrees of freedom are given 1/2 kT of KE each, but after
c     subtracting out system center-of-mass motion and equilibrating
c     into rotational degrees of freedom, this will end up distributing
c     into nrdof degrees of freedom, so tmp is scaled:
c***********************************************************************
      if (nmol .eq. 1) then
         do 95 idim = 1, 3
            do 90 imind = 1, numatm(1)
               vel(iatmol(1,imind),idim) = 0.0d0
 90         continue
 95      continue
         return
      endif
      tmpnew = tmp * nmol / (nmol - 1) * nrdof / (3 * nmol - 3)
      do 120 imol = 1, nmol
         imty = molty(imol)
         rmmass = 0.0d0
         do 100 imaind = 1, nmamol(imty)
            imind = imamol(imty,imaind)
            rmmass = rmmass + atmass(iatype(imty,imind))
 100     continue
         rmminv = 1.0d0 / rmmass
         do 110 idim = 1, 3
            if (lspare) then
               ranvel = sqrt(ranb1 * rmminv) * cos(ranb2)
               lspare = .false.
            else
               ranb1 = -2.0d0 * log(gran(jseed)) * k * tmpnew
               ranb2 = 2.0d0 * pi * gran(jseed)
               ranvel = sqrt(ranb1 * rmminv) * sin(ranb2)
               lspare = .true.
            endif
            do 105 ifind = 1, nframe(imty)
               vel(iatmol(imol,ifind),idim) = ranvel
 105        continue
            tmom(idim) = tmom(idim) + ranvel * rmmass
 110     continue
 120  continue
      tmom(1) = tmom(1) / symass
      tmom(2) = tmom(2) / symass
      tmom(3) = tmom(3) / symass
      do 135 imol = 1, nmol
         imty = molty(imol)
         do 130 ifind = 1, nframe(imty)
            iatom = iatmol(imol,ifind)
            vel(iatom,1) = vel(iatom,1) - tmom(1)
            vel(iatom,2) = vel(iatom,2) - tmom(2)
            vel(iatom,3) = vel(iatom,3) - tmom(3)
 130     continue
 135  continue
      return
      end
c
c***********************************************************************
c     gettmp calculates the atomic and molecular temperatures of the 
c     system....sjs 1/18/93
c***********************************************************************
c
      subroutine gettmp()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
c***********************************************************************
c     Avoid division by zero.  nmdof will be zero for unimolecular
c     systems.  nrdof will be zero for (a) single-atom systems, or (b)
c     nonperiodic systems of a single rigid molecule.
c***********************************************************************
      if (nrdof .eq. 0) then
         rketmp = 0.d0
      endif
      if (nmdof .eq. 0) then
         rcmtmp = 0.d0
      endif
c$$$      if (nmol .eq. 1) then
c$$$         gettmp = 0
c$$$         return
c$$$      endif
c***********************************************************************
c     The atomic temperature is calculated from the KE of all of the 
c     real degrees of freedom.  Drude KE is removed, charge KE was
c     never included.
c***********************************************************************
      rtfac = 2.d0 / (nrdof * k)
c$$$      gettmp = ake * rtfac
      rketmp = (ake - dke) * rtfac
c***********************************************************************
c     The molecular temperature depends on the molecular KE, which uses
c     only center-of-mass velocities:
c***********************************************************************
      rmke = 0.d0
      do 120 imol = 1,nmol
         imty = molty(imol)
         cmpx = 0.d0
         cmpy = 0.d0
         cmpz = 0.d0
         rmolms = 0.d0
         do 100 imaind = 1, nmamol(imty)
            imind = imamol(imty,imaind)
            iaty = iatype(imty,imind)
            cmpx = cmpx + atmass(iaty) * vel(iatmol(imol,imind),1)
            cmpy = cmpy + atmass(iaty) * vel(iatmol(imol,imind),2)
            cmpz = cmpz + atmass(iaty) * vel(iatmol(imol,imind),3)
            rmolms = rmolms + atmass(iaty)
  100    continue
         rmke = rmke + 0.5d0 * (cmpx ** 2 + cmpy ** 2 + cmpz ** 2) /
     $        rmolms
  120 continue
      rtfac = 2.d0 / (nmdof * k)
      rcmtmp = rmke * rtfac
      return
      end
c
c***********************************************************************
c     rattlx does a velocity step of velocity Verlet
c     (v = v + dt / 2 * F), then a position step (x = x + dt * v) and
c     RATTLES the positions until they satisfy the bond constraints.
c     See Andersen, J Comp Phys 52, 24 (1983) for the RATTLE algorithm.
c     The subroutine is passed the timestep size, the force array, and
c     a mask for atom types, so that it may be used in generic RESPA
c     breakups....sjs 3/9/95
c***********************************************************************
c
      subroutine rattlx(dt, idist, imass)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
      logical okflag
      real*8  rij(3), rijold(3),
     $     oldpos(maxatm,3)
c
c***********************************************************************
c     Remember the old atom positions, because the constraint forces
c     will act along the old bonds:
c***********************************************************************
      do 210 imol = 1, nmol
         imty = molty(imol)
         do 200 imind = 1, nframe(imty)
            iaty = iatype(imty,imind)
            if (masmsk(iaty,imass)) then
               iatom = iatmol(imol,imind)
               oldpos(iatom,1) = pos(iatom,1)
               oldpos(iatom,2) = pos(iatom,2)
               oldpos(iatom,3) = pos(iatom,3)
            endif
 200     continue
 210  continue
c***********************************************************************
c     Update the velocities for a half-timestep and the positions for
c     a full timestep using the unconstrained forces:
c***********************************************************************
      call vverlv(dt, idist, imass)
      call vverlx(dt, imass)
c***********************************************************************
c     Loop over molecules.  For each one, loop over constraints, fixing
c     the velocities and positions with the constraint forces until
c     everything is okay.
c***********************************************************************
      do 420 imol = 1, nmol
         imty = molty(imol)
         iriter = 0
 400     continue
         okflag = .true.
         iriter = iriter + 1
         do 410 iconi = 1, ncons(imty)
            imind1 = icons(imty,iconi,1)
            imind2 = icons(imty,iconi,2)
            iaty1 = iatype(imty,imind1)
            iaty2 = iatype(imty,imind2)
            if (masmsk(iaty1,imass) .and. masmsk(iaty2,imass)) then
               iatom1 = iatmol(imol,imind1)
               iatom2 = iatmol(imol,imind2)
               rij(1) = pos(iatom1,1) - pos(iatom2,1)
               rij(2) = pos(iatom1,2) - pos(iatom2,2)
               rij(3) = pos(iatom1,3) - pos(iatom2,3)
               tr2dif = rij(1) * rij(1) + rij(2) * rij(2) +
     $              rij(3) * rij(3) - r2cons(imty,iconi)
               if (dabs(tr2dif) .gt. rattol) then
                  okflag = .false.
                  rijold(1) = oldpos(iatom1,1) - oldpos(iatom2,1)
                  rijold(2) = oldpos(iatom1,2) - oldpos(iatom2,2)
                  rijold(3) = oldpos(iatom1,3) - oldpos(iatom2,3)
                  temp = rij(1) * rijold(1) + rij(2) * rijold(2) +
     $                 rij(3) * rijold(3)
                  scalpc = tr2dif / (2.d0 * dt * 
     $                 (atmass(iaty1) + atmass(iaty2)) * temp)
                  temp = scalpc * atmass(iaty2)
                  vel(iatom1,1) = vel(iatom1,1) - temp * rijold(1)
                  vel(iatom1,2) = vel(iatom1,2) - temp * rijold(2)
                  vel(iatom1,3) = vel(iatom1,3) - temp * rijold(3)
                  temp = temp * dt
                  pos(iatom1,1) = pos(iatom1,1) - temp * rijold(1)
                  pos(iatom1,2) = pos(iatom1,2) - temp * rijold(2)
                  pos(iatom1,3) = pos(iatom1,3) - temp * rijold(3)
                  temp = scalpc * atmass(iaty1)
                  vel(iatom2,1) = vel(iatom2,1) + temp * rijold(1)
                  vel(iatom2,2) = vel(iatom2,2) + temp * rijold(2)
                  vel(iatom2,3) = vel(iatom2,3) + temp * rijold(3)
                  temp = temp * dt
                  pos(iatom2,1) = pos(iatom2,1) + temp * rijold(1)
                  pos(iatom2,2) = pos(iatom2,2) + temp * rijold(2)
                  pos(iatom2,3) = pos(iatom2,3) + temp * rijold(3)
c***********************************************************************
c     This pressure virial is off by a factor of two from the one
c     calculated in the old rattl1 subroutine.  I don't know which is
c     wrong...
c***********************************************************************
                  if (ioiter) then
                     wvirrx = wvirrx + 2.d0 * atmass(iaty1) * 
     $                    atmass(iaty2) * scalpc / dt * 
     $                    (rijold(1) * rijold(1) + 
     $                    rijold(2) * rijold(2) + rijold(3) * rijold(3))
                  endif
               endif
            endif
 410     continue
         if (.not. okflag) then
            go to 400
         endif
 420  continue
      call tagup(imass)
      call foldr(imass)
      return
      end
c
c***********************************************************************
c     rattlv does a velocity step of velocity Verlet 
c     (v = v + dt / 2 * F) and RATTLEs the velocities until they
c     satisfy the (time derivative of the) bond constraints.  See
c     Andersen, J Comp Phys 52, 24 (1983) for the RATTLE algorithm.
c     ...sjs 3/9/95
c***********************************************************************
c
      subroutine rattlv(dt, idist, imass)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
      logical okflag
      real*8  rij(3)
c
c***********************************************************************
c     Update the velocities for a half-timestep using the unconstrained
c     forces:
c***********************************************************************
      call vverlv(dt, idist, imass)
c***********************************************************************
c     Loop over molecules.  For each one, loop over the constraints, 
c     fixing the velocities with the constraint forces until everything
c     is okay:
c***********************************************************************
      do 320 imol = 1, nmol
         imty = molty(imol)
 300     continue
         okflag = .true.
         do 310 iconi = 1, ncons(imty)
            iconi1 = icons(imty,iconi,1)
            iconi2 = icons(imty,iconi,2)
            iaty1 = iatype(imty,iconi1)
            iaty2 = iatype(imty,iconi2)
            if (masmsk(iaty1,imass) .and. masmsk(iaty2,imass)) then
               iatom1 = iatmol(imol,iconi1)
               iatom2 = iatmol(imol,iconi2)
               rij(1) = pos(iatom1,1) - pos(iatom2,1)
               rij(2) = pos(iatom1,2) - pos(iatom2,2)
               rij(3) = pos(iatom1,3) - pos(iatom2,3)
               trdotv = rij(1) * (vel(iatom1,1) - vel(iatom2,1)) +
     $              rij(2) * (vel(iatom1,2) - vel(iatom2,2)) +
     $              rij(3) * (vel(iatom1,3) - vel(iatom2,3))
               if (abs(trdotv) .gt. rattol) then
                  okflag = .false.
                  scalpc = trdotv / ((atmass(iaty1) + atmass(iaty2)) * 
     $                 r2cons(imty,iconi))
                  temp = scalpc * atmass(iaty2)
                  vel(iatom1,1) = vel(iatom1,1) - temp * rij(1)
                  vel(iatom1,2) = vel(iatom1,2) - temp * rij(2)
                  vel(iatom1,3) = vel(iatom1,3) - temp * rij(3)
                  temp = scalpc * atmass(iaty1)
                  vel(iatom2,1) = vel(iatom2,1) + temp * rij(1)
                  vel(iatom2,2) = vel(iatom2,2) + temp * rij(2)
                  vel(iatom2,3) = vel(iatom2,3) + temp * rij(3)
c***********************************************************************
c     This pressure virial is off by a factor of two from the one
c     calculated in the old rattl2 subroutine.  I don't know which is
c     wrong...
c***********************************************************************
                  if (ioflg .and. 
     $                 nioint * ((nciter + 1) / nioint) .eq. nciter + 1)
     $                 then
                     wvirrv = wvirrv - 2.d0 * atmass(iaty1) * 
     $                    atmass(iaty2) * scalpc / dt * 
     $                    (rij(1) * rij(1) + rij(2) * rij(2) + 
     $                    rij(3) * rij(3))
                  endif
               endif
            endif
 310     continue
         if (.not. okflag) then
            go to 300
         endif
 320  continue
      return
      end
c
c***********************************************************************
c     foldr folds the pos() array of positions back into the main box,
c     storing them in fldpos().  All atoms in a given molecule must be 
c     in the same box, however....sjs 9/14/93
c***********************************************************************
c     changed so that the primary box is now (-L/2,L/2) instead of 
c     (0,L).  This is more convenient for slab calcs, so that the center
c     of the slab can be at z=0....sjs 4/23/96
c***********************************************************************
c
      subroutine foldr(imass)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      do 210 imol = 1, nmol
         imty = molty(imol)
         iaty = iatype(imty,1)
         if (masmsk(iaty,imass)) then
            ioxy = iatmol(imol,1)
            dx = anint(pos(ioxy,1) / boxsiz(1)) * boxsiz(1)
            dy = anint(pos(ioxy,2) / boxsiz(2)) * boxsiz(2)
            dz = anint(pos(ioxy,3) / boxsiz(3)) * boxsiz(3)
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
c     vverlv does the velocity part of velocity Verlet, updating the
c     velocities by half a timestep....sjs 2/20/95
c***********************************************************************
c
      subroutine vverlv(dt, idist, imass)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
      dtby2 = 0.5d0 * dt
      do 210 imol = 1, nmol
         imty = molty(imol)
         do 200 ifind = 1, nframe(imty)
            iaty = iatype(imty,ifind)
            if (masmsk(iaty,imass)) then
               iatom = iatmol(imol,ifind)
               vel(iatom,1) = vel(iatom,1) + dtby2 / atmass(iaty) *
     $              (force(iatom,1,idist,imass) + 
     $              force(iatom,1,idist,imixed))
               vel(iatom,2) = vel(iatom,2) + dtby2 / atmass(iaty) *
     $              (force(iatom,2,idist,imass) +
     $              force(iatom,2,idist,imixed))
               vel(iatom,3) = vel(iatom,3) + dtby2 / atmass(iaty) *
     $              (force(iatom,3,idist,imass) +
     $              force(iatom,3,idist,imixed))
            endif
 200     continue
 210  continue
      return
      end 
c
c***********************************************************************
c     vverlx does the position part of velocity Verlet, updating the
c     positions by a full timestep using the half-step velocities.  The
c     timestep size and atom mask are passed in as an arguments so that 
c     this routine can be used in various types of RESPA breakups.
c     ...sjs 2/20/95
c***********************************************************************
c
      subroutine vverlx(dt, imass)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
      do 210 imol = 1, nmol
         imty = molty(imol)
         do 200 ifind = 1, nframe(imty)
            iaty = iatype(imty,ifind)
            if (masmsk(iaty,imass)) then
               iatom = iatmol(imol,ifind)
               pos(iatom,1) = pos(iatom,1) + dt * vel(iatom,1)
               pos(iatom,2) = pos(iatom,2) + dt * vel(iatom,2)
               pos(iatom,3) = pos(iatom,3) + dt * vel(iatom,3)
            endif
 200     continue
 210  continue
c***********************************************************************
c     pull the tag atoms up to the frame, and fold everybody back into
c     the simulation box:
c***********************************************************************
      call tagup(imass)
      call foldr(imass)
c***********************************************************************
c     set the flags to tell other subroutines which atoms have moved:
c***********************************************************************
      if (imass .eq. ilight) then
         movlt = .true.
      else if (imass .eq. imixed) then
         movlt = .true.
         movhvy = .true.
      else if (imass .eq. iheavy) then
         movhvy = .true.
      endif
      return
      end
c
c***********************************************************************
c     qvvrlv does the velocity part of velocity Verlet for the charge 
c     degrees of freedom, updating the velocities by half a timestep.
c     ...sjs 2/20/95
c***********************************************************************
c
      subroutine qvvrlv(dt, idist, imass)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
      dtby2 = 0.5d0 * dt
      do 210 imol = 1, nmol
         imty = molty(imol)
         do 200 iqind = 1, nfqmol(imty)
            imind = ifqmol(imty,iqind)
            iaty = iatype(imty,imind)
            if (masmsk(iaty,imass)) then
               iatom = iatmol(imol,imind)
               qvel(iatom) = qvel(iatom) + dtby2 * qminv * 
     $              (qforce(iatom,idist,imass) + 
     $              qforce(iatom,idist,imixed))
            endif
 200     continue
 210  continue
      return
      end
c
c***********************************************************************
c     qvvrlx does the position part of velocity Verlet for the charge
c     degrees of freedom, updating the positions by a full timestep.
c     ...sjs 2/20/95
c***********************************************************************
c
      subroutine qvvrlx(dt, imass)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
      do 210 imol = 1, nmol
         imty = molty(imol)
         do 200 iqind = 1, nfqmol(imty)
            imind = ifqmol(imty,iqind)
            iaty = iatype(imty,imind)
            if (masmsk(iaty,imass)) then
               iatom = iatmol(imol,imind)
               q(iatom) = q(iatom) + dt * qvel(iatom)
            endif
 200     continue
 210  continue
c***********************************************************************
c     some charges have changed, so the atoms have "moved"
c***********************************************************************
      if (imass .eq. ilight) then
         movlt = .true.
      else if (imass .eq. imixed) then
         movlt = .true.
         movhvy = .true.
      else if (imass .eq. iheavy) then
         movhvy = .true.
      endif
      return
      end
c
c***********************************************************************
c     brmono builds a bunch of molecule-sized position arrays (r, r2,
c     1/r, dx, dy, dz) for intramolecular distances in a given molecule.
c     ...sjs 6/17/94
c***********************************************************************
c     try rewriting this with BLAS routines.
c***********************************************************************
c
      subroutine brmono(imol)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      do 210 imind1 = 1, numatm(imol)
         iatom1 = iatmol(imol,imind1)
         do 200 imind2 = 1, imind1 - 1
            iatom2 = iatmol(imol,imind2)
            dxmat(imind1,imind2) = pos(iatom2,1) - pos(iatom1,1)
            dxmat(imind2,imind1) = -dxmat(imind1,imind2)
            dymat(imind1,imind2) = pos(iatom2,2) - pos(iatom1,2)
            dymat(imind2,imind1) = -dymat(imind1,imind2)
            dzmat(imind1,imind2) = pos(iatom2,3) - pos(iatom1,3)
            dzmat(imind2,imind1) = -dzmat(imind1,imind2)
            r2mat(imind1,imind2) = 
     $           dxmat(imind1,imind2) * dxmat(imind1,imind2) +
     $           dymat(imind1,imind2) * dymat(imind1,imind2) +
     $           dzmat(imind1,imind2) * dzmat(imind1,imind2)
            r2mat(imind2,imind1) = r2mat(imind1,imind2)
            rmat(imind1,imind2) = sqrt(r2mat(imind1,imind2))
            rmat(imind2,imind1) = rmat(imind1,imind2)
            rinmat(imind1,imind2) = 1.0d0 / rmat(imind1,imind2)
            rinmat(imind2,imind1) = rinmat(imind1,imind2)
 200     continue
         dxmat(imind1,imind1) = 0.0d0
         dymat(imind1,imind1) = 0.0d0
         dzmat(imind1,imind1) = 0.0d0
         r2mat(imind1,imind1) = 0.0d0
         rmat(imind1,imind2) = 0.0d0
         rinmat(imind1,imind2) = 0.0d0
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
c
      subroutine brhead(imol1, imol2, dxinc, dyinc, dzinc)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      iatom1 = iatmol(imol1,1)
      iatom2 = iatmol(imol2,1)
      dxmat(1,1) = pos(iatom2,1) - pos(iatom1,1)
      dymat(1,1) = pos(iatom2,2) - pos(iatom1,2)
      dzmat(1,1) = pos(iatom2,3) - pos(iatom1,3)
      if (pbflag) then
         dxinc = -1 * dnint(dxmat(1,1) / boxsiz(1)) * boxsiz(1)
         dyinc = -1 * dnint(dymat(1,1) / boxsiz(2)) * boxsiz(2)
         dzinc = -1 * dnint(dzmat(1,1) / boxsiz(3)) * boxsiz(3)
      else
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
c
      subroutine brpair(imol1, imol2, dxinc, dyinc, dzinc)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c     
      imty1 = molty(imol1)
      imty2 = molty(imol2)
      iatom1 = iatmol(imol1,1)
      iaty1 = iatype(imty1,1)
      do 200 imind2 = 2, numatm(imol2)
         iatom2 = iatmol(imol2,imind2)
         iaty2 = iatype(imty2,imind2)
         if ((isqat(iaty1) .and. isqat(iaty2)) .or. 
     $        (isljat(iaty1) .and. isljat(iaty2))) then
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
            if ((isqat(iaty1) .and. isqat(iaty2)) .or.
     $           (isljat(iaty1) .and. isljat(iaty2))) then
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
c     fdist distributes the forces on all of the tag atoms onto the
c     frame....sjs 6/15/94
c***********************************************************************
c
      subroutine fdist(getlt, gethvy, getnr, getfar)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      logical getlt, gethvy, getnr, getfar
c
      do 220 imol = 1, nmol
         imty = molty(imol)
         do 210 itind = nframe(imty) + 1, nframe(imty) + ntag(imty)
            iatom = iatmol(imol,itind)
            do 200 ifind = 1, nframe(imty)
               iatom2 = iatmol(imol,ifind)
               if (getlt .and. movlt) then
                  if (getnr) then
                     force(iatom2,1,inear,ilight) = 
     $                    force(iatom2,1,inear,ilight) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,1,inear,ilight)
                     force(iatom2,2,inear,ilight) = 
     $                    force(iatom2,2,inear,ilight) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,2,inear,ilight)
                     force(iatom2,3,inear,ilight) = 
     $                    force(iatom2,3,inear,ilight) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,3,inear,ilight)
                  endif
                  if (getfar) then
                     force(iatom2,1,ifar,ilight) = 
     $                    force(iatom2,1,ifar,ilight) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,1,ifar,ilight)
                     force(iatom2,2,ifar,ilight) = 
     $                    force(iatom2,2,ifar,ilight) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,2,ifar,ilight)
                     force(iatom2,3,ifar,ilight) = 
     $                    force(iatom2,3,ifar,ilight) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,3,ifar,ilight)
                  endif
               endif
               if ((getlt .or. gethvy) .and. (movlt .or. movhvy)) then
                  if (getnr) then
                     force(iatom2,1,inear,imixed) = 
     $                    force(iatom2,1,inear,imixed) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,1,inear,imixed)
                     force(iatom2,2,inear,imixed) = 
     $                    force(iatom2,2,inear,imixed) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,2,inear,imixed)
                     force(iatom2,3,inear,imixed) = 
     $                    force(iatom2,3,inear,imixed) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,3,inear,imixed)
                  endif
                  if (getfar) then
                     force(iatom2,1,ifar,imixed) = 
     $                    force(iatom2,1,ifar,imixed) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,1,ifar,imixed)
                     force(iatom2,2,ifar,imixed) = 
     $                    force(iatom2,2,ifar,imixed) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,2,ifar,imixed)
                     force(iatom2,3,ifar,imixed) = 
     $                    force(iatom2,3,ifar,imixed) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,3,ifar,imixed)
                  endif
               endif
               if (gethvy .and. movhvy) then
                  if (getnr) then
                     force(iatom2,1,inear,iheavy) = 
     $                    force(iatom2,1,inear,iheavy) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,1,inear,iheavy)
                     force(iatom2,2,inear,iheavy) = 
     $                    force(iatom2,2,inear,iheavy) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,2,inear,iheavy)
                     force(iatom2,3,inear,iheavy) = 
     $                    force(iatom2,3,inear,iheavy) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,3,inear,iheavy)
                  endif
                  if (getfar) then
                     force(iatom2,1,ifar,iheavy) = 
     $                    force(iatom2,1,ifar,iheavy) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,1,ifar,iheavy)
                     force(iatom2,2,ifar,iheavy) = 
     $                    force(iatom2,2,ifar,iheavy) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,2,ifar,iheavy)
                     force(iatom2,3,ifar,iheavy) = 
     $                    force(iatom2,3,ifar,iheavy) + 
     $                    tagcof(imty,itind,ifind) * 
     $                    force(iatom,3,ifar,iheavy)
                  endif
               endif
 200        continue
 210     continue
 220  continue
      return
      end
c
c***********************************************************************
c     tagup moves all the tag atoms to their proper place, relative to
c     the frame atoms....sjs 6/16/94
c***********************************************************************
c     
      subroutine tagup(imass)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      do 220 imol = 1, nmol
         imty = molty(imol)
         do 210 itind = nframe(imty) + 1, nframe(imty) + ntag(imty)
            if (masmsk(iatype(imty,itind),imass)) then
               iatom = iatmol(imol,itind)
               iatom1 = iatmol(imol,1)
               pos(iatom,1) = pos(iatom1,1)
               pos(iatom,2) = pos(iatom1,2)
               pos(iatom,3) = pos(iatom1,3)
               do 200 ifind = 2, nframe(molty(imol))
                  iatom2 = iatmol(imol,ifind)
                  pos(iatom,1) = pos(iatom,1) + 
     $                 tagcof(imty,itind,ifind) * 
     $                 (pos(iatom2,1) - pos(iatom1,1))
                  pos(iatom,2) = pos(iatom,2) + 
     $                 tagcof(imty,itind,ifind) *
     $                 (pos(iatom2,2) - pos(iatom1,2))
                  pos(iatom,3) = pos(iatom,3) + 
     $                 tagcof(imty,itind,ifind) *
     $                 (pos(iatom2,3) - pos(iatom1,3))
 200           continue
            endif
 210     continue
 220  continue
      return
      end
c
c***********************************************************************
c     getcom finds the center of mass of the cluster (should only be
c     called when not using periodic b.c.)....sjs 2/15/95
c***********************************************************************
c
      subroutine getcom()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
      rcom(1) = 0.0d0
      rcom(2) = 0.0d0
      rcom(3) = 0.0d0
c***********************************************************************
c     find the center of mass and store it in rcom:
c***********************************************************************
      do 210 imol = 1, nmol
         imty = molty(imol)
         do 200 imaind = 1, nmamol(imty)
            iatom = iatmol(imol,imamol(imty,imaind))
            rcom(1) = rcom(1) + pos(iatom,1)
            rcom(2) = rcom(2) + pos(iatom,2)
            rcom(3) = rcom(3) + pos(iatom,3)
 200     continue
 210  continue
      rcom(1) = rcom(1) / symass
      rcom(2) = rcom(2) / symass
      rcom(3) = rcom(3) / symass
      continue
      end
c
