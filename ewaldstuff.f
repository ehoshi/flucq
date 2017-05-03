c***********************************************************************
c     getkvc sets up arrays of k-vectors and associated calculations to
c     be used in Ewald summation....sjs 8/30/93
c***********************************************************************
c
      subroutine getkvc()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
c***********************************************************************
c     Initialize some stuff:
c***********************************************************************
      ewkmx2 = ewlkmx * ewlkmx
      rkpfac = 0.25d0 / ewlkp2
      numkvc = 0
      twopi = 2.0d0 * pi
      rntok(1) = twopi / boxsiz(1)
      rntok(2) = twopi / boxsiz(2)
      rntok(3) = twopi / boxsiz(3)
      volfac = 4.0d0 * pi / vol
c***********************************************************************
c     Get the largest possible lattice vectors in each direction:
c***********************************************************************
      ixend = int(ewlkmx / rntok(1))
      if (cubflg) then
         iyend = ixend
         izend = ixend
         nbig = ixend
      else
         iyend = int(ewlkmx / rntok(2))
         izend = int(ewlkmx / rntok(3))
         nbig = max(ixend,iyend,izend)
         if (nbig .gt. maxk1d) then
            write(iuout, *) 'getkvc: largest k vector component is ',
     $           nbig, ' but maxk1d is only ', maxk1d
            stop
         endif
      endif
c***********************************************************************
c     Counting only one of any inversion-symmetric pair, find all
c     vectors smaller than ewlkmx, except the zero vector:
c***********************************************************************
      do 220 nx = 0, izend
         rkx = rntok(1) * nx
         rk21 = rkx * rkx
         if (nx .eq. 0) then
            iystrt = 0
         else
            iystrt = -iyend
         endif
         do 210 ny = iystrt, iyend
            rky = rntok(2) * ny
            rk22 = rk21 + rky * rky
            if (rk22 .gt. ewkmx2) go to 210
            if (ny .eq. 0 .and. nx .eq. 0) then
               izstrt = 1
            else
               izstrt = -izend
            endif
            do 200 nz = izstrt, izend
               rkz = rntok(3) * nz
               rk23 = rk22 + rkz * rkz
               if (rk23 .gt. ewkmx2) go to 200
               numkvc = numkvc + 1
               if (numkvc .gt. maxkvc) then
                  write(iuout, *) 'getkvc: more than ', maxkvc, 
     $                 ' k-vectors found'
                  stop
               endif
               nvec(numkvc,1) = nx
               nvec(numkvc,2) = ny
               nvec(numkvc,3) = nz
               rkvec(numkvc,1) = rkx
               rkvec(numkvc,2) = rky
               rkvec(numkvc,3) = rkz
               r1kfac(numkvc) = rk23 * rkpfac
               rkfac(numkvc) = volfac / rk23 * exp(-r1kfac(numkvc))
c$$$               rkfac(numkvc) = volfac / rk23 * exp(-rk23 * rkpfac)
               r2kfac(numkvc) = 2.0d0 * rkfac(numkvc)
  200       continue
  210    continue
  220 continue
      write(iuout, '(1x,a,i4,a)') 'Using ', numkvc, ' k-vectors.'
      return
      end
c
c***********************************************************************
c     ewlsrf calculates the surface energy term in the Ewald sum, and
c     the forces on the charges and atoms that result from that term.
c     The subroutine dipols must be called before this one, to set up
c     the smu array....sjs 9/7/93
c     This routine should be called only for a very large Ewald sphere
c     surrounded by vacuum.  For a surrounding conductor (better for
c     dielectric constant calcs), do not call this routine.
c     ...sjs 9/30/93
c***********************************************************************
c     Fixed up so that it uses a mask of atom types, for use with
c     RESPA....sjs 3/16/95
c***********************************************************************
c     Needs a pressure calculation...sjs 6/17/96
c***********************************************************************
c
      subroutine ewlsrf(getlt, gethvy, gete)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      logical getlt, gethvy, gete

c
      real*8  fcpcv(3,nmass)
c
      pepc = srffac
c***********************************************************************
c     calculate the molecular and system dipoles for the masked atoms:
c***********************************************************************
      call dipols(getlt, gethvy)
c***********************************************************************
c     calculate the mass-split PE:
c***********************************************************************
      if (gete) then
         if (getlt .and. movlt) then
            qpe(ilight) = qpe(ilight) + pepc * 
     $           (smu(1,ilight) * smu(1,ilight) +
     $           smu(2,ilight) * smu(2,ilight) +
     $           smu(3,ilight) * smu(3,ilight))
         endif
         if ((getlt .or. gethvy) .and. (movlt .or. movhvy)) then
            qpe(imixed) = qpe(imixed) + pepc * 2.d0 * 
     $           (smu(1,ilight) * smu(1,iheavy) +
     $           smu(2,ilight) * smu(2,iheavy) +
     $           smu(3,ilight) * smu(3,iheavy))
         endif
         if (gethvy .and. movhvy) then
            qpe(iheavy) = qpe(iheavy) + pepc *
     $           (smu(1,iheavy) * smu(1,iheavy) +
     $           smu(2,iheavy) * smu(2,iheavy) +
     $           smu(3,iheavy) * smu(3,iheavy))
         endif
         if (dbflag) then
            if (getlt .and. movlt) then
               ebit(4,ilight) = ebit(4,ilight) + pepc * 
     $              (smu(1,ilight) * smu(1,ilight) +
     $              smu(2,ilight) * smu(2,ilight) +
     $              smu(3,ilight) * smu(3,ilight))
            endif
            if ((getlt .or. gethvy) .and. (movlt .or. movhvy)) then
               ebit(4,imixed) = ebit(4,imixed) + pepc * 2.d0 *
     $              (smu(1,ilight) * smu(1,iheavy) +
     $              smu(2,ilight) * smu(2,iheavy) +
     $              smu(3,ilight) * smu(3,iheavy))
            endif
            if (gethvy .and. movlt) then
               ebit(4,iheavy) = ebit(4,iheavy) + pepc *
     $              (smu(1,iheavy) * smu(1,iheavy) +
     $              smu(2,iheavy) * smu(2,iheavy) +
     $              smu(3,iheavy) * smu(3,iheavy))
            endif
         endif
      endif
cmov      if (bniter) then
cmov         do 100 imol1 = 1, nmol
cmov            qepre(imol1,imol1) = qepre(imol1,imol1) + pepc *
cmov     $           rmskmu(imol1,1) * rmu(imol1,1) +
cmov     $           rmskmu(imol1,2) * rmu(imol1,2) +
cmov     $           rmskmu(imol1,3) * rmu(imol1,3)
cmov            do 90 imol2 = 1, imol1 - 1
cmov               qepre(imol1,imol2) = qepre(imol1,imol2) + 2.d0 * pepc *
cmov     $              rmskmu(imol1,1) * rmu(imol2,1) +
cmov     $              rmskmu(imol1,2) * rmu(imol2,2) +
cmov     $              rmskmu(imol1,3) * rmu(imol2,3)
cmov 90         continue
cmov 100     continue
cmov      endif
c***********************************************************************
c     The system dipole is split by mass:
c***********************************************************************
      fcpc = pepc
      if (getlt .and. movlt) then
         fcpcv(1,ilight) = fcpc * smu(1,ilight)
         fcpcv(2,ilight) = fcpc * smu(2,ilight)
         fcpcv(3,ilight) = fcpc * smu(3,ilight)
      endif
      fcpcv(1,imixed) = fcpc * smu(1,imixed)
      fcpcv(2,imixed) = fcpc * smu(2,imixed)
      fcpcv(3,imixed) = fcpc * smu(3,imixed)
      if (gethvy .and. movhvy) then
         fcpcv(1,iheavy) = smu(1,imixed) - smu(1,ilight)
         fcpcv(2,iheavy) = smu(2,imixed) - smu(2,ilight)
         fcpcv(3,iheavy) = smu(3,imixed) - smu(3,ilight)
      endif
c***********************************************************************
c     Loop over atoms to calculate the forces and electric fields:
c***********************************************************************
      do 205 imol = 1, nmol
         imty = molty(imol)
         do 200 iqind = 1, nqmol(imty)
            imind = iqmol(imty,iqind)
            iaty = iatype(imty,imind)
            iatom = iatmol(imol,imind)
c***********************************************************************
c     The electric field is not split by mass:
c***********************************************************************
            if (dsiter) then
               field(iatom,1) = field(iatom,1) - fcpcv(1,imixed)
               field(iatom,2) = field(iatom,2) - fcpcv(2,imixed)
               field(iatom,3) = field(iatom,3) - fcpcv(3,imixed)
            endif
c***********************************************************************
c     Don't bother calculating pieces of the force array that we don't
c     care about or that haven't changed:
c***********************************************************************
            if (masmsk(iaty,ilight) .and. getlt .and. movlt) then
               force(iatom,1,ifar,ilight) = force(iatom,1,ifar,ilight) -
     $              q(iatom) * fcpcv(1,ilight)
               force(iatom,2,ifar,ilight) = force(iatom,2,ifar,ilight) -
     $              q(iatom) * fcpcv(2,ilight)
               force(iatom,3,ifar,ilight) = force(iatom,3,ifar,ilight) -
     $              q(iatom) * fcpcv(3,ilight)
            endif
            if ((getlt .or. gethvy) .and. (movlt .or. movhvy)) then
               if (masmsk(iaty,ilight)) then
                  imopp = iheavy
               else
                  imopp = ilight
               endif
               force(iatom,1,ifar,imixed) = force(iatom,1,ifar,imixed) -
     $              q(iatom) * fcpcv(1,imopp)
               force(iatom,2,ifar,imixed) = force(iatom,2,ifar,imixed) -
     $              q(iatom) * fcpcv(2,imopp)
               force(iatom,3,ifar,imixed) = force(iatom,3,ifar,imixed) -
     $              q(iatom) * fcpcv(3,imopp)
            endif
            if (masmsk(iaty,iheavy) .and. gethvy .and. movhvy) then
               force(iatom,1,ifar,iheavy) = force(iatom,1,ifar,iheavy) -
     $              q(iatom) * fcpcv(1,iheavy)
               force(iatom,2,ifar,iheavy) = force(iatom,2,ifar,iheavy) -
     $              q(iatom) * fcpcv(2,iheavy)
               force(iatom,3,ifar,iheavy) = force(iatom,3,ifar,iheavy) -
     $              q(iatom) * fcpcv(3,iheavy)
            endif
cmov            if (ioiter) then
cmov               wvir = wvir - q(iatom) *
cmov     $              (fldpos(iatom,1) * fcpcv(1,imixed) -
cmov     $              fldpos(iatom,3) * fcpcv(2,imixed) -
cmov     $              fldpos(iatom,2) * fcpcv(3,imixed))
cmov            endif
c$$$            if (ioiter) then
c$$$               wvir = wvir - fldpos(iatom,1) * fcx -
c$$$     $              fldpos(iatom,2) * fcy - fldpos(iatom,3) * fcz
c$$$            endif
 200     continue
 205  continue
c***********************************************************************
c     calculate the contribution to the electronegativity:
c***********************************************************************
      do 215 imol = 1, nmol
         imty = molty(imol)
         do 210 iqind = 1, nfqmol(imty)
            iatom = iatmol(imol,ifqmol(imty,iqind))
c***********************************************************************
c     Don't bother calculating pieces of the chi array that we don't
c     care about or that haven't changed:
c***********************************************************************
            if (masmsk(iaty,ilight) .and. getlt .and. movlt) then
               chi(iatom,ifar,ilight) = chi(iatom,ifar,ilight) +
     $              fldpos(iatom,1) * fcpcv(1,ilight) +
     $              fldpos(iatom,2) * fcpcv(2,ilight) +
     $              fldpos(iatom,3) * fcpcv(3,ilight)
            endif
            if ((getlt .or. gethvy) .and. (movlt .or. movhvy)) then
               if (masmsk(iaty,ilight)) then
                  imopp = iheavy
               else
                  imopp = ilight
               endif
               chi(iatom,ifar,imixed) = chi(iatom,ifar,imixed) +
     $              fldpos(iatom,1) * fcpcv(1,imopp) +
     $              fldpos(iatom,2) * fcpcv(2,imopp) +
     $              fldpos(iatom,3) * fcpcv(3,imopp)
            endif
            if (masmsk(iaty,iheavy) .and. gethvy .and. movhvy) then
               chi(iatom,ifar,iheavy) = chi(iatom,ifar,iheavy) +
     $              fldpos(iatom,1) * fcpcv(1,iheavy) +
     $              fldpos(iatom,2) * fcpcv(2,iheavy) +
     $              fldpos(iatom,3) * fcpcv(3,iheavy)
            endif
 210     continue
 215  continue
c***********************************************************************
c     calculate the charge Hessian for use in solving for the charges
c     explicitly.  this is always done with no mass masks, so there is
c     no mass split:
c***********************************************************************
      if (qsiter) then
         do 235 imol1 = 1, nmol
            imty1 = molty(imol1)
            do 230 iqind1 = 1, nqmol(imty1)
               imind1 = iqmol(imty1,iqind1)
               iaty1 = iatype(imty1,imind1)
               iatom1 = iatmol(imol1,imind1)
               qmat(iatom1,iatom1) = qmat(iatom1,iatom1) + pepc * 
     $              (fldpos(iatom1,1) * fldpos(iatom1,1) +
     $              fldpos(iatom1,2) * fldpos(iatom1,2) +
     $              fldpos(iatom1,3) * fldpos(iatom1,3))
               do 218 iqind2 = 1, iqind1 - 1
                  imind2 = iqmol(imty1,iqind2)
                  iaty2 = iatype(imty1,imind2)
                  iatom2 = iatmol(imol1,imind2)
                  qmatrm = pepc *
     $                 (fldpos(iatom1,1) * fldpos(iatom2,1) +
     $                 fldpos(iatom1,2) * fldpos(iatom2,2) +
     $                 fldpos(iatom1,3) * fldpos(iatom2,3))
                  qmat(iatom1,iatom2) = qmat(iatom1,iatom2) +
     $                 qmatrm
                  qmat(iatom2,iatom1) = qmat(iatom2,iatom1) +
     $                 qmatrm
 218           continue
               do 225 imol2 = 1, imol1 - 1
                  imty2 = molty(imol2)
                  do 220 iqind2 = 1, nqmol(imty2)
                     imind2 = iqmol(imty2,iqind2)
                     iaty2 = iatype(imty2,imind2)
                     iatom2 = iatmol(imol2,imind2)
                     qmatrm = pepc *
     $                    (fldpos(iatom1,1) * fldpos(iatom2,1) +
     $                    fldpos(iatom1,2) * fldpos(iatom2,2) +
     $                    fldpos(iatom1,3) * fldpos(iatom2,3))
                     qmat(iatom1,iatom2) = qmat(iatom1,iatom2) +
     $                    qmatrm
                     qmat(iatom2,iatom1) = qmat(iatom2,iatom1) +
     $                    qmatrm
 220              continue
 225           continue
 230        continue
 235     continue
      endif
      return
      end
c
c***********************************************************************
c     ewlslf calculates the self energy correction in the Ewald sum,
c     and the forces on the charges that result from that term.
c     ...sjs 9/7/93
c***********************************************************************
c
      subroutine ewlslf(getlt, gethvy, gete)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      logical getlt, gethvy, gete
c
c***********************************************************************
c     initialize some stuff:
c***********************************************************************
      sqsl = 0.d0
      sqsh = 0.d0
      sqs = 0.0d0
      consts = r2kbrp * epsinv
c***********************************************************************
c     calculate the sum of squared charges, split by mass:
c***********************************************************************
      do 205 imol = 1, nmol
         imty = molty(imol)
         do 200 iqind = 1, nqmol(imty)
            imind = iqmol(imty,iqind)
            iaty = iatype(imty,imind)
            iatom = iatmol(imol,imind)
            if (masmsk(iaty,ilight)) then
               sqsl = sqsl + q(iatom) * q(iatom)
            else
               sqsh = sqsh + q(iatom) * q(iatom)
            endif
 200     continue
 205  continue
      sqs = sqsl + sqsh
c***********************************************************************
c     calculate the PE:
c***********************************************************************
      if (gete) then
         if (getlt .and. movlt) then
            qpe(ilight) = qpe(ilight) - consts * sqsl
         endif
         if (gethvy .and. movhvy) then
            qpe(iheavy) = qpe(iheavy) - consts * sqsh
         endif
         if (dbflag) then
            if (getlt .and. movlt) then
               ebit(5,ilight) = ebit(5,ilight) - consts * sqsl
            endif
            if (gethvy .and. movhvy) then
               ebit(5,iheavy) = ebit(5,iheavy) - consts * sqsh
            endif
         endif
      endif
cmov      if (bniter) then
cmov         do 208 imol = 1, nmol
cmov            imty = molty(imol)
cmov            do 206 iqind = 1, nqmol(imty)
cmov               imind = iqmol(imty,iqind)
cmov               iaty = iatype(imty,imind)
cmov               if (mask(iaty)) then
cmov                  iatom = iatmol(imol,imind)
cmov                  qepre(imol,imol) = qepre(imol,imol) - 
cmov     $                 2.0d0 * consts * q(iatom) * q(iatom)
cmov               endif
cmov 206        continue
cmov 208     continue
cmov      endif
c***********************************************************************
c     calculate contributions to the electronegativity, dV/dq:
c***********************************************************************
      do 215 imol = 1, nmol
         imty = molty(imol)
         do 210 iqind = 1, nfqmol(imty)
            imind = ifqmol(imty,iqind)
            iaty = iatype(imty,imind)
            iatom = iatmol(imol,imind)
            if (masmsk(iaty,ilight) .and. getlt .and. movlt) then
               chi(iatom,ifar,ilight) = chi(iatom,ifar,ilight) -
     $              consts * q(iatom)
            else if (masmsk(iaty,iheavy) .and. gethvy .and. movhvy) then
               chi(iatom,ifar,iheavy) = chi(iatom,ifar,iheavy) -
     $              consts * q(iatom)
            endif
 210     continue
 215  continue
c***********************************************************************
c     calculate the charge Hessian for use in charge-solving.  This is
c     always done without mass masks:
c***********************************************************************
      if (qsiter) then
         do 225 imol = 1, nmol
            imty = molty(imol)
            do 220 iqind = 1, nqmol(imty)
               imind = iqmol(imty,iqind)
               iaty = iatype(imty,imind)
               iatom = iatmol(imol,imind)
               qmat(iatom,iatom) = qmat(iatom,iatom) - consts
 220        continue
 225     continue
      endif
      return
      end
c
c***********************************************************************
c     ewlrcp calculates the reciprocal-space term in the Ewald sum.
c     About how it works:  (acreal,acimag) contains the sum of 
c     exp(i k . r) for all atoms and a given k.  The complex 
c     exponentials are evaluated using the coskx and sinkx values set up
c     in gtcos()....sjs 9/10/93
c***********************************************************************
c     updated to use an atom type mask for compatibility with RESPA.
c     This only recalculates the forces and energies for terms that
c     involve a subset of the atoms, but does recalculate all forces
c     and energy terms that this subset participates in....sjs 3/14/95
c***********************************************************************
c     exp(i k r_i) is now stored in a huge array from one timestep to
c     the next, so that this routine is cheaper when the heavy atoms
c     have not moved since the last timestep....sjs 6/26/95
c***********************************************************************
c     Note:  this part of the Ewald potential includes some of the
c     central cell interaction (1/r - erfc(kr)/r).  This could be
c     short-ranged (for the RESPA split), even though this subroutine
c     should be called with the long-range force array.  The correction
c     for this is taken care of in qqjter....sjs 3/23/95
c***********************************************************************
c
      subroutine ewlrcp(getlt, gethvy, gete)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c     
      logical getlt, gethvy, gete
c
      logical yflip, zflip, yeqvz
c
c***********************************************************************
c     precalculate cos(n k . r) and sin(n k . r) for many n and all r:
c***********************************************************************
      call gtcos()
c***********************************************************************
c     loop over k-vectors, calculating exp(i k . r) from the cos and
c     sin components, storing them and also accumulating them in 
c     mass-split totals (only for those atoms that have moved):
c***********************************************************************
      do 300 ivec = 1, numkvc
         acreal = 0.d0
         acimag = 0.d0
         subacr = 0.d0
         subaci = 0.d0
         ackdre = 0.d0
         ackdim = 0.d0
         nxp1 = nvec(ivec,1) + 1
         if (nvec(ivec,2) .lt. 0) then
            yflip = .true.
            nyp1 = -nvec(ivec,2) + 1
         else
            yflip = .false.
             nyp1 = nvec(ivec,2) + 1
         endif
         if (nvec(ivec,3) .lt. 0) then
            zflip = .true.
            nzp1 = -nvec(ivec,3) + 1
         else
            zflip = .false.
            nzp1 = nvec(ivec,3) + 1
         endif
         yeqvz = yflip .eqv. zflip
         pepc = r2kfac(ivec) * epsinv
         do 253 imol = 1, nmol
            imty = molty(imol)
            do 250 iqind = 1, nqmol(imty)
               imind = iqmol(imty,iqind)
               iaty = iatype(imty,imind)
               iatom = iatmol(imol,imind)
               ihead = iatmol(imol,1)
c***********************************************************************
c     if this atom has moved, calculate exp(i k . r) and keep it in a 
c     big array.
c***********************************************************************
               if ((masmsk(iaty,ilight) .and. movlt) .or.
     $              (masmsk(iaty,iheavy) .and. movhvy)) then
                  if (yeqvz) then
                     ccss23 = 
     $                    coskx(iatom,nyp1,2) * coskx(iatom,nzp1,3) -
     $                    sinkx(iatom,nyp1,2) * sinkx(iatom,nzp1,3)
                     cssc23 = 
     $                    coskx(iatom,nyp1,2) * sinkx(iatom,nzp1,3) +
     $                    sinkx(iatom,nyp1,2) * coskx(iatom,nzp1,3)
                  else
                     ccss23 = 
     $                    coskx(iatom,nyp1,2) * coskx(iatom,nzp1,3) +
     $                    sinkx(iatom,nyp1,2) * sinkx(iatom,nzp1,3)
                     cssc23 = 
     $                    coskx(iatom,nyp1,2) * sinkx(iatom,nzp1,3) -
     $                    sinkx(iatom,nyp1,2) * coskx(iatom,nzp1,3)
                  endif
                  if (zflip) then
                     sgnds1 = sinkx(iatom,nxp1,1)
                     sgndc1 = -coskx(iatom,nxp1,1)
                  else
                     sgnds1 = -sinkx(iatom,nxp1,1)
                     sgndc1 = coskx(iatom,nxp1,1)
                  endif
                  eikrre(iatom,ivec) = coskx(iatom,nxp1,1) * ccss23 + 
     $                 sgnds1 * cssc23
                  eikrim(iatom,ivec) = sgndc1 * cssc23 +
     $                 sinkx(iatom,nxp1,1) * ccss23
               endif
c***********************************************************************
c     accumulate the sum of q * exp(i k . r) terms into mass-split
c     arrays:
c***********************************************************************
               qeikrr = q(iatom) * eikrre(iatom,ivec)
               qeikri = q(iatom) * eikrim(iatom,ivec)
               acreal = acreal + qeikrr
               acimag = acimag + qeikri
               if (masmsk(iaty,ilight)) then
                  subacr = subacr + qeikrr
                  subaci = subaci + qeikri
               endif
c***********************************************************************
c     accumulate q (k . d) exp(i k . r) for use in calculating the
c     pressure.  d is zero for a flexible atom and r-rhead for a rigid
c     one.
c***********************************************************************
               if (ioiter .and. gete) then
                  if (ncons(imty) .gt. 0) then
                     deex = fldpos(iatom,1) - fldpos(ihead,1)
                     deey = fldpos(iatom,2) - fldpos(ihead,2)
                     deez = fldpos(iatom,3) - fldpos(ihead,3)
c$$$                     tx = 0.d0
c$$$                     ty = 0.d0
c$$$                     tz = 0.d0
c$$$                     do 249 imaint = 1, nmamol(imty)
c$$$                        imindt = imamol(imty,imaint)
c$$$                        iatomt = iatmol(imol,imint)
c$$$                        tx = tx + fldpos(iatomt,1) * comcof(imty,imindt)
c$$$                        ty = ty + fldpos(iatomt,2) * comcof(imty,imindt)
c$$$                        tz = tz + fldpos(iatomt,3) * comcof(imty,imindt)
c$$$  249                continue
                     deex = fldpos(iatom,1) - tx
                     deey = fldpos(iatom,2) - ty
                     deez = fldpos(iatom,3) - tz
                     rkdee = deex * rkvec(ivec,1) + 
     $                    deey * rkvec(ivec,2) + deez * rkvec(ivec,3)
                     ackdre = ackdre + qeikrr * rkdee
                     ackdim = ackdim + qeikri * rkdee
                  endif
               endif
 250        continue
 253     continue
         supacr = acreal - subacr
         supaci = acimag - subaci
         scacre = pepc * acreal
         scacim = pepc * acimag
         scsubr = pepc * subacr
         scsubi = pepc * subaci
         scsupr = scacre - scsubr
         scsupi = scacim - scsubi
c***********************************************************************
c     calculate the mass-split PE:
c***********************************************************************
         if (gete) then
            if (getlt .and. movlt) then
               qpe(ilight) = qpe(ilight) + 
     $              subacr * scsubr + subaci * scsubi
            endif
            if ((getlt .or. gethvy) .and. (movlt .or. movhvy)) then
               qpe(imixed) = qpe(imixed) + 
     $              2.d0 * (subacr * scsupr + subaci * scsupi)
            endif
            if (gethvy .and. movhvy) then
               qpe(iheavy) = qpe(iheavy) + 
     $              supacr * scsupr + supaci * scsupi
            endif
            if (dbflag) then
               if (getlt .and. movlt) then
                  ebit(6,ilight) = ebit(6,ilight) + 
     $                 subacr * scsubr + subaci * scsubi
               endif
               if ((getlt .or. gethvy) .and. (movlt .or. movhvy)) then
                  ebit(6,imixed) = ebit(6,imixed) +
     $                 2.d0 * (subacr * scsupr + subaci * scsupi)
               endif
               if (gethvy .and. movhvy) then
                  ebit(6,iheavy) = ebit(6,iheavy) +
     $                 supacr * scsupr + supaci * scsupi
               endif
            endif
         endif
cmov         if (bniter) then
cmovc***********************************************************************
cmovc     for other pair energies E_ij = E_ji and only one (or half of both)
cmovc     is included in the system energy.  Here E_ij is the 
cmovc     interaction of atom/molecule i with j and all of its images. 
cmovc     E_ji is the interaction of j with i and all of its images, so it's
cmovc     not totally symmetric, but E_ij is still equal to E_ji, and half
cmovc     of the sum is included in the system energy.  The pair energy
cmovc     that gets put in qepre is one full interaction, the same as the
cmovc     contribution to the system PE.  For the E_ii term, what is put
cmovc     into qepre is again the full interaction energy, but this is twice
cmovc     what is put into the system PE.  This k-space part of the Ewald
cmovc     sum includes not only the image cell interactions, but also some
cmovc     part of the central cell, direct interactions (the part not 
cmovc     included when the direct interactions get multiplied by an erfc).
cmovc     This part of the central cell energy is double-counted, and this
cmovc     is taken care of when those interactions are added in.
cmovc***********************************************************************
cmov            do 258 imol1 = 1, nmol
cmov               imty1 = molty(imol1)
cmov               do 257 iqind1 = 1, nqmol(imty1)
cmov                  imind1 = iqmol(imty1,iqind1)
cmov                  iaty1 = iatype(imty1,imind1)
cmov                  iatom1 = iatmol(imol1,imind1)
cmov                  if (mask(iaty1)) then
cmov                     qepre(imol1,imol1) = qepre(imol1,imol1) + 
cmov     $                    2.0d0 * pepc * 
cmov     $                    (qeikrr(iatom1) * qeikrr(iatom1) +
cmov     $                    qeikri(iatom1) * qeikri(iatom1))
cmov                  endif
cmov                  do 254 iqind2 = 1, iqind1 - 1
cmov                     imind2 = iqmol(imty1,iqind2)
cmov                     iaty2 = iatype(imty1,imind2)
cmov                     iatom2 = iatmol(imol2,imind2)
cmov                     if (mask(iaty1) .or. mask(iaty2)) then
cmov                        qepre(imol1,imol1) = qepre(imol1,imol1) +
cmov     $                       4.0d0 * pepc *
cmov     $                       (qeikrr(iatom1) * qeikrr(iatom2) +
cmov     $                       qeikri(iatom1) * qeikri(iatom2))
cmov                     endif
cmov 254              continue
cmov                  do 256 imol2 = 1, imol1 - 1
cmov                     imty2 = molty(imol2)
cmov                     do 255 iqind2 = 1, nqmol(imty2)
cmov                        imind2 = iqmol(imty2,iqind2)
cmov                        iaty2 = iatype(imty2,imind2)
cmov                        iatom2 = iatmol(imol2,imind2)
cmov                        if (mask(iaty1) .or. mask(iaty2)) then
cmov                           qepre(imol1,imol2) = qepre(imol1,imol2) +
cmov     $                          2.0d0 * pepc * 
cmov     $                          (qeikrr(iatom1) * qeikrr(iatom2) +
cmov     $                          qeikri(iatom1) * qeikri(iatom2))
cmov                        endif
cmov 255                 continue
cmov 256              continue
cmov 257           continue
cmov 258        continue
cmov         endif
c***********************************************************************
c     calculate the force on each atom.  
c***********************************************************************
         do 265 imol = 1, nmol
            imty = molty(imol)
            do 260 iqind = 1, nqmol(imty)
               imind = iqmol(imty,iqind)
               iaty = iatype(imty,imind)
               iatom = iatmol(imol,imind)
c***********************************************************************
c     The force is split b/w contributions from light and heavy atoms:
c***********************************************************************
               fcpcl = eikrim(iatom,ivec) * scsubr - 
     $              eikrre(iatom,ivec) * scsubi
               fcpch = eikrim(iatom,ivec) * scsupr - 
     $              eikrre(iatom,ivec) * scsupi
c***********************************************************************
c     Which force subarray each piece contributes to depends on whether
c     this atom (iatom) is light or heavy:
c***********************************************************************
               if (masmsk(iaty,ilight)) then
                  fcpcm = fcpch
                  fcpch = 0.d0
               else
                  fcpcm = fcpcl
                  fcpcl = 0.d0
               endif
               fcpc = fcpcl + fcpcm + fcpch
               if (dsiter) then
                  fcx = rkvec(ivec,1) * fcpc
                  fcy = rkvec(ivec,2) * fcpc
                  fcz = rkvec(ivec,3) * fcpc
                  field(iatom,1) = field(iatom,1) + fcx
                  field(iatom,2) = field(iatom,2) + fcy
                  field(iatom,3) = field(iatom,3) + fcz
               endif
               fcpcx = q(iatom) * rkvec(ivec,1)
               fcpcy = q(iatom) * rkvec(ivec,2)
               fcpcz = q(iatom) * rkvec(ivec,3)
c***********************************************************************
c     Don't bother updating any pieces of the force array that we don't
c     care about or that haven't changed:
c***********************************************************************
               if (getlt .and. movlt) then
                  force(iatom,1,ifar,ilight) = 
     $                 force(iatom,1,ifar,ilight) + fcpcl * fcpcx
                  force(iatom,2,ifar,ilight) = 
     $                 force(iatom,2,ifar,ilight) + fcpcl * fcpcy
                  force(iatom,3,ifar,ilight) = 
     $                 force(iatom,3,ifar,ilight) + fcpcl * fcpcz
               endif
               if ((getlt .or. gethvy) .and. 
     $              (movlt .or. movhvy)) then
                  force(iatom,1,ifar,imixed) = 
     $                 force(iatom,1,ifar,imixed) + fcpcm * fcpcx
                  force(iatom,2,ifar,imixed) = 
     $                 force(iatom,2,ifar,imixed) + fcpcm * fcpcy
                  force(iatom,3,ifar,imixed) = 
     $                 force(iatom,3,ifar,imixed) + fcpcm * fcpcz
               endif
               if (gethvy .and. movhvy) then
                  force(iatom,1,ifar,iheavy) = 
     $                 force(iatom,1,ifar,iheavy) + fcpch * fcpcx
                  force(iatom,2,ifar,iheavy) = 
     $                 force(iatom,2,ifar,iheavy) + fcpch * fcpcy
                  force(iatom,3,ifar,iheavy) = 
     $                 force(iatom,3,ifar,iheavy) + fcpch * fcpcz
               endif
  260       continue
 265     continue
c***********************************************************************
c     Don't pay attention to mass masks when evalutating the virial --
c     pressure should only be calculated with all of the forces on.
c***********************************************************************
c     The virial as calculated here contributes positively to the 
c     pressure (i.e., no more sign changes).  Hard to say rigorously
c     how each term contributes to the pressure, although after summing
c     over k-space the first bit is negative and the second bit should
c     average to zero.
c***********************************************************************
c     There are no spline-induced forces because the Ewald potential is
c     not splined.
c***********************************************************************
         if (ioiter .and. gete) then
            wvirk = wvirk + (0.5d0 - r1kfac(ivec)) * 
     $           (scacre * acreal + scacim * acimag)
            wvirkd = wvirkd + ackdre * scacim - ackdim * scacre
         endif
c***********************************************************************
c     calculate the contribution to the electronegativity:
c***********************************************************************
         do 275 imol = 1, nmol
            imty = molty(imol)
            do 270 iqind = 1, nfqmol(imty)
               imind = ifqmol(imty,iqind)
               iaty = iatype(imty,imind)
               iatom = iatmol(imol,imind)
c***********************************************************************
c     The charge force is split b/w contributions from light and heavy 
c     atoms:
c***********************************************************************
               chipcl = eikrre(iatom,ivec) * scsubr +
     $              eikrim(iatom,ivec) * scsubi
               chipch = eikrre(iatom,ivec) * scsupr +
     $              eikrim(iatom,ivec) * scsupi
c***********************************************************************
c     Which chi subarray each piece contributes to depends on whether
c     this atom (iatom) is light or heavy:
c***********************************************************************
               if (masmsk(iaty,ilight)) then
                  chipcm = chipch
                  chipch = 0.d0
               else
                  chipcm = chipcl
                  chipcl = 0.d0
               endif
c***********************************************************************
c     Don't bother updating any pieces of the force array that we don't
c     care about or that haven't changed:
c***********************************************************************
               if (getlt .and. movlt) then
                  chi(iatom,ifar,ilight) = chi(iatom,ifar,ilight) +
     $                 chipcl
               endif
               if ((getlt .or. gethvy) .and. 
     $              (movlt .or. movhvy)) then
                  chi(iatom,ifar,imixed) = chi(iatom,ifar,imixed) + 
     $                 chipcm
               endif
               if (gethvy .and. movhvy) then
                  chi(iatom,ifar,iheavy) = chi(iatom,ifar,iheavy) +
     $                 chipch
               endif
 270        continue
 275     continue
c***********************************************************************
c     calculate the charge Hessian for use in solving for the charges.
c     this should always be done with no mass masking:
c***********************************************************************
         if (qsiter) then
            do 295 imol1 = 1, nmol
               imty1 = molty(imol1)
               do 290 iqind1 = 1, nqmol(imty1)
                  imind1 = iqmol(imty1,iqind1)
                  iaty1 = iatype(imty1,imind1)
                  iatom1 = iatmol(imol1,imind1)
                  qmat(iatom1,iatom1) = qmat(iatom1,iatom1) +
     $                 pepc * 
     $                 (eikrre(iatom1,ivec) * eikrre(iatom1,ivec) +
     $                 eikrim(iatom1,ivec) * eikrim(iatom1,ivec))
                  do 278 iqind2 = 1, iqind1 - 1
                     imind2 = iqmol(imty1,iqind2)
                     iaty2 = iatype(imty1,imind2)
                     iatom2 = iatmol(imol1,imind2)
                     qmatrm = pepc * 
     $                    (eikrre(iatom1,ivec) * eikrre(iatom2,ivec) +
     $                    eikrim(iatom1,ivec) * eikrim(iatom2,ivec))
                     qmat(iatom1,iatom2) = qmat(iatom1,iatom2) +
     $                    qmatrm
                     qmat(iatom2,iatom1) = qmat(iatom2,iatom1) +
     $                    qmatrm
 278              continue
                  do 285 imol2 = 1, imol1 - 1
                     imty2 = molty(imol2)
                     do 280 iqind2 = 1, nqmol(imty2)
                        imind2 = iqmol(imty2,iqind2)
                        iaty2 = iatype(imty2,imind2)
                        iatom2 = iatmol(imol2,imind2)
                        qmatrm = pepc * 
     $                       (eikrre(iatom1,ivec) * 
     $                       eikrre(iatom2,ivec) +
     $                       eikrim(iatom1,ivec) * eikrim(iatom2,ivec))
                        qmat(iatom1,iatom2) = 
     $                       qmat(iatom1,iatom2) + qmatrm
                        qmat(iatom2,iatom1) = 
     $                       qmat(iatom2,iatom1) + qmatrm
 280                 continue
 285              continue
 290           continue
 295        continue
         endif
 300  continue
      end
c
c***********************************************************************
c     gtcos() calculates cos(n*r(i,j)) and sin(n*r(i,j) for all n from 
c     -nbig to nbig, and all x-, y-, and z-positions of all atoms.  This
c     is done using deMoivre's identity: 
c                (cos nx + i sin nx) = (cos x + i sin x)^n
c     It is called every timestep from inside ewlrcp.  These are used in
c     solving for exp(i*k.r(i))....sjs 9/10/93
c***********************************************************************
c     updated for light/heavy RESPA split, to only recalculate for atoms
c     that have moved....sjs 6/26/95
c***********************************************************************
c     minor correction: different nbig factors could speed up the inner
c     loop for non-cubic boxes.
c***********************************************************************
c
      subroutine gtcos()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      do 305 imol = 1, nmol
         imty = molty(imol)
         do 300 iqind = 1, nqmol(imty)
            imind = iqmol(imty,iqind)
            iaty = iatype(imty,imind)
            if ((masmsk(iaty,ilight) .and. movlt) .or.
     $           (masmsk(iaty,iheavy) .and. movhvy)) then
               iatom = iatmol(imol,imind)
               scalex = rntok(1) * fldpos(iatom,1)
               scaley = rntok(2) * fldpos(iatom,2)
               scalez = rntok(3) * fldpos(iatom,3)
               coskx(iatom,2,1) = cos(scalex)
               coskx(iatom,2,2) = cos(scaley)
               coskx(iatom,2,3) = cos(scalez)
               sinkx(iatom,2,1) = sin(scalex)
               sinkx(iatom,2,2) = sin(scaley)
               sinkx(iatom,2,3) = sin(scalez)
               do 200 n = 2, nbig
                  np1 = n + 1
                  coskx(iatom,np1,1) = 
     $                 coskx(iatom,n,1) * coskx(iatom,2,1) -
     $                 sinkx(iatom,n,1) * sinkx(iatom,2,1)
                  coskx(iatom,np1,2) = 
     $                 coskx(iatom,n,2) * coskx(iatom,2,2) -
     $                 sinkx(iatom,n,2) * sinkx(iatom,2,2)
                  coskx(iatom,np1,3) = 
     $                 coskx(iatom,n,3) * coskx(iatom,2,3) -
     $                 sinkx(iatom,n,3) * sinkx(iatom,2,3)
                  sinkx(iatom,np1,1) = 
     $                 coskx(iatom,n,1) * sinkx(iatom,2,1) +
     $                 sinkx(iatom,n,1) * coskx(iatom,2,1)
                  sinkx(iatom,np1,2) = 
     $                 coskx(iatom,n,2) * sinkx(iatom,2,2) +
     $                 sinkx(iatom,n,2) * coskx(iatom,2,2)
                  sinkx(iatom,np1,3) = 
     $                 coskx(iatom,n,3) * sinkx(iatom,2,3) +
     $                 sinkx(iatom,n,3) * coskx(iatom,2,3)
  200          continue
            endif
 300     continue
 305  continue
      return
      end

