c***********************************************************************
c     Function J calculates the r-dependent part of the Coulomb energy
c     for two partial charges.  The value returned is in units of Jmd.
c     For distances less than dlimr, the arrays calculated in Jinit() 
c     contain the shielded potential.  For longer distances, 1/r does 
c     just fine...sjs 5/1/92
c     If there are cutoff restrictions, they should be checked before
c     calling J()...sjs 6/3/92
c     Added the cutJ cutoff....sjs 7/2/93
c     Now decides whether to used 1/r, coulomb integral, or scaled
c     Coulomb integral based on the Jtype argument...sjs 7/3/93
c***********************************************************************
      real*8 function J(imty1, imty2, imind1, imind2, r, Jtype)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
c***********************************************************************
c     Case 1: J = 1/r
c***********************************************************************
      if (Jtype .eq. i1byr .or. r .gt. dlimr) then
         J = epsinv / r
c***********************************************************************
c     Case 2: for intramolecular interactions in rigid molecules, J is 
c     precalculated and stored in trmj
c***********************************************************************
      else if (Jtype .eq. Jisol) then
         if (imty1 .ne. imty2) then
            write(iuout, *) 'J: intermolecular bond distances aren''t ',
     $           'constant - can''t use Jisol'
            stop
         endif
         J = trmj(imty1, imind1, imind2)
c***********************************************************************
c     Case 3: general J(r), stored in Jdat array, evaluate using canned
c     cubic spline interpolation routine.  remember that
c     atom types must be in the right order when referencing Jdat:
c***********************************************************************
      else if (Jtype .eq. Jreg) then
c$$$         write(iuout, *) 'J: Jdat is not set'
c$$$         stop
         iaty1 = iatype(imty1, imind1)
         iaty2 = iatype(imty2, imind2)
         if (iaty1 .lt. iaty2) then
            iiaty1 = iaty2
            iiaty2 = iaty1
         else
            iiaty1 = iaty1
            iiaty2 = iaty2
         endif
c$$$         bindis = dlimr / (numbin - 1)
c$$$         ibin = int(r / bindis)
c$$$         J = Jdat(ibin + 1,i1,i2) + (r - ibin * bindis) / bindis
c$$$     $        * (Jdat(ibin + 2,i1,i2) - Jdat(ibin + 1,i1,i2))
c$$$         J = Jdat(ibin + 1,iiaty1,iiaty2) + (r - ibin * bindis) / 
c$$$     $        bindis * (Jdat(ibin + 2,iiaty1,iiaty2) - 
c$$$     $        Jdat(ibin + 1,iiaty1,iiaty2))
         init = 4
         npt = 1
         call dcsint(rvec, Jdat(1,iiaty1,iiaty2), 
     $        splcj(1,1,iiaty1,iiaty2), numbin, init, r, J, npt)
      else
         write(iuout, *) 'J: unknown J type'
         stop
      endif
      return
      end
c
c***********************************************************************
c     Function dJdr calculates the derivative of the J function
c     calculated in J().  For distances less than dlimr, the derivative
c     is calculated numerically from the Jdat arrays.  For longer
c     distances, it's simply -1/r^2...sjs 5/1/92
c     There's no check comparing r to rcut here - that should be done
c     before calling this routine...sjs 6/2/92
c     Accepts a Jtype argument to select function type, as in J(), and
c     knows about cutJ....sjs 7/3/93
c***********************************************************************
c     Minor correction:  make the derivative continuous by blending,
c     like in J().
c     minor correction: make this faster by accepting an rinv arg.
c***********************************************************************
      real*8 function dJdr(imty1, imty2, imind1, imind2, r, Jtype)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
c***********************************************************************
c     Case 1: dJ/dr = -1/r^2
c***********************************************************************
      if (Jtype .eq. i1byr .or. r .gt. dlimr) then
         dJdr = -epsinv / (r * r)
c***********************************************************************
c     Case 2: for intramolecular interactions in rigid molecules, dJ/dr
c     is precalculated and stored in trmdj
c***********************************************************************
      else if (Jtype .eq. Jisol) then
         if (imty1 .ne. imty2) then
            write(iuout, *) 'dJdr: intermolecular distances aren''t ',
     $           'constant - can''t use Jisol'
            stop
         endif
         J = trmdj(imty1,imind1,imind2)
c***********************************************************************
c     Case 3: general dJ/dr, stored in dJdat array.  evaluate using the
c     canned cubic spline interpolation routine.  remember that
c     atom types must be in the right order when referencing Jdat:
c***********************************************************************
      else if (Jtype .eq. Jreg) then
c$$$         write(iuout, *) 'dJdr: Jdat is not set and dJdat doens''t ',
c$$$     $        'exist yet'
c$$$         stop
         iaty1 = iatype(imty1,imind1)
         iaty2 = iatype(imty2,imind2)
         if (iaty1 .lt. iaty2) then
            iiaty1 = iaty2
            iiaty2 = iaty1
         else
            iiaty1 = iaty1
            iiaty2 = iaty2
         endif
c$$$         bindis = dlimr / (numbin - 1)
c$$$         ibin = int(r / bindis)
c$$$         dJdr = dJdat(ibin + 1,iiaty1,iiaty2) + (r - ibin * bindis) / 
c$$$     $        bindis * (dJdat(ibin + 2,iiaty1,iiaty2) - 
c$$$     $        dJdat(ibin + 1,iiaty1,iiaty2))
c$$$         dJdr = (Jdat(ibin + 2, i1, i2) - Jdat(ibin + 1, i1, i2)) /
c$$$     $        bindis
         init = 4
         npt = 1
         call dcsint(rvec, dJdat(1,iiaty1,iiaty2), 
     $        splcdj(1,1,iiaty1,iiaty2), numbin, init, r, dJdr, npt)
      else
         write(iuout, *) 'dJdr: unknown J type'
         stop
      endif
      return
      end
c
c***********************************************************************
c     quaddj calculates dJ/dr using finite difference, with J(r)
c     calculated from quadj (ie. by quadrature)....sjs 6/6/95
c***********************************************************************
c     This could just as easily be done by quadrature, but it's not.
c***********************************************************************
c
      function quaddj(iatyp1, iatyp2, r)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
      dr = 0.00001
      if (r .eq. 0.d0) then
         quaddj = 0.d0
      else
         quaddj = (quadj(iatyp1, iatyp2, r + dr) -
     $        quadj(iatyp1, iatyp2, r - dr)) / 2.d0 / dr
      endif
      return
      end
c
c***********************************************************************
c     quadj calculates the J(r) function using quadrature.
c     ...sjs 6/2/95
c***********************************************************************
c
      function quadj(iatyp1, iatyp2, r)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
c***********************************************************************
c     Figure out what kinds of atoms these are, and what size s orbitals
c     they need:
c***********************************************************************
      id1 = ideaty(iatyp1)
      id2 = ideaty(iatyp2)
      if (id1 .le. 2) then
         ishel1 = 1
      else if (id1 .le. 10) then
         ishel1 = 2
      else if (id1 .le. 18) then
         ishel1 = 3
      else if (id1 .eq. 35) then
c        (pretend Br is in row 3)
         ishel1 = 3
      else
         write(iuout, *) 'quadj: unknown atom type'
         stop
      endif
      if (id2 .le. 2) then
         ishel2 = 1
      else if (id2 .le. 10) then
         ishel2 = 2
      else if (id2 .le. 18) then
         ishel2 = 3
      else if (id2 .eq. 35) then
c        (pretend Br is in row 3)
         ishel2 = 3
      else
         write(iuout, *) 'quadj: unknown atom type'
         stop
      endif
c***********************************************************************
c     Make sure the smaller one comes first:
c***********************************************************************
      if (ishel2 .lt. ishel1) then
         iishl1 = ishel2
         iishl2 = ishel1
         iiaty1 = iatyp2
         iiaty2 = iatyp1
      else
         iishl1 = ishel1
         iishl2 = ishel2
         iiaty1 = iatyp1
         iiaty2 = iatyp2
      endif
      zeta1 = trmq(iiaty1,2)
      zeta2 = trmq(iiaty2,2)
c***********************************************************************
c     loop over k values from 0 to infinity, and integrate the 1-d
c     integral to get J(r):
c***********************************************************************
      dk = rjkmax / njk
      rk = 0.d0
      quadj = 0.5d0 * rjkern(iishl1, iishl2, zeta1, zeta2, r, rk)
      do 200 ik = 1, njk - 1
         rk = dk * ik
         quadj = quadj + rjkern(iishl1, iishl2, zeta1, zeta2, r, rk)
  200 continue
      rk = rjkmax
      quadj = quadj + 0.5d0 * 
     $     rjkern(iishl1, iishl2, zeta1, zeta2, r, rk)
      quadj = quadj * 2.d0 / pi * dk * epsinv
      return
      end
c
c***********************************************************************
c     rjkern calculates the kernel of the integral needed to get a J(r)
c     by quadrature.  The functions used depend on the principal quantum
c     number of the Slater orbitals being used....sjs 6/2/95
c***********************************************************************
c
      function rjkern(in1, in2, zeta1, zeta2, r, rk)
c
      include 'implic'
      include 'genpar'
      include 'qpar'
      include 'commons'
c
      temp = r * rk
      if (temp .eq. 0.d0) then
         rjkern = 1.d0
      else
         rjkern = sin(temp) / temp
      endif
      temp = rk / zeta1
      rf = 1.d0 / (1.d0 + (0.5d0 * temp) ** 2) ** 2
      rfac = temp ** 2 * sqrt(rf)
      if (in1 .eq. 1) then
         rjkern = rjkern * rf
      else if (in1 .eq. 2) then
         rjkern = rjkern * rf * (1 + rfac * (-3.d0 / 4.d0 + 
     $        1.d0 / 8.d0 * rfac))
      else if (in1 .eq. 3) then
         rjkern = rjkern * rf * (1 + rfac * (-11.d0 / 6.d0 +
     $        rfac * (17.d0 / 16.d0 + rfac * (-1.d0 / 4.d0 +
     $        1.d0 / 48.d0 * rfac))))
      else
         write(iuout, *) 'rjkern: illegal Slater n of ', in1
         stop
      endif
      temp = rk / zeta2
      rf = 1.d0 / (1.d0 + (0.5d0 * temp) ** 2) ** 2
      rfac = temp ** 2 * sqrt(rf)
      if (in2 .eq. 1) then
         rjkern = rjkern * rf
      else if (in2 .eq. 2) then
         rjkern = rjkern * rf * (1 + rfac * (-3.d0 / 4.d0 +
     $        1.d0 / 8.d0 * rfac))
      else if (in2 .eq. 3) then
         rjkern = rjkern * rf * (1 + rfac * (-11.d0 / 6.d0 +
     $        rfac * (17.d0 / 16.d0 + rfac * (-1.d0 / 4.d0 +
     $        1.d0 / 48.d0 * rfac))))
      else
         write(iuout, *) 'rjkern: illegal Slater n of ', in2
         stop
      endif
      return
      end
c
c***********************************************************************
c     analj calculates the J(r) function between atoms of type iatpy1
c     and iatyp2, separated by a distance r.  This is done analytically,
c     using Joel
c     Bader's formulae.  These are usually stored in a lookup table, or
c     else only certain values are kept.  The returned value has units
c     of Jmd/e^2.  Care must be taken to call this function with
c     known ids (currently 1 and 8) and r >= 0, otherwise it will gladly
c     return garbage....sjs 7/27/94
c***********************************************************************
c     Still need to add 1s-1s and 2s-2s interactions b/w different
c     zetas.
c***********************************************************************
c
      function analj(iatyp1, iatyp2, r)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
c***********************************************************************
c     Figure out what kinds of atoms these are, and what size s orbital
c     they need:
c***********************************************************************
      id1 = ideaty(iatyp1)
      if (id1 .eq. 1) then
         ishel1 = 1
      else if (id1 .eq. 8) then
         ishel1 = 2
      else if (id1 .eq. 17) then
         ishel1 = 3
      else
         write(iuout, *) 'analj: unknown atom type'
         stop
      endif
      id2 = ideaty(iatyp2)
      if (id2 .eq. 1) then
         ishel2 = 1
      else if (id2 .eq. 8) then
         ishel2 = 2
      else if (id2 .eq. 17) then
         ishel2 = 3
      else
         write(iuout, *) 'analj: unknown atom type'
         stop
      endif
c***********************************************************************
c     Make sure the smaller one comes first:
c***********************************************************************
      if (ishel2 .lt. ishel1) then
         iishl1 = ishel2
         iishl2 = ishel1
         iiaty1 = iatyp2
         iiaty2 = iatyp1
      else
         iishl1 = ishel1
         iishl2 = ishel2
         iiaty1 = iatyp1
         iiaty2 = iatyp2
      endif
c***********************************************************************
c     Solve for J(r):
c***********************************************************************
      zeta1 = trmq(iiaty1,2)
      zeta2 = trmq(iiaty2,2)
      zetar1 = zeta1 * r
      zetar2 = zeta2 * r
c***********************************************************************
c     1s-2s:
c***********************************************************************
      if (iishl1 .eq. 1 .and. iishl2 .eq. 2) then
         if (r .eq. 0.0d0) then
            zeta12 = zeta1 * zeta1
            zeta13 = zeta12 * zeta1
            zeta14 = zeta13 * zeta1
            tnum = zeta1 * zeta2 * (zeta14 + zeta2 * 
     $           (5. * zeta13 + zeta2 * (10. * zeta12 + zeta2 * 
     $           (10. * zeta1 + 2. * zeta2))))
            tdenom = 2. * (zeta1 + zeta2) ** 5
         else
            zetr12 = zetar1 * zetar1
            zetr14 = zetr12 * zetr12
            zetr22 = zetar2 * zetar2
            zetr24 = zetr22 * zetr22
            zetr26 = zetr24 * zetr22
            texpx = exp(2. * zetar1)
            texpy = exp(2. * zetar2)
            texpxy = texpx * texpy * 6. * (zetr22 - zetr12) ** 5
            ta = texpx * zetr14 * (zetr26 * (-84. + zetar2 * 
     $           (-63. + zetar2 * (-18. - 2. * zetar2))) + zetr12 * 
     $           (zetr24 * (60. + zetar2 * (99. + zetar2 * 
     $           (42. + 6. * zetar2))) + zetr12 * (zetr22 * 
     $           (-30. + zetar2 * (-45. + zetar2 * 
     $           (-30. - 6. * zetar2))) + zetr12  * (6. + zetar2 * 
     $           (9. + zetar2 * (6. + 2. * zetar2))))))
            tb = texpy * zetr26 * (zetr14 * (24. + 6. * zetar1) + 
     $           zetr22 * (30. * zetr12 + zetr22 * 
     $           (-6. - 6. * zetar1)))
            tc = texpxy
            tnum = texpx * zetr14 * (zetr26 * (-84. + zetar2 * 
     $           (-63. + zetar2 * (-18. - 2. * zetar2))) + zetr12 * 
     $           (zetr24 * (60. + zetar2 * (99. + zetar2 * 
     $           (42. + 6. * zetar2))) + zetr12 * (zetr22 * 
     $           (-30. + zetar2 * (-45. + zetar2 * 
     $           (-30. - 6. * zetar2))) + zetr12  * (6. + zetar2 * 
     $           (9. + zetar2 * (6. + 2. * zetar2)))))) + 
     $           texpy * zetr26 * (zetr14 * (24. + 6. * zetar1) + 
     $           zetr22 * (30. * zetr12 + zetr22 * 
     $           (-6. - 6. * zetar1))) + texpxy
            tdenom = texpxy * r
         endif
c***********************************************************************
c     1s-1s:
c***********************************************************************
      else if (iishl1 .eq. 1 .and. iishl2 .eq. 1) then
         if (r .eq. 0.0d0) then
            tnum = 5. * zeta1
            tdenom = 8.
         else
            texpfc = 24. * exp(2. * zetar1)
            tnum = -24. + texpfc + 
     $           zetar1 * (-33. + zetar1 * (-18. - 4. * zetar1))
            tdenom = texpfc * r
         endif
c***********************************************************************
c     2s-2s:
c***********************************************************************
      else if (iishl1 .eq. 2 .and. iishl2 .eq. 2) then
         if (r .eq. 0.0d0) then
            tnum = 93. * zeta1
            tdenom = 256.
         else
            texpfc = 80640. * exp(2. * zetar1)
            tnum = -80640. + texpfc + zetar1 * (-131985. + zetar1 * 
     $           (-102690. + zetar1 * (-49980. + zetar1 * (-16800. + 
     $           zetar1 * (-4032. + zetar1 * (-672. - 64. * zetar1))))))
            tdenom = texpfc * r
         endif
      endif
      analj = epsinv * tnum / tdenom
      return
      end
c
c***********************************************************************
c     analdj calculates dJ(r)/dr for J(r) between atoms imid1 and imid2
c     in molecules of types imty1 and imty2.  Their separation is taken
c     from the rmat array.  It does so analytically, using Joel
c     Bader's formulae.  These are usually stored in a lookup table, or
c     else only certain values are kept.  The returned value has units
c     of Jmd/e^2/A.  Care must be taken to call this function with
c     known ids (currently 1 and 8) and r >= 0, otherwise it will gladly
c     return garbage....sjs 7/27/94
c***********************************************************************
c     Still need to add 1s-1s and 2s-2s interactions b/w different
c     zetas.
c***********************************************************************
c
      function analdj(iatyp1, iatyp2, r)
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
c***********************************************************************
c     Figure out what kinds of atoms these are, and what size s orbital
c     they need:
c***********************************************************************
      id1 = ideaty(iatyp1)
      if (id1 .eq. 1) then
         ishel1 = 1
      else if (id1 .eq. 8) then
         ishel1 = 2
      endif
      id2 = ideaty(iatyp2)
      if (id2 .eq. 1) then
         ishel2 = 1
      else if (id2 .eq. 8) then
         ishel2 = 2
      endif
c***********************************************************************
c     Make sure the smaller one comes first:
c***********************************************************************
      if (ishel2 .lt. ishel1) then
         iishl1 = ishel2
         iishl2 = ishel1
         iiaty1 = iatyp2
         iiaty2 = iatyp1
      else
         iishl1 = ishel1
         iishl2 = ishel2
         iiaty1 = iatyp1
         iiaty2 = iatyp2
      endif
c***********************************************************************
c     Solve for J(r):
c***********************************************************************
      zeta1 = trmq(iiaty1,2)
      zeta2 = trmq(iiaty2,2)
      zetar1 = zeta1 * r
      zetar2 = zeta2 * r
c***********************************************************************
c     1s-2s:
c***********************************************************************
      if (iishl1 .eq. 1 .and. iishl2 .eq. 2) then
         write(iuout, *) 'analdj: 1s-2s dJ/dr is broken'
         stop
         if (r .eq. 0.0d0) then
            tnum = 0.0d0
            tdenom = 1.0d0
            tdnum = 0.0d0
            tddnom = 0.0d0
         else
            zetr12 = zetar1 * zetar1
            zetr13 = zetr12 * zetar1
            zetr14 = zetr12 * zetr12
            zetr22 = zetar2 * zetar2
            zetr23 = zetr22 * zetar2
            zetr24 = zetr22 * zetr22
            zetr25 = zetr24 * zetar2
            zetr26 = zetr24 * zetr22
            texpx = exp(2. * zetar1)
            texpy = exp(2. * zetar2)
            texpxy = texpx * texpy
            ta = zetr14 * (zetr26 * (-84. + zetar2 * 
     $           (-63. + zetar2 * (-18. - 2. * zetar2))) + zetr12 * 
     $           (zetr24 * (60. + zetar2 * (99. + zetar2 * 
     $           (42. + 6. * zetar2))) + zetr12 * (zetr22 * 
     $           (-30. + zetar2 * (-45. + zetar2 * 
     $           (-30. - 6. * zetar2))) + zetr12  * (6. + zetar2 * 
     $           (9. + zetar2 * (6. + 2. * zetar2))))))
            tb = zetr26 * (zetr14 * (24. + 6. * zetar1) + 
     $           zetr22 * (30. * zetr12 + zetr22 * 
     $           (-6. - 6. * zetar1)))
            tc = 6.d0 * (zetr22 - zetr12) ** 5
            tnum = texpx * ta + texpy * tb + tc * texpxy
            tnum = texpx * zetr14 * (zetr26 * (-84. + zetar2 * 
     $           (-63. + zetar2 * (-18. - 2. * zetar2))) + zetr12 * 
     $           (zetr24 * (60. + zetar2 * (99. + zetar2 * 
     $           (42. + 6. * zetar2))) + zetr12 * (zetr22 * 
     $           (-30. + zetar2 * (-45. + zetar2 * 
     $           (-30. - 6. * zetar2))) + zetr12  * (6. + zetar2 * 
     $           (9. + zetar2 * (6. + 2. * zetar2)))))) + 
     $           texpy * zetr26 * (zetr14 * (24. + 6. * zetar1) + 
     $           zetr22 * (30. * zetr12 + zetr22 * 
     $           (-6. - 6. * zetar1))) + tc * texpxy
            tdnum = texpx * (zetr13 * (zetr26 * zeta1 * (-336.d0 + 
     $           zetar2 * (-252.d0 + zetar2 * (-72.d0 - 8.d0 *
     $           zetar2))) + zetar1 * (zetr25 * zeta2 * (-504.d0 +
     $           zetar2 * (-441.d0 + zetar2 * (-144.d0 - 18.d0 * 
     $           zetar2))) + zetar1 * (zetr24 * zeta1 * (360.d0 + 
     $           zetar2 * (594.d0 + zetar2 * (252.d0 + 36.d0 * 
     $           zetar2))) + zetar1 * (zetr23 * zeta2 * (240.d0 + 
     $           zetar2 * (495.d0 + zetar2 * (252.d0 + 42.d0 * 
     $           zetar2))) + zetar1 * (zetr22 * zeta1 * (-240.d0 +
     $           zetar2 * (-360.d0 + zetar2 * (-240.d0 - 48.d0 *
     $           zetar2))) + zetar1 * (zetar2 * zeta2 * (-60.d0 + 
     $           zetar2 * (-135.d0 + zetar2 * (-120.d0 - 30.d0 * 
     $           zetar2))) + zetar1 * (zeta1 * (60.d0 + zetar2 * 
     $           (90.d0 + zetar2 * (60.d0 + 20.d0 * zetar2))) + 
     $           zetar1 * (zeta2 * (9.d0 + zetar2 * (12.d0 + 6.d0 *
     $           zetar2)))))))))) + 2.d0 * zeta1 * ta) + 
     $           texpy * (zeta1 * zetr26 * (zetr13 * (96.d0 + 30.d0 *
     $           zetar1) + zetr22 * (60.d0 * zetar1 - 6.d0 * zetr22)) +
     $           zeta2 * zetr25 * (zetr14 * (144.d0 + 36.d0 * zetar1) +
     $           zetr22 * (240.d0 * zetr12 + zetr22 * (-60.d0 - 60.d0 *
     $           zetar1))) + 2.d0 * zeta2 * tb) + texpxy * ((2.d0 *
     $           zeta1 + 2.d0 * zeta2) * tc + 30.d0 * (2.d0 * zeta1 *
     $           zetar1 - 2.d0 * zeta2 * zetar2) * 
     $           (zetr22 - zetr12) ** 4)
            tdenom = tc * texpxy * r
            tddnom = texpxy * tc * (2.d0 * (zetar1 + zetar2) + 11)
         endif
c***********************************************************************
c     1s-1s:
c***********************************************************************
      else if (iishl1 .eq. 1 .and. iishl2 .eq. 1) then
         if (r .eq. 0.0d0) then
            tnum = 0.0d0
            tdenom = 1.0d0
            tdnum = 0.0d0
            tddnom = 0.0d0
         else
            texpfc = 24. * exp(2. * zetar1)
            tnum = -24. + texpfc + 
     $           zetar1 * (-33. + zetar1 * (-18. - 4. * zetar1))
            tdenom = texpfc * r
            tdnum = zeta1 * (2.0d0 * texpfc - 33.0 + zetar1 * 
     $           (-36.0d0 - 12.0d0 * zetar1))
            tddnom = texpfc * (2.0d0 * zetar1 + 1.0d0)
         endif
c***********************************************************************
c     2s-2s:
c***********************************************************************
      else if (iishl1 .eq. 2 .and. iishl2 .eq. 2) then
         if (r .eq. 0.0d0) then
            tnum = 0.0d0
            tdenom = 1.0d0
            tdnum = 0.0d0
            tddnom = 0.0d0
         else
            texpfc = 80640. * exp(2. * zetar1)
            tnum = -80640. + texpfc + zetar1 * (-131985. + zetar1 * 
     $           (-102690. + zetar1 * (-49980. + zetar1 * (-16800. + 
     $           zetar1 * (-4032. + zetar1 * (-672. - 64. * zetar1))))))
            tdenom = texpfc * r
            tdnum = zeta2 * (2.0d0 * texpfc - 131985.d0 + zetar1 *
     $           (-205380.d0 + zetar1 * (-149940.d0 + zetar1 *
     $           (-67200.d0 + zetar1 * (-20160.d0 + zetar1 *
     $           (-4032.d0 - 448.d0 * zetar1))))))
            tddnom = texpfc * (2.d0 * zetar1 + 1)
         endif
      endif
      analdj = epsinv * (tdnum * tdenom - tnum * tddnom) / tdenom ** 2
      return
      end
c
c***********************************************************************
c     Jinit has been completely rewritten.  It still initializes the
c     Jdat stuff for all necessary pairs of atom types.  But it no
c     longer checks for a pre-existing data file, it always recalculates
c     the values it needs.  And instead of using Friesner's code for 
c     doing the contracted Gaussian stuff, it uses Joel Bader's 
c     equations to do the J(r) calc.  It is only called if either interJ
c     or intraJ is Jreg.  Jdat(ibin,i1,i2) must still be called with 
c     i1 > i2.  A factor of 1/(4*pi*epsilon) is multiplied into the Jdat
c     array....sjs 7/27/94
c***********************************************************************
c     another rewrite.  Instead of using Joel's analytical expressions
c     for J(r), they are done by quadrature.  This is because I only 
c     have Joel's expressions for 1s/1s, 1s/2s, and 2s/2s, and Steve
c     Rick's 1-d integral representation is more general, a little
c     easier to extend, and is worked out for combinations with 3s.
c     All of the J(r) are now done by quadrature, even when for the
c     cases which I know analytically, for consistency....sjs 6/6/95
c***********************************************************************
c     
      subroutine Jinit()
c
      include 'implic'
      include 'qpar'
      include 'genpar'
      include 'commons'
c
      character name*6
c
      write(iuout, *) 'Initializing J(r) data...'
c***********************************************************************
c     Set up the array of r points at which to evaluate J(r) and dJ/dr:
c***********************************************************************
      do 190 ibin = 1, numbin
         rvec(ibin) = (ibin - 1) * dlimr / (numbin - 1)
  190 continue
c***********************************************************************
c     loop over pairs of atom types
c***********************************************************************
      do 220 iatyp1 = 1, naty
         if (.not. isqat(iatyp1)) then
            go to 220
         endif
         id1 = ideaty(iatyp1)
         do 210 iatyp2 = iatyp1, 1, -1
            if (.not. isqat(iatyp2)) then
               go to 210
            endif
            id2 = ideaty(iatyp2)
            write(name, '(''J'', 2i2.2, ''.'')') id1, id2
            write(filnam, '(4a)') datdir(1:index(datdir, ' ')-1), '/', 
     .           name, 
     $           suffix(1:index(suffix, ' ')-1)
            open(iutmp, file = filnam, err = 9990, iostat = ioval)
c$$$            dr = 0.00001
c***********************************************************************
c     Calculate J(r) and dJ/dr from a 1-d integral expression and
c     finite difference:
c***********************************************************************
            do 200 ibin = 1,numbin
               Jdat(ibin,iatyp1,iatyp2) = 
     $              quadj(iatyp1, iatyp2, rvec(ibin))
               dJdat(ibin,iatyp1,iatyp2) = 
     $              quaddj(iatyp1, iatyp2, rvec(ibin))
               rJtemp = quadj(iatyp1, iatyp2, rdist)
               if (rdist .eq. 0.d0) then
                  rJtdr = 0.d0
               else
                  rJtdr = (quadj(iatyp1, iatyp2, rdist + dr) - 
     $                 quadj(iatyp1, iatyp2, rdist - dr)) / 2.d0 / dr
               endif
               write(iutmp, *) rvec(ibin), Jdat(ibin,iatyp1,iatyp2),
     $              dJdat(ibin,iatyp1,iatyp2)
 200        continue
            close(iutmp)
c***********************************************************************
c     Now that we have J(r) on a grid of points, fit it with a cubic
c     spline so we can have a smoothly-varying function.  The derivative
c     at the endpoints is known analytically.
c***********************************************************************
            init = 3
            npt = 0
            splcj(1,1,iatyp1,iatyp2) = 0.d0
            splcj(2,1,iatyp1,iatyp2) = -epsinv / dlimr ** 2
            call dcsint(rvec, Jdat(1,iatyp1,iatyp2), 
     $           splcj(1,1,iatyp1,iatyp2), numbin, init, rjunk, rjunk2,
     $           npt)
c***********************************************************************
c     We also have dJ/dr on a grid of points, so fit it with a cubic
c     spline as well.  The derivative at the large-r endpoint is known,
c     but not at r=0 (unless iatyp1 = iatyp2), so finite difference is
c     used.
c***********************************************************************
            dr = 0.0001
            zero = 0.d0
            splcdj(1,1,iatyp1,iatyp2) = 2.d0 * 
     $           (quadj(iatyp1, iatyp2, dr) -
     $           quadj(iatyp1, iatyp2, zero)) / dr ** 2
            splcdj(2,1,iatyp1,iatyp2) = 2.d0 * epsinv / dlimr ** 3
            init = 3
            npt = 0
            call dcsint(rvec, dJdat(1,iatyp1,iatyp2),
     $           splcdj(1,1,iatyp1,iatyp2), numbin, init, rjunk, rjunk2,
     $           npt)
 210     continue
 220  continue
      write(iuout, *) '                      ...finished.'
      return
 9990 write(iuout, *) 'Jinit: error opening file ', filnam
      call ioerr(ioval)
      stop
      end
c
