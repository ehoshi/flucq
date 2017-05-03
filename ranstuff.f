c***********************************************************************
c     qran is a (quick) normal, linear congruential random number  
c     generator with "well-chosen" constants taken from Numerical
c     Recipes.  Call it with a positive integer seed, it'll change that
c     seed and return a number in [0,1]...sjs 3/18/92
c***********************************************************************
c     Note:  Numerical Recipes doesn't SAVE the variables, which it
c     should.
c***********************************************************************
      function qran(iseed)
c
      include 'implic'
c
      parameter (imod = 134456, imult = 3877, iadd = 29573)
c
      save
c
      real*8 qran
c
      iseed = mod(imult * iseed + iadd, imod)
      qran = dble(iseed) / dble(imod)
      return
      end
c***********************************************************************
c     gran is a (good) more reliable random number generator based on 3
c     linear congruential cycles, with shuffling.  It was taken straight
c     from Numerical Recipes....sjs 3/18/92
c***********************************************************************
c     Note:  Numerical Recipes doesn't SAVE the variables, which it 
c     should.
c***********************************************************************
      function gran(iseed)
c
      include 'implic'
c
      logical init
c
      parameter (isize = 97)
      parameter (imod1 = 259200, imult1 = 7141, iadd1 = 54773, 
     $           div1 = 1.d0 / imod1)
      parameter (imod2 = 134456, imult2 = 8121, iadd2 = 28411, 
     $           div2 = 1.d0 / imod2)
      parameter (imod3 = 243000, imult3 = 4561, iadd3 = 51349)
c
      save
c
      dimension store(isize)
c
      data init/.true./
c
      if (iseed .lt. 0 .or. init) then
         init = .false.
         iseed1 = mod(iadd1 - iseed, imod1)
         iseed1 = mod(imult1 * iseed1 + iadd1, imod1)
         iseed2 = mod(iseed1, imod2)
         iseed1 = mod(imult1 * iseed1 + iadd1, imod1)
         iseed3 = mod(iseed1, imod3)
         do 11 islot = 1, isize
            iseed1 = mod(imult1 * iseed1 + iadd1, imod1)
            iseed2 = mod(imult2 * iseed2 + iadd2, imod2)
            store(islot) = (dble(iseed1) + dble(iseed2) * div2) * div1
 11      continue
c         iseed = 1
      endif
      iseed1 = mod(imult1 * iseed1 + iadd1, imod1)
      iseed2 = mod(imult2 * iseed2 + iadd2, imod2)
      iseed3 = mod(imult3 * iseed3 + iadd3, imod3)
      islot = 1 + (isize * iseed3) / imod3
      gran = store(islot)
      store(islot) = (dble(iseed1) + dble(iseed2) * div2) * div1
      return
      end




