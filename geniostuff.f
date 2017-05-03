c***********************************************************************
c     irdump simply writes out a particular number of elements from a
c     real*4 (important!) array that it is passed....sjs 9/14/92
c***********************************************************************
c
      subroutine irdump(iunit, rbuff, num)
c
      include 'implic'
      include 'genpar'
c
      real*4 rbuff(num)
c
      write(iunit, err = 9990, iostat = ioval) rbuff
      return
 9990 write(iuout, *) 'irdump: error writing to unit ', iunit
      call ioerr(ioval)
      stop
      end
c
c***********************************************************************
c     iddump simply writes out a particular number of elements from a
c     real*8 (important!) array that it is passed....sjs 9/21/94
c***********************************************************************
c
      subroutine iddump(iunit, dbuff, num)
c
      include 'implic'
      include 'genpar'
c
      real*8 dbuff(num)
c
      write(iunit, err = 9990, iostat = ioval) dbuff
      return
 9990 write(iuout, *) 'iddump: error writing to unit ', iunit
      call ioerr(ioval)
      stop
      end
c
c***********************************************************************
c     iidump simply writes out a particular number of elecments from a
c     integer*4 (important!) array that it is passed....sjs 7/25/95
c***********************************************************************
c
      subroutine iidump(iunit, ibuff, num)
c
      include 'implic'
      include 'genpar'
c
      integer*4 ibuff(num)
c
      write(iunit, err = 9990, iostat = ioval) ibuff
      return
 9990 write(iuout, *) 'iidump: error writing to unit ', iunit
      call ioerr(ioval)
      stop
      end
c
c***********************************************************************
c     iurdmp does the reverse of irdump - it reads a particular number
c     of real*4 (important!) elements into the array that it is passed.
c     (Simply reading the array name causes a silly warning message on 
c     the IBMs.)  It will crash on error, but return a -1 on EOF.
c     ...sjs 9/21/94
c***********************************************************************
c     modified to not crash on error....sjs 2/8/17
c***********************************************************************
c
      integer*4 function iurdmp(iunit, rbuff, num)
c
      include 'implic'
      include 'genpar'
c
      real*4 rbuff(num)
c
c      read(iunit, err = 9990, iostat = ioval) rbuff
      read(iunit, iostat = ioval) rbuff
      iurdmp = ioval
      return
 9990 write(iuout, *) 'iurdmp: error reading from unit ', iunit
      call ioerr(ioval)
      stop
      end
c
c***********************************************************************
c     iuddmp does the reverse of iddump - it reads a particular number
c     of real*8 (important!) elements into the array that it is passed.
c     (Simply reading the array name causes a silly warning message on 
c     the IBMs.).  It will crash on a read error, and return -1 on EOF.
c     ...sjs 9/17/92
c***********************************************************************
c
      integer*4 function iuddmp(iunit, dbuff, num)
c
      include 'implic'
      include 'genpar'
c
      real*8 dbuff(num)
c
      read(iunit, err = 9990, iostat = ioval) dbuff
      iuddmp = ioval
      return
 9990 write(iuout, *) 'iuddmp: error reading from unit ', iunit
      call ioerr(ioval)
      stop
      end
c
c***********************************************************************
c     iuidmp does the reverse of iidump - it reads a particular number
c     of integer*4 (important!) elements into the array that it is
c     passed.  (Simply reading the array name causes a silly warning
c     message on the IBMs.)  It will crash on a read error, and return
c     -1 on EOF....sjs 7/25/95
c***********************************************************************
c
      integer*4 function iuidmp(iunit, ibuff, num)
c
      include 'implic'
      include 'genpar'
c     
      integer*4 ibuff(num)
c
      read(iunit, err = 9990, iostat = ioval) ibuff
      iuidmp = ioval
      return
 9990 write(iuout, *) 'iuidmp: error reading from unit ', iunit
      call ioerr(ioval)
      stop
      end
c
c***********************************************************************
c     ioerr is an IBM AIX-specific (and terminal) routine designed to
c     interpret file I/O error conditions, since most are coded with
c     ERR labels.  The input variable is the value returned in IOSTAT
c     in the failed I/O operation...sjs 4/14/93
c***********************************************************************
c
      subroutine ioerr(ioval)
c
      include 'implic'
c
      if (ioval .eq. -1) then
         write(iuout, *) 'IOSTAT = ', ioval, ', End of file'
      else if (ioval .eq. -2) then
         write(iuout, *) 'IOSTAT = ', ioval, ', End of internal file'
      else if (ioval .eq. 6) then
         write(iuout, *) 'IOSTAT = ', ioval, 
     $        ', File not found and STATUS=OLD specified in OPEN'
      else if (ioval .eq. 7) then
         write(iuout, *) 'IOSTAT = ', ioval, 
     $        ', Incorrect format of list-directed input in file'
      else if (ioval .eq. 8) then
         write(iuout, *) 'IOSTAT = ', ioval, ', Incorrect format of ',
     $        'list-directed input in internal file'
      else if (ioval .eq. 9) then
         write(iuout, *) 'IOSTAT = ', ioval, 
     $        ', Data item too long for the internal file'
      else if (ioval .eq. 10) then
         write(iuout, *) 
     $        'IOSTAT = ', ioval, ', Read error on direct file'
      else if (ioval .eq. 11) then
         write(iuout, *) 
     $        'IOSTAT = ', ioval, ', Write error on direct file'
      else if (ioval .eq. 12) then
         write(iuout, *) 'IOSTAT = ', ioval,
     $        ', Read error on sequential file'
      else if (ioval .eq. 13) then
         write(iuout, *) 'IOSTAT = ', ioval,
     $        ', Write error on sequential file'
      else if (ioval .eq. 14) then
         write(iuout, *) 'IOSTAT = ', ioval, ', Error opening a file'
      else if (ioval .eq. 15) then
         write(iuout, *) 'IOSTAT = ', ioval, 
     $        ', Permanent I/O error encountered on a file'
      else if (ioval .eq. 16) then
         write(iuout, *) 'IOSTAT = ', ioval,
     $        ', Record number invalid for direct I/O'
      else if (ioval .eq. 17) then
         write(iuout, *) 'IOSTAT = ', ioval,
     $        ', I/O statement not allowed for direct file'
      else if (ioval .eq. 18) then
         write(iuout, *) 'IOSTAT = ', ioval,
     $        ', Direct I/O statement on a file not open'
      else if (ioval .eq. 19) then
         write(iuout, *) 'IOSTAT = ', ioval,
     $        ', Unformatted I/O on formatted file'
      else if (ioval .eq. 20) then
         write(iuout, *) 'IOSTAT = ', ioval,
     $        ', Formatted I/O on unformatted file'
      else if (ioval .eq. 21) then
         write(iuout, *) 'IOSTAT = ', ioval, 
     $        ', Sequential I/O on direct file'
      else if (ioval .eq. 22) then
         write(iuout, *) 'IOSTAT = ', ioval, 
     $        ', Direct I/O on sequential file'
      else if (ioval .eq. 23) then
         write(iuout, *) 'IOSTAT = ', ioval, 
     $        ', File connected to another unit'
      else if (ioval .eq. 24) then
         write(iuout, *) 'IOSTAT = ', ioval,
     $        ', OPEN specifiers don''t match file attributes'
      else if (ioval .eq. 25) then
         write(iuout, *) 'IOSTAT = ', ioval, 
     $        ', RECL not given on OPEN for direct file'
      else if (ioval .eq. 26) then
         write(iuout, *) 'IOSTAT = ', ioval, 
     $        ', Negative record length in OPEN'
      else if (ioval .eq. 27) then
         write(iuout, *) 'IOSTAT = ', ioval, 
     $        ', OPEN ACCESS specifier invalid'
      else if (ioval .eq. 28) then
         write(iuout, *) 'IOSTAT = ', ioval, 
     $        ', OPEN FORMAT specifier invalid'
      else if (ioval .eq. 31) then
         write(iuout, *) 'IOSTAT = ', ioval,
     $        ', OPEN FILE specifier invalid'
      else if (ioval .eq. 35) then
         write(iuout, *) 'IOSTAT = ', ioval, ', Recursive I/O operation'
      else if (ioval .eq. 36) then
         write(iuout, *) 'IOSTAT = ', ioval, ', Invalid unit number'
      else if (ioval .eq. 38) then
         write(iuout, *) 'IOSTAT = ', ioval, ', REWIND error'
      else if (ioval .eq. 39) then
         write(iuout, *) 'IOSTAT = ', ioval, ', ENDFILE error'
      else if (ioval .eq. 40) then
         write(iuout, *) 'IOSTAT = ', ioval, ', BACKSPACE error'
      else if (ioval .eq. 84) then
         write(iuout, *) 'IOSTAT = ', ioval, 
     $        ', NAMELIST name not found in file'
      else if (ioval .eq. 85) then
         write(iuout, *) 'IOSTAT = ', ioval,
     $        ', NAMELIST name not found in internal file'
      else if (ioval .eq. 93) then
         write(iuout, *) 'IOSTAT = ', ioval,
     $        ', I/O statement not allowed on error unit'
      else
         write(iuout, *) 'error type unknown'
      endif
      stop
      end
c
