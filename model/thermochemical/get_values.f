!      program test; call tester; end program test

      subroutine get_values_tester
      implicit none

      integer          :: i
      double precision :: value1 , value2 , value3
       character(len=81)     :: line

      open (unit=10,file='BURCAT.THR',status='old')
      do i = 1,9765
        read (10,'(a80)') line
      end do

      write (*,'(1x,a80)') line

      call get_values(line,value1,value2,value3)

      write (*,*) value1
      write (*,*) value2
      write (*,*) value3

      stop
      end subroutine get_values_tester


      subroutine get_values(line,value1,value2,value3)

      implicit none

      character(*)    :: line
      character(len=80)     :: tokens(80)
      logical          :: bSpace
      integer          :: nTokens , start , i
      double precision :: value1 , value2 , value3

      bSpace = .false.
      nTokens = 0
      start = 1

      do i = 1,len(line)
        if (line(i:i).eq.char(9) .or. line(i:i).eq.char(32)) then
           if (.not.bSpace) then
              bSpace = .true.
               nTokens = nTokens + 1
              tokens(nTokens) = line(start:i-1)
           end if
        else
           if (bSpace) start = i
           bSpace = .false.
        end if
      end do

      if (start .ne. len(line)) nTokens = nTokens - 1

      if (nTokens .ge. 4) then
         read (tokens(nTokens  ),'(f16.8)') value3
         read (tokens(nTokens-1),'(f16.8)',err=111) value2
         read (tokens(nTokens-2),'(f16.8)') value1
         return
 111     continue
         read (tokens(nTokens-3),'(f16.8)') value1
         read (tokens(nTokens-2),'(f16.8)') value2
       end if


       return
       end
