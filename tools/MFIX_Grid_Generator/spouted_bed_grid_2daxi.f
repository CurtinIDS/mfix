!
!******************************************************************
!
!  Simple MFIX Grid Generator for a Spouted Bed with Stretched Grids
!
!	A very simple program to generate grid spacings and sloping section boundary
!	information for MFIX.DAT
!
!       Please use this as a basis and perform any additional customizations as
!       needed and this module does not come with any guarantees
!
!	When using stretched grids only FOUP works in MFIX
!
!  Author: Sreekanth Pannala  (May-16-06)
!
!******************************************************************
!
!

	program gridgen
	implicit real*8 (a-h,o-z)

!	Local variables - statically allocated but need to made allocatable

	dimension y(0:500), x(0:500), dy(0:500), dx(0:500)
	dimension dx1(10), xlen(10), i1(10), i2(10), irev(10), &
        xbegin(10), xend(10)
	dimension dy1(10), ylen(10), j1(10), j2(10), jrev(10), &
        ybegin(10), yend(10)
        integer one

!	Opening files for writing out dx, dy and sloping section information
!       The dxdat, dydat and slope files info can be directly imported to a
!       valid mfix.dat file to simulate the geometry specified in this file

	open(unit=15,file='dydat.data')
	open(unit=16,file='dxdat.data')
	open(unit=17,file='slope.data')

        pi = 4.0d0*atan(1.0d0)

        angle = 30.0d0 ! 1/2 the included angle of the sloping wall
        slope = 1.0d0/(tan(pi*angle/180.0d0))
        dxmin = 0.475d0/2.0d0 ! minimum dx you would like at the finest level


!	Define patches in x direction with various geometric information
!       irev - stretching direction (1- forward and 2 - reverse)

!       Inner pipe - X-coordinates

	nxelem = 2
        i1_1 = 0
        i2_1 = 2
	dx1(1) = dxmin
	xbegin(1) = 0.0d0
	xend(1) = 0.475d0
	i1(1) = i1_1
	i2(1) = i2_1
	irev(1) = 1

	dx1(2) = dxmin
        i1_2 = 2
        i2_2 = 4
	xbegin(2) = 0.475d0
	xend(2) = 0.95d0
	i1(2) = i1_2
	i2(2) = i2_2
	irev(2) = 1

!	Call geometric to return actual geometric information using geomtric stretching

	do nx = 1,nxelem

	   call geometric(dx1(nx), str, i1(nx),i2(nx), xbegin(nx), xend(nx), irev(nx), x, dx)

	   do i = i1(nx), i2(nx)
	      write(*,*) i,x(i)
	   enddo

	enddo

!	Define patches in y direction with various geometric information

        nyelem = 1
        j1_1 = 0
        j2_1 = 20
        dy1(1) = slope*dxmin ! minimum dy you would like at the finest level
        ybegin(1) = 0.0d0
        yend(1) = slope*(7.6d0-0.95d0)
        j1(1) = j1_1
        j2(1) = j2_1
        jrev(1) = 1


!	Call geometric to return actual geometric information using geometric stretching

        do ny = 1,nyelem

	   call geometric(dy1(ny), str, j1(ny),j2(ny), ybegin(ny), yend(ny), jrev(ny), y, dy)

	   do j = j1(ny), j2(ny)
	      write(*,*) 'j1', j,y(j)
	   enddo

        enddo

!	Define more patches in y direction with various geomtric information - can be
!       generalized and club with above

        j1_2 = j2_1
        j2_2 = 60
	nyelem = 1
        ybegin(1) = slope*(7.6d0-0.95d0)
        yend(1) = 90.0d0
        j1(1) = j1_2
        j2(1) = j2_2
        dy1(1) =  (y(j1(1)) - y(j1(1)-1))*str
        jrev(1) = 1

        do ny = 1,nyelem

	   call geometric(dy1(ny), str, j1(ny),j2(ny), ybegin(ny), yend(ny), jrev(ny), y, dy)

	   do j = j1(ny), j2(ny)
	      write(*,*) 'j2', j,y(j),str
	   enddo

        enddo

	iprint = 5 ! Number of BCs per line

!	Determine the x co-ordinates of the sloping section based on the y co-ordinates.

        idiff = j1_1-i2_2	! difference in i and j index at the start of the sloping section

	do i = j1_1, j1_2

	   x(i-idiff) = (y(i))/slope+xend(2)
	   dx(i-idiff) = x(i-idiff)-x(i-idiff-1)

	   write(*,*) 'sloping section', i, x(i-idiff), y(i)

	end do

!	Write out sloping section information in MFIX format

	islope_st = j1_1
        islope_end = j1_2 - 1

        ibc_offset = 10 - islope_st

	nseg = (islope_end-islope_st+1)/iprint
	nrem = mod(islope_end-islope_st+1,5)
	n1 = islope_st
	n2 = islope_st + iprint - 1

        ymin = 0.0d0
        zero = 0.0d0
        one = 1

	do n = 1,nseg

	   write(17,990) n1+ibc_offset, (x(i-idiff),i=n1,n2)
	   write(17,991) n1+ibc_offset, (x(i-idiff+1),i=n1,n2)
	   write(17,992) n1+ibc_offset, (ymin,i=n1,n2)
	   write(17,993) n1+ibc_offset, (y(i+1),i=n1,n2)
	   write(17,994) n1+ibc_offset
	   write(17,1093) n1+ibc_offset, (zero,i=n1,n2)
	   write(17,1094) n1+ibc_offset, (zero,i=n1,n2)
	   write(17,1095) n1+ibc_offset, (zero,i=n1,n2)
	   write(17,1096) n1+ibc_offset, (one,i=n1,n2)
	   write(17,1097) n1+ibc_offset, (zero,i=n1,n2)

	   n1 = n1 + iprint
	   n2 = n2 + iprint

	enddo

	if(nrem.ne.0) then

	   write(17,990) n1+ibc_offset, (x(i-idiff),i=n1,islope_end)
	   write(17,991) n1+ibc_offset, (x(i-idiff+1),i=n1,islope_end)
	   write(17,992) n1+ibc_offset, (ymin,i=n1,islope_end)
	   write(17,993) n1+ibc_offset, (y(i+1),i=n1,islope_end)
	   write(17,1093) n1+ibc_offset, (zero,i=n1,islope_end)
	   write(17,1094) n1+ibc_offset, (zero,i=n1,islope_end)
	   write(17,1095) n1+ibc_offset, (zero,i=n1,islope_end)
	   write(17,1096) n1+ibc_offset, (one,i=n1,islope_end)
	   write(17,1097) n1+ibc_offset, (zero,i=n1,islope_end)

	   if(nrem.eq.1) then
	      write(17,995) n1+ibc_offset
	   else if(nrem.eq.2) then
	      write(17,996) n1+ibc_offset
	   else if(nrem.eq.3) then
	      write(17,997) n1+ibc_offset
	   else if(nrem.eq.4) then
	      write(17,998) n1+ibc_offset
	   endif

	endif

!	Write out dx information in MFIX format

	do i = 1,j2_1-idiff

	   write(16,1000) i-1,dx(i)
	   write(*,*) 'dx', i-1, dx(i)

	enddo

!	Write out dy information in MFIX format

	do i = 1,j2_2

	   dy(i) = y(i)-y(i-1)
	   write(15,1001) i-1,dy(i)

	enddo

 990	FORMAT('#'/2X,'BC_X_w(',I3,')',5X,'=',2X, 5G12.5)
 991	FORMAT(2X,'BC_X_e(',I3,')',5X,'=',2X, 5G12.5)
 992	FORMAT(2X,'BC_Y_s(',I3,')',5X,'=',2X, 5G12.5)
 993	FORMAT(2X,'BC_Y_n(',I3,')',5X,'=',2X, 5G12.5)
 994	FORMAT(2X,'BC_TYPE(',I3,')',4X,'=',2X, 5(3x,"'NSW'",4x)/'#')
 995	FORMAT(2X,'BC_TYPE(',I3,')',4X,'=',2X, 1(3x,"'NSW'",4x)/'#')
 996	FORMAT(2X,'BC_TYPE(',I3,')',4X,'=',2X, 2(3x,"'NSW'",4x)/'#')
 997	FORMAT(2X,'BC_TYPE(',I3,')',4X,'=',2X, 3(3x,"'NSW'",4x)/'#')
 998	FORMAT(2X,'BC_TYPE(',I3,')',4X,'=',2X, 4(3x,"'NSW'",4x)/'#')
 1093	FORMAT(2X,'BC_Uw_s(',I3,',1)',2X,'=',2X, 5G12.5)
 1094	FORMAT(2X,'BC_Vw_s(',I3,',1)',2X,'=',2X, 5G12.5)
 1095	FORMAT(2X,'BC_Ww_s(',I3,',1)',2X,'=',2X, 5G12.5)
 1096	FORMAT(2X,'BC_JJ_PS(',I3,')',3X,'=',2X, 5I10)
 1097	FORMAT(2X,'BC_Thetaw_m(',I3,',1)',1X,'=',2X, 5G10.5)
 1000   FORMAT(2X,'DX (',I3,')',15X,'=',2X, G12.5)
 1001   FORMAT(2X,'DY (',I3,')',15X,'=',2X, G12.5)
	stop
	end

!
!       Subroutine which determines geometric stretch factor given initial & final locations
!       and number of points in any direction
!
	subroutine geometric(dx1, str, i1,i2, xbegin, xend, irev, x, dx)
        implicit real*8 (a-h,o-z)
	dimension x(0:500), dx(0:500)

        x(i1) = xbegin
        x(i2) = xend
        xlen1 = xend-xbegin

        str = 1.2d0
        rr = (xlen1-dx1)/dx1

! 	determine the stretch parameter for the geometric stretching using Newton's iteration

        do n = 1,100

	   str1 = str - (str**(i2-i1)-(rr+1)*str+rr) &
	   /((i2-i1)*str**(i2-(i1+1))-(rr+1))

	   write(*,*) 'stretch',str, str1

	   if(abs(str1-str).lt.1.0d-12) then
	      str = str1
	      write(*,*) 'converged,x1'
	      go to 102
	   endif

	   str = str1

        enddo

 102    continue

! 	Reverse the stretching depending on user's requirement

	if(irev.eq.1) then

	   do i = i1+1,i2

	      if(i.eq.i1+1) then
		 x(i) = x(i-1)+dx1
		 dx(i1) = dx1
	      else
		 x(i) = x(i-1)+str*(x(i-1)-x(i-2))
	      endif
	      dx(i) = x(i) - x(i-1)

	   enddo

        else

	   do i = i2-1, i1, -1

	      if(i.eq.i2-1) then
		 x(i) = x(i+1)-dx1
		 dx(i2) = dx1
	      else
		 x(i) = x(i+1)-str*(x(i+2)-x(i+1))
	      endif
	      dx(i) = x(i) - x(i-1)

	   enddo

        endif


	return
	end

