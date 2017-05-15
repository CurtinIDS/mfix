program main
implicit none

  integer, parameter  :: dp = selected_real_kind(p=13,r=200)
  integer, parameter  :: ncases=5 ! max no. of cases
  integer, parameter  :: neq=12
  integer             :: ndata    ! no. of cases for which data is available

  integer :: i, j

  real(dp),dimension(ncases,neq)  :: l1de
  real(dp),dimension(ncases,neq)  :: l2de
  real(dp),dimension(ncases,neq)  :: linfde
  real(dp),dimension(neq)         :: ooal2, ooalinf ! temporary variables

  integer   :: io_status
  logical   :: eofflag=.FALSE.

  real(dp)  :: fsmall=1.0e-20_dp

  ! read de_norms data
  open(11,file='de_norms_collected.dat',status='unknown')
  do i=1,ncases+1
    do j=1,4
      read(11,*,iostat=io_status) ! read initial text lines
      if(IS_IOSTAT_END(io_status)) then
        ndata = i-1
        eofflag=.TRUE.
        exit
      endif
    enddo
    if(eofflag) exit

    read(11,*) (l1de(i,j),j=1,neq)
    read(11,*) (l2de(i,j),j=1,neq)
    read(11,*) (linfde(i,j),j=1,neq)
  enddo
  close(11)

  ! check ndata for validity
  if((ndata.le.1).or.(ndata.gt.5)) then
    write(*,*) "Inside ooa_test.f95. Check input data file."
  endif

  open(21,file='de_l2.dat',status='unknown')
  write(21,*) 'variables="h""pg""ug""vg""wg""us""vs""ws""tg""ts""epg""rops""ths"'
  do i = 1,ndata
    write(21,200) 2**(ncases-i),(l2de(i,j),j=1,neq)
  end do
  close(21)

  open(22,file='de_linf.dat',status='unknown')
  write(22,*) 'variables="h""pg""ug""vg""wg""us""vs""ws""tg""ts""epg""rops""ths"'
  do i = 1,ndata
    write(22,200) 2**(ncases-i),(linfde(i,j),j=1,neq)
  end do
  close(22)

  open(23,file='ooa_l2.dat',status='unknown')
  write(23,*) 'variables="h""pg""ug""vg""wg""us""vs""ws""tg""ts""epg""rops""ths"'
  do i = 2,ndata
    do j = 1,neq
      ooal2(j) = log((l2de(i-1,j)+fsmall)/(l2de(i,j)+fsmall))/log(2.0_dp)
    end do
    write(23,201) 2**(ncases-i),(ooal2(j),j=1,neq)
  end do
  close(23)

  open(24,file='ooa_linf.dat',status='unknown')
  write(24,*) 'variables="h""pg""ug""vg""wg""us""vs""ws""tg""ts""epg""rops""ths"'
  do i = 2,ndata
    do j = 1,neq
      ooalinf(j) = log((linfde(i-1,j)+fsmall)/(linfde(i,j)+fsmall))/log(2.0_dp)
    end do
    write(24,201) 2**(ncases-i),(ooalinf(j),j=1,neq)
  end do
  close(24)

200 Format (I3,30Es24.16)
201 Format (I3,30F8.3)

end program main
