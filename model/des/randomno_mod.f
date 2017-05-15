! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Random Number Generation Utilities                     C
!  Purpose: Removed from interpolation mod and added built-in random
!           number routines instead of Pope's
!                                                                      C
!                                                                      C
!  Author: Sreekanth Pannala and Rahul Garg           Date: 23-Oct-08  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
    MODULE randomno

    USE constant

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: uni_rno, nor_rno


    CONTAINS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      SUBROUTINE UNI_RNO(Y)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      double precision, intent(out), dimension(:) :: y
      double precision rmean, variance, sigma
      integer i, nsize
!-----------------------------------------------

      nsize = size(y(:))

      call init_random_seed
      call random_number(y)

      rmean = sum(y(:))/nsize

!      write(*,*) 'Generating Uniform Random Variables for size', nsize
!      write(*,*) 'mean', rmean

      variance = 0.0
      do i = 1, nsize
!         write(20,*) i, y(i)
         variance = variance + (y(i)-rmean)**2
      end do

      close(20)

      variance = variance/nsize
      sigma = sqrt(variance)

!      write(*,*) 'sigma', sigma

      RETURN
      END SUBROUTINE UNI_RNO



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      SUBROUTINE NOR_RNO(Y, mean, sigma)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      double precision, intent(out), dimension(:) :: y
      double precision mean, sigma

      double precision lmean, lvariance, lsigma
      double precision x(2), w
      integer i, nsize, n
! no. of times this routine has been called
      integer, save :: COUNTER = 0
! so all components are written
      integer fileunit
!-----------------------------------------------
      COUNTER = COUNTER + 1
      fileunit = 20 + COUNTER

      nsize = size(y(:))

      call init_random_seed

      do i = 1, ceiling(real(nsize/2.0))
         do n = 1,100000
            call random_number(x)
            x = 2.0 * x - 1.0
            w = x(1)**2 + x(2)**2
            if(w.lt.1.0) exit
         end do

         w = sqrt( (-2.0 * log( w ) ) / w )
         y(2*i-1) = x(1) * w * sigma + mean
         if(2*i.lt.nsize) y(2*i) = x(2) * w * sigma + mean
      end do

      lmean = sum(y(:))/nsize

      !write(*,'(7X,A)') 'Generating Normal Random Variables'
      !write(*,'(7X,A,2X,ES15.5)') 'specified mean =', mean
      !write(*,'(7X,A,2X,ES15.5)') 'computed mean =', lmean

      !write(fileunit,'(A)') 'FROM NOR_RNO'
! specific to the call from init_particles_jn
      !write(fileunit,'(A,I5,A)') 'FOR DIRECTION = ', &
      !   COUNTER, ' where (1=X,2=Y,3=Z)'
      !write(fileunit,'(5X,A,5X,A)') 'particle no.', 'velocity component'

      lvariance = 0.0
      do i = 1, nsize
         !write(fileunit,'(I10,5X,ES15.5)') i, y(i)
         lvariance = lvariance + (y(i)-lmean)**2
      end do

      !close(fileunit)

      lvariance = lvariance/nsize
      lsigma = sqrt(lvariance)

      !write(*,'(7X,A,2X,ES15.5)') 'specified sigma =', sigma
      !write(*,'(7X,A,2X,ES15.5)') 'computed sigma =', lsigma

      RETURN
      END SUBROUTINE NOR_RNO



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      SUBROUTINE init_random_seed

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER              :: isize,idate(8)
      INTEGER,ALLOCATABLE  :: iseed(:)
!-----------------------------------------------

      CALL DATE_AND_TIME(VALUES=idate)
      CALL RANDOM_SEED(SIZE=isize)
      ALLOCATE( iseed(isize) )
      CALL RANDOM_SEED(GET=iseed)
      iseed = iseed * (idate(8)-500) ! idate(8) contains millisecond
      CALL RANDOM_SEED(PUT=iseed)

      DEALLOCATE( iseed )

      END SUBROUTINE init_random_seed


    END MODULE randomno

