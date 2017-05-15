

      MODULE usr


        Use param
        Use param1


!
!       Include user-defined variables in this module.  To access the variables
!       from a subroutine add the statement "Use usr".  If allocatable arrays
!       are defined in this module allocate them in usr0.  To turn on the
!       user defined subroutines (usr0, usr1, and usr2) set call_usr to true in
!       mfix.dat.
!
!
!                      Sherwood number
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: N_Sh
!
!                      Function of ash-layer void fraction
      DOUBLE PRECISION f_EP_A

      END MODULE usr
