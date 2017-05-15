

      MODULE usr


        Use param
        Use param1


!
!       Declare the user-defined namelist variables (usrnlst.inc) in this module.
!       Also Include user-defined variables in this module.  To access the
!       variables from a subroutine add the statement "Use usr".  If allocatable
!       arrays are defined in this module allocate them in usr0.  To turn on the
!       user defined subroutines (usr0, usr1, and usr2) set call_usr to true in
!       mfix.dat.
!
!                       a dummy variable listed in usrnlst.inc
        DOUBLE PRECISION DUMMY_DP

! set .true. for entering tecplot output routine
        logical                               :: tecplot_output=.false.

! set .true. if 2D flow in XY plane only
        logical                               :: tec_no_k = .false.

! set .true. for writing out tecplot data in a traditional cell-centered
! format (BLOCK format):
        logical                               :: tec_output_block = .false.

! set .true. for writing out tecplot data at cell-center locations 
! (POINT format)
        logical                               :: tec_output_point = .false.

! set .true. if raw output (where calculated) is desired
        logical                               :: tec_output_raw = .false.

! temporary array for discretization error: de = p_g - mms_p_g ; etc.
        double precision, allocatable         :: de_ep_g(:)
        double precision, allocatable         :: de_p_g(:)
        double precision, allocatable         :: de_u_g(:)
        double precision, allocatable         :: de_v_g(:)
        double precision, allocatable         :: de_w_g(:)
        double precision, allocatable         :: de_t_g(:)
        double precision, allocatable         :: de_rop_s(:)
        double precision, allocatable         :: de_u_s(:)
        double precision, allocatable         :: de_v_s(:)
        double precision, allocatable         :: de_w_s(:)
        double precision, allocatable         :: de_t_s(:)
        double precision, allocatable         :: de_theta_m(:)


! norms of discretization errors
        double precision, allocatable         :: lnorms_ep_g(:)
        double precision, allocatable         :: lnorms_p_g(:)
        double precision, allocatable         :: lnorms_u_g(:)
        double precision, allocatable         :: lnorms_v_g(:)
        double precision, allocatable         :: lnorms_w_g(:)
        double precision, allocatable         :: lnorms_t_g(:)
        double precision, allocatable         :: lnorms_rop_s(:)
        double precision, allocatable         :: lnorms_u_s(:)
        double precision, allocatable         :: lnorms_v_s(:)
        double precision, allocatable         :: lnorms_w_s(:)
        double precision, allocatable         :: lnorms_t_s(:)
        double precision, allocatable         :: lnorms_theta_m(:)


! pressure shift value
        double precision                      :: delta_p_g


      END MODULE usr
