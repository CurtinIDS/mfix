

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

! set .true. if 2D flow in XY plane only
        logical                               :: tec_no_k = .true.

! set .true. for entering tecplot output routine
        logical                               :: tecplot_output=.true.

! set .true. for writing out tecplot data in a traditional cell-centered
! format (BLOCK format):
        logical                               :: tec_output_block = .false.

! set .true. for writing out the gas temerature profile:
        logical                               :: t_g_profile = .true.

! exact solutions      
        double precision, allocatable         :: t_g_ex(:)

! norms of discretization errors        
        double precision, allocatable         :: lnorms_t_g(:)          

! x, y, z coordinates of the top-right corner of a cell.
! used to find the node locations in the mesh  
        double precision, allocatable         :: xtr(:)
        double precision, allocatable         :: ytr(:)
        double precision, allocatable         :: ztr(:)     

! temporary array for discretization error: de = t_g - t_g_ex ; etc.        
        double precision, allocatable         :: de_t_g(:)

      END MODULE usr
