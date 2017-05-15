MODULE utilities

  IMPLICIT NONE

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  function: mfix_isnan                                                !
!  Purpose: check whether argument is NAN                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      LOGICAL FUNCTION mfix_isnan(x)

! Dummy arguments
!---------------------------------------------------------------------//
      double precision, intent(in) :: x

! Local variables
!---------------------------------------------------------------------//
      CHARACTER(LEN=80) :: notnumber
!---------------------------------------------------------------------//

      mfix_isnan = .False.
      WRITE(notnumber,*) x
! To check for NaN's in x, see if x (a real number) contain a letter "N"
! "n" or symbol "?", in which case it is a NaN (Not a Number)

      IF(INDEX(notnumber,'?') > 0 .OR.     &
         INDEX(notnumber,'n') > 0 .OR.     &
         INDEX(notnumber,'N') > 0 ) THEN
        mfix_isnan = .TRUE.
         RETURN
      ENDIF

      RETURN
      END FUNCTION mfix_isnan

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  function: MAX_VEL_INLET                                             C
!  Purpose: Find maximum velocity at inlets.                           C
!                                                                      C
!  Author: S. Benyahia                                Date: 26-AUG-05  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION MAX_VEL_INLET()

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_defined, bc_type_enum
      use bc, only: mass_inflow, p_inflow
      use bc, only: bc_plane
      use bc, only: bc_k_b, bc_k_t, bc_j_s, bc_j_n, bc_i_w, bc_i_e
      use compar, only: dead_cell_at
      use fldvar, only: u_g, v_g, w_g, u_s, v_s, w_s
      use functions, only: funijk, is_on_mype_owns
      use functions, only: im_of, jm_of, km_of
      use mpi_utility, only: global_all_max
      use param, only: dimension_bc
      use param1, only: zero, small_number
      use physprop, only: mmax
      use run, only: units
      use toleranc, only: max_allowed_vel, max_inlet_vel, max_inlet_vel_fac
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: L, I, J, K, IJK, M
      DOUBLE PRECISION :: maxVEL
!---------------------------------------------------------------------//

! initializing
      maxVEL = ZERO

      DO L = 1, DIMENSION_BC

         IF(.NOT.(BC_DEFINED(L).AND.(BC_TYPE_ENUM(L)==MASS_INFLOW .OR. &
            BC_TYPE_ENUM(L) == P_INFLOW))) CYCLE

         DO K = BC_K_B(L), BC_K_T(L)
         DO J = BC_J_S(L), BC_J_N(L)
         DO I = BC_I_W(L), BC_I_E(L)

            IF (.NOT.IS_ON_myPE_OWNS(I,J,K)) CYCLE
            IF (DEAD_CELL_AT(I,J,K)) CYCLE
            IJK = FUNIJK(I,J,K)

            SELECT CASE (BC_PLANE(L))
            CASE ('E'); maxVEL = max(maxVEL,abs(U_G(IJK)))
            CASE ('N'); maxVEL = max(maxVEL,abs(V_G(IJK)))
            CASE ('T'); maxVEL = max(maxVEL,abs(W_G(IJK)))
            CASE ('W'); maxVEL = max(maxVEL,abs(U_G(IM_OF(IJK))))
            CASE ('S'); maxVEL = max(maxVEL,ABS(V_G(JM_OF(IJK))))
            CASE ('B'); maxVEL = max(maxVEL,abs(W_G(KM_OF(IJK))))
            END SELECT

            DO M=1,MMAX
               SELECT CASE (BC_PLANE(L))
               CASE ('E'); maxVEL = max(maxVEL,abs(U_s(IJK,M)))
               CASE ('N'); maxVEL = max(maxVEL,abs(V_s(IJK,M)))
               CASE ('T'); maxVEL = max(maxVEL,abs(W_s(IJK,M)))
               CASE ('W'); maxVEL = max(maxVEL,abs(U_s(IM_OF(IJK),M)))
               CASE ('S'); maxVEL = max(maxVEL,abs(V_s(JM_OF(IJK),M)))
               CASE ('B'); maxVEL = max(maxVEL,abs(W_s(KM_OF(IJK),M)))
               END SELECT
            ENDDO

         ENDDO
         ENDDO
         ENDDO
      ENDDO

      CALL GLOBAL_ALL_MAX(maxVEL, MAX_VEL_INLET)

! If no inlet velocity is specified, use an upper limit defined in
! toleranc_mod.f
      IF(MAX_VEL_INLET <= SMALL_NUMBER) THEN
         MAX_VEL_INLET = MAX_ALLOWED_VEL
         IF(UNITS == 'SI') MAX_INLET_VEL = 1D-2*MAX_ALLOWED_VEL
      ELSE
! Scale the value using a user defined scale factor
         MAX_VEL_INLET = 100.0d0*MAX_INLET_VEL_FAC*MAX_VEL_INLET
      ENDIF

      RETURN
      END FUNCTION MAX_VEL_INLET


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function: CHECK_VEL_BOUND()                                         C
!  Purpose: Check velocities upper bound to be less than speed of      C
!           sound                                                      C
!                                                                      C
!  Author: S. Benyahia                                Date: 25-AUG-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      LOGICAL FUNCTION CHECK_VEL_BOUND ()

! Modules
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
      use discretelement, only: des_continuum_coupled, des_continuum_hybrid
      use fldvar, only: ep_g, u_g, v_g, w_g
      use fldvar, only: ep_s, u_s, v_s, w_s
      use functions, only: fluid_at 
      use indices, only: i_of, j_of, k_of
      use mpi_utility, only: global_all_or
      use physprop, only: mmax
      use toleranc, only: max_inlet_vel

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: M
! Indices
      INTEGER :: IJK
      LOGICAL :: ALL_IS_ERROR
!---------------------------------------------------------------------//

!!$omp   parallel do private(IJK)
! initializing
      CHECK_VEL_BOUND = .FALSE.
      ALL_IS_ERROR    = .FALSE.

LOOP_FLUID : DO IJK = IJKSTART3, IJKEND3

         IF (FLUID_AT(IJK)) THEN
            IF(ABS(U_G(IJK)) > MAX_INLET_VEL .OR. &
               ABS(V_G(IJK)) > MAX_INLET_VEL .OR. &
               ABS(W_G(IJK)) > MAX_INLET_VEL) THEN
               CHECK_VEL_BOUND = .TRUE.
               WRITE(*,1000) MAX_INLET_VEL, I_OF(IJK), J_OF(IJK), K_OF(IJK), &
                             EP_g(IJK), U_G(IJK), V_G(IJK), W_G(IJK)
               EXIT LOOP_FLUID
            ENDIF

            IF (.NOT.DES_CONTINUUM_COUPLED .OR. DES_CONTINUUM_HYBRID) THEN
               DO M = 1, MMAX
                 IF(ABS(U_S(IJK,M)) > MAX_INLET_VEL .OR. &
                    ABS(V_S(IJK,M)) > MAX_INLET_VEL .OR. &
                    ABS(W_S(IJK,M)) > MAX_INLET_VEL) THEN
                   CHECK_VEL_BOUND = .TRUE.
                   WRITE(*,1010) MAX_INLET_VEL, I_OF(IJK), J_OF(IJK), K_OF(IJK), M, &
                                 EP_s(IJK, M), U_S(IJK,M), V_S(IJK,M), W_S(IJK,M)
                   EXIT LOOP_FLUID
                 ENDIF
               ENDDO
            ENDIF   ! end if(.not.des_continuum_coupled or des_continuum_hybrid)
         ENDIF

      ENDDO LOOP_FLUID

      CALL GLOBAL_ALL_OR(CHECK_VEL_BOUND, ALL_IS_ERROR)
      IF(ALL_IS_ERROR) CHECK_VEL_BOUND = .TRUE.

      RETURN
 1000 FORMAT(1X,'Message from: CHECK_VEL_BOUND',/&
            'WARNING: velocity higher than maximum allowed velocity: ', &
            G12.5, '(to change this adjust the scale factor MAX_INLET_VEL_FAC)'/&
            'in this cell: ','I = ',I4,2X,' J = ',I4,2X,' K = ',I4, /&
            '  ','Epg = ', G12.5, 'Ug = ', G12.5, 'Vg = ', G12.5, 'Wg = ', G12.5)
 1010 FORMAT(1X,'Message from: CHECK_VEL_BOUND',/&
            'WARNING: velocity higher than maximum allowed velocity: ', &
            G12.5,/&
            'in this cell: ','I = ',I4,2X,' J = ',I4,2X,' K = ',I4,' M = ',I4, /&
            '  ','Eps = ', G12.5,'Us = ', G12.5, 'Vs = ', G12.5, 'Ws = ', G12.5)

      END FUNCTION CHECK_VEL_BOUND


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Function: SEEK_COMMENT                                              !
!  Purpose: determine if (and where) a comment character appears       !
!           in a data input line                                       !
!     The function SEEK_COMMENT returns the index to where a comment   !
!     character was found in the input data line.  Equals MAXCOL + 1   !
!     if no-comment characters in the line.                            !
!                                                                      !
!  Author: P.Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      INTEGER FUNCTION SEEK_COMMENT (LINE, MAXCOL)
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! input data line
      CHARACTER(len=*) :: LINE
! maximum column of input data line to search
      INTEGER :: MAXCOL

! Local parameters
!---------------------------------------------------------------------//
! the number of designated comment characters
      INTEGER, PARAMETER :: DIM_COMMENT = 2

! Local variables
!---------------------------------------------------------------------//
! loop indicies
      INTEGER :: L, L2
! the comment characters
      CHARACTER, DIMENSION(DIM_COMMENT) :: COMMENT_CHAR
!---------------------------------------------------------------------//
      DATA COMMENT_CHAR/'#', '!'/

      DO L = 1, MAXCOL
         DO L2 = 1, DIM_COMMENT
            IF (LINE(L:L) == COMMENT_CHAR(L2)) THEN
               SEEK_COMMENT = L
               RETURN
            ENDIF
         END DO
      END DO
      SEEK_COMMENT = MAXCOL + 1
      RETURN
      END FUNCTION SEEK_COMMENT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Function: SEEK_END                                                  !
!  Purpose: determine where trailing blanks begin in a line            !
!     The function SEEK_END returns the index to where the last        !
!     character was found in the input data line.  Equals MAXCOL       !
!     if no trailing blank characters in the line.                     !
!                                                                      !
!  Author: P.Nicoletti, M. Syamlal                    Date: 7-AUG-92   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      INTEGER FUNCTION SEEK_END (LINE, MAXCOL)
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! maximum column of input data line to search
      INTEGER :: MAXCOL
! input data line
      CHARACTER :: LINE*(*)

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: L
!---------------------------------------------------------------------//
      SEEK_END = 0
      DO L = 1, MAXCOL
         IF (LINE(L:L) /= ' ') SEEK_END = L
      END DO
      RETURN
      END FUNCTION SEEK_END


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Function: LINE_TOO_BIG                                              !
!  Purpose: return an error condition if input data is located past    !
!           column MAXCOL in the data input file                       !
!     The function LINE_TOO_BIG returns a value greater than 0 to      !
!     indicate an error condition (data passed column MAXCOL in LINE)  !
!                                                                      !
!  Author: P.Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      INTEGER FUNCTION LINE_TOO_BIG (LINE, LINE_LEN, MAXCOL)
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! input data line
      CHARACTER(LEN=*) :: LINE
! length of input data line
      INTEGER :: LINE_LEN
! maximum column that non-blank charcaters are
! are in the input data line
      INTEGER :: MAXCOL

! Local variables
!---------------------------------------------------------------------//
! Loop index
      INTEGER :: L
!---------------------------------------------------------------------//
      DO L = MAXCOL + 1, LINE_LEN
         IF (LINE(L:L) /= ' ') THEN
            LINE_TOO_BIG = L
            RETURN
         ENDIF
      ENDDO
      LINE_TOO_BIG = 0
      RETURN
      END FUNCTION LINE_TOO_BIG


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Function: BLANK_LINE                                                 !
! Author: P. Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
! Purpose: Return .TRUE. if a line contains no input or only spaces.   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      LOGICAL FUNCTION BLANK_LINE (line)
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      CHARACTER :: LINE*(*)

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: L
!---------------------------------------------------------------------//
      BLANK_LINE = .FALSE.
      DO L=1, len(line)
         IF(line(L:L)/=' ' .and. line(L:L)/='    ')RETURN
      ENDDO
      BLANK_LINE = .TRUE.
      RETURN
      END FUNCTION BLANK_LINE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: bound_x                                                 !
!  Purpose: bound the values of x array                                !
!                                                                      !
!  Author: M. Syamlal                                 Date: 15-SEP-98  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE BOUND_X(ARRAY, IJKMAX2)

! Modules
!---------------------------------------------------------------------//
      USE param, only: dimension_3
      USE param1, only: one, zero
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Maximum dimension
      INTEGER :: IJKMAX2
! Array
      DOUBLE PRECISION :: Array(DIMENSION_3)
!--------------------------------------------------------------------//
      IF (IJKMAX2 > 0) THEN
         ARRAY(:) = MIN(ONE,MAX(ZERO,ARRAY(:)))
      ENDIF
      RETURN
      END SUBROUTINE BOUND_X

END MODULE utilities
