!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutinee: SET_PARAMETERS                                         !
!  Purpose: Set parameters used in array allocations.                  !
!                                                                      !
!  Author: J.Musser                                  Date: 17-APR-14   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_PARAMETERS

! Global Varialbes:
!---------------------------------------------------------------------//
! Number of TFM solids
      use physprop, only: MMAX
! Number of DES solids
      use discretelement, only: DES_MMAX
! Number of species
      use physprop, only: NMAX
! Number of scalar equations
      use scalars, only: NSCALAR

! Domain indices.
      use geometry, only: IMAX3, JMAX3, KMAX3, IJKMAX3

! Rank specific decompositions.
      use compar, only: IJKSIZE3_ALL
      USE compar, only: iStart3, iEnd3, iStart4, iEnd4
      USE compar, only: jStart3, jEnd3, jStart4, jEnd4
      USE compar, only: kStart3, kEnd3, kStart4, kEnd4

! Global Parameters:
!---------------------------------------------------------------------//
! Total number of solids phases
      use param, only: DIMENSION_M
! Maximum number of species.
      use param, only: DIMENSION_N_g ! Gas
      use param, only: DIMENSION_N_s ! Solids
! Maximum number of user-defined scalars
      use param, only: DIMENSION_SCALAR
      use param, only: DIM_SCALAR2
! Dimension for the upper triangle of an MxM matrix
      use param1, only: DIMENSION_LM
! Total number of species
      use param1, only: DIMENSION_N_all
! Parameter constants.
      use param1, only: UNDEFINED_I
! Axis decomposition
      USE param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K
      USE param, only: DIMENSION_3, DIMENSION_4
      USE param, only: DIMENSION_3G, DIMENSION_3L, DIMENSION_3P

! MPI-Domain decompoint and rank flags.
      use compar, only: myPE
! Flag for POST_MFIX
      use cdist, only: bDoing_postmfix


      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop index
      INTEGER :: M, MMAX_TOT
!......................................................................!

! The total number of solids.
      DIMENSION_M = MAX(1, merge(MMAX,MMAX+DES_MMAX,bDoing_postmfix))

! Number of gas phase species.
      DIMENSION_N_g = merge(NMAX(0), 1, NMAX(0) /= UNDEFINED_I)

! Number of solids phase species. (Max over all phases)
      DIMENSION_N_s = 1

      MMAX_TOT = merge(MMAX, MMAX+DES_MMAX, bDoing_postMFIX)
      DO M = 1, MMAX_TOT
         IF(NMAX(M) /= UNDEFINED_I) &
            DIMENSION_N_s = max(DIMENSION_N_s, NMAX(M))
      ENDDO

! Max number of species over all phases
      DIMENSION_N_all = max(DIMENSION_N_g, DIMENSION_N_s)

! Size of MxM upper triangular matrix.
      DIMENSION_LM = (DIMENSION_M * (DIMENSION_M-1)/2)+1

! Number of scalar equations.
      DIMENSION_SCALAR = NSCALAR
      DIM_SCALAR2 = 2*NSCALAR



! Set DIMENSION_x variables.
      DIMENSION_I = IMAX3
      DIMENSION_J = JMAX3
      DIMENSION_K = KMAX3

      DIMENSION_3 = (kEnd3-kStart3+1)*(jEnd3-jStart3+1)*(iEnd3-iStart3+1)
      DIMENSION_4 = (kEnd4-kStart4+1)*(jEnd4-jStart4+1)*(iEnd4-iStart4+1)

      DIMENSION_3G = IJKMAX3            ! Global IJK array
      DIMENSION_3L = IJKSIZE3_ALL(myPE) ! Local IJK array
      DIMENSION_3P = merge(DIMENSION_3, 1, .NOT.bDoing_PostMFIX)


      RETURN
      END SUBROUTINE SET_PARAMETERS

