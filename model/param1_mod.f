! -*- f90 -*-
      MODULE param1

! Number of single precision .SPx files
      INTEGER, PARAMETER :: N_SPX = 11

! Maximum number of cell classes
      INTEGER, PARAMETER :: MAX_CLASS = 1000000
! Maximum number of corner cells
      INTEGER, PARAMETER :: MAX_NCORN = 4000

! Dimension for the upper triangle of an MxM matrix
      INTEGER :: DIMENSION_LM
! Total number of species
      INTEGER :: DIMENSION_N_all

! Parameters for testing if user input was specified.
      DOUBLE PRECISION, PARAMETER :: UNDEFINED = 9.87654321D31
      INTEGER, PARAMETER :: UNDEFINED_I = 987654321
      CHARACTER, PARAMETER :: UNDEFINED_C = ' '

! Cutoffs for large and small numbers
      DOUBLE PRECISION, PARAMETER :: LARGE_NUMBER = 1.0D32
      DOUBLE PRECISION, PARAMETER :: SMALL_NUMBER = 1.0D-15

! ZERO, HALF, ONE
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0
      DOUBLE PRECISION, PARAMETER :: HALF = 0.5d0
      DOUBLE PRECISION, PARAMETER :: ONE  = 1.0d0

   CONTAINS
      SUBROUTINE FILLER1
         IMPLICIT NONE
         ! empty subroutine so param is accessible from pymfix
      END SUBROUTINE FILLER1

      END MODULE param1
