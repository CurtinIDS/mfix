!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: read_thermochemical_mod.f                              C
!  Purpose: Read thermochemical data                                   C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: None                                C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

#include "version.inc"

MODULE read_thermochemical

  IMPLICIT NONE

! Full path to Burcat and Ruscic database
  CHARACTER(len=142) :: THERM = 'BURCAT.THR'

CONTAINS

!      Program Test; CALL Read_Therm_tester; END Program Test

      SUBROUTINE READ_Therm_tester
      Implicit none
      DOUBLE PRECISION Ahigh(7), Alow(7)
      DOUBLE PRECISION Thigh, Tlow, Tcom, MW
      DOUBLE PRECISION Cp1, h1, h2, Hf298oR
      CHARACTER(LEN=132) :: PATH
      CHARACTER(LEN=18) :: SPECIES
      integer funit, IER
      CHARACTER(len=255) FILENAME
      LOGICAL LocalCopy

      SPECIES = 'CH4'
      PATH = '.'
      funit = 5

      INQUIRE(FILE=TRIM(THERM),EXIST=LocalCopy)
      IF(LocalCopy)Then
        OPEN(CONVERT='BIG_ENDIAN',UNIT=funit,FILE=TRIM(THERM))
      ELSE
        FILENAME = './BURCAT.THR'
        OPEN(CONVERT='BIG_ENDIAN',UNIT=funit,FILE=TRIM(FILENAME), ERR=500)
      ENDIF

 !      Call Read_Therm(PATH, 'N2', Thigh, Tlow, Tcom, Ahigh, Alow, Hf298oR)
      Call Read_Therm(funit, SPECIES, Thigh, Tlow, Tcom, MW, Ahigh, &
         Alow, Hf298oR, IER)
      IF(IER /= 0) GOTO 200

      print *, SPECIES
      print *, Thigh, Tlow, Tcom, MW, Hf298oR*1.987207

!      print *, Hf298oR
!      T = 300
!      DO i = 1, 12
!        Cp1 = calc_CpoR(T, Thigh, Tlow, Tcom, Ahigh, Alow)*1.987207
!        T = T + 100
!        print *, T, Cp1
!      ENDDO

!      Cp1 = calc_CpoR(8D2, Thigh, Tlow, Tcom, Ahigh, Alow)*1.987207
!      h1 = calc_H0oR(4D2, Thigh, Tlow, Tcom, Ahigh, Alow)*1.987207
!      h2 = calc_H0oR(12D2, Thigh, Tlow, Tcom, Ahigh, Alow)*1.987207
      print *, Cp1, h1, h2
      CLOSE(UNIT=funit)
      ERROR_STOP
200   PRINT *, 'READ_Therm_tester: Species ', &
         TRIM(SPECIES), ' not found in Database!'
      ERROR_STOP
500   PRINT *, 'READ_Therm_tester: Cannot Open file ', TRIM(THERM), '!'
      PRINT *, 'Check path or copy mfix/model/thermochemical/', &
         TRIM(THERM), ' into run directory'
      ERROR_STOP
      END Subroutine READ_Therm_tester

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv  C
!                                                                        C
!     Module name: READ_Therm()                                          C
!     Purpose: Read Thermo coefficients from Burcat and Ruscic (2005)    C
!     Author: M. Syamlal                                 Date: 30-SEP-05 C
!                                                                        C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  C
!
      SUBROUTINE READ_Therm(funit, Sp, Thigh, Tlow, Tcom, MW, Ahigh, &
         Alow, Hf298oR, IER)
!
!-----------------------------------------------
!     M o d u l e s
!-----------------------------------------------

      IMPLICIT NONE

      CHARACTER(*) SP

!     holds one line in the input file
      CHARACTER(len=80) :: LINE_STRING

      CHARACTER(len=18) :: SPECIES, ss
      INTEGER i, funit, IER
      DOUBLE PRECISION Ahigh(7), Alow(7), Hf298oR
      DOUBLE PRECISION Thigh, Tlow, Tcom, MW


!     Associate a simple species name with that in BURCAT.THR file
      INTEGER, PARAMETER :: Max = 3
      CHARACTER(LEN=18), DIMENSION(2,Max) :: SPECIES_ALIAS
      SPECIES_ALIAS = RESHAPE ( (/ &
!        common name             BURCAT.THR name
!        123456789012345678     123456789012345678
        'O2                ',  'O2 REF ELEMENT    ', &
        'N2                ',  'N2  REF ELEMENT   ', &
        'CH4               ',  'CH4   ANHARMONIC  ' &
       /), (/2,Max/))


      IER = 0
      SPECIES = SP
      DO i = 1, MAX
        IF(TRIM(SP) == TRIM(SPECIES_ALIAS(1,I)))THEN
          SPECIES = SPECIES_ALIAS(2,I)
        ENDIF
      ENDDO


      LINE_STRING = '                '
      DO WHILE(LINE_STRING(1:11) /= 'THERMO DATA')
        READ(UNIT=funit, FMT='(A)',ERR=100,END=100)LINE_STRING
      END DO

      ss = '                 '
      call trimTab(SPECIES)
      DO WHILE(TRIM(ss) /= TRIM(SPECIES))
        READ(UNIT=funit, FMT='(A)',ERR=100,END=100)LINE_STRING
        ss = LINE_STRING(1:18)
        call trimTab(ss)

      END DO
!      print *, LINE_STRING


      call get_values(LINE_STRING, Tlow, Thigh, MW)

! Tcom is almost always 1000K, however there are a few species where
! this value is too high and causes a problem (e.g., liquid water).
! Therefore, set Tcom = Thigh when Thigh < 1000K.
      Tcom = min(1.0d3, Thigh)
      READ(UNIT=funit, FMT='(5E15.0)',ERR=300,END=300)Ahigh(1:5)
      READ(UNIT=funit, FMT='(5E15.0)',ERR=300,END=300)Ahigh(6:7), Alow(1:3)
      READ(UNIT=funit, FMT='(5E15.0)',ERR=300,END=300)Alow(4:7), Hf298oR



      RETURN
      ! species not found or THERMO DATA not found!
100   IER = 1
      RETURN

300   PRINT *, 'READ_Therm: Error reading coefficients for Species ', &
         TRIM(LINE_STRING(1:18))
      ERROR_STOP

      END SUBROUTINE READ_Therm

!**********************************************************************!
! Function: calc_CpoR                                                  !
! Purpose: Evaluate the polynomial form of the specific heat.          !
!                                                                      !
!**********************************************************************!
      PURE DOUBLE PRECISION FUNCTION  calc_CpoR(T, M, N)

! Polynomial coefficients
      use physprop, only: Ahigh  ! for T in [Tcom, Thigh]
      use physprop, only: Alow   ! for T in [Tlow, Tcom)
      use physprop, only: Thigh  ! Upper bound of use
      use physprop, only: Tlow   ! Lower bound of use
      use physprop, only: Tcom   ! Switch from low to high coeffs

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Evaluation temperaure (K)
      DOUBLE PRECISION, intent(in) :: T
! Phase index.
      INTEGER, intent(in) :: M
! Species index.
      INTEGER, intent(in) :: N

! Local Variables:
!-----------------------------------------------------------------------
! Bounded temperature.
      DOUBLE PRECISION :: xT

! Bound the temperature to the valid region.
      xT = max(min(T,Thigh(M,N)),Tlow(M,N))

! Evaluate the polynomial form.
      IF(T < Tcom(M,N))THEN
        calc_CpoR = calc_CpoR0(xT, Alow(1:5,M,N))
      ELSE
        calc_CpoR = calc_CpoR0(xT, Ahigh(1:5,M,N))
      ENDIF

      RETURN
      END Function calc_CpoR


!**********************************************************************!
! Function: calc_CpoR0                                                 !
! Purpose: Evaluate the polynomial form of the specific heat.          !
!                                                                      !
!**********************************************************************!
      PURE DOUBLE PRECISION FUNCTION  calc_CpoR0(T, A)

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Evaluation temperaure (K)
      DOUBLE PRECISION, intent(in) :: T
! Polynomial coefficients.
      DOUBLE PRECISION, intent(in) :: A(1:5)

! Evaluate the polynomial.
      calc_CpoR0 = (((A(5)*T +A(4))*T + A(3))*T + A(2))*T + A(1)

      RETURN
      END Function calc_CpoR0


!**********************************************************************!
! Function: calc_ICpoR                                                 !
! Purpose: Integrate the polynomial form of the specific heat.         !
!                                                                      !
!**********************************************************************!
      DOUBLE PRECISION FUNCTION calc_ICpoR(T, M, N, IER)

      use physprop, only: Ahigh
      use physprop, only: Thigh
      use physprop, only: ICpoR_h
      use physprop, only: Alow
      use physprop, only: Tlow
      use physprop, only: ICpoR_l
      use physprop, only: Tcom

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Evaluation temperaure (K)
      DOUBLE PRECISION, intent(in) :: T
! Phase index.
      INTEGER, intent(in) :: M
! Species index.
      INTEGER, intent(in) :: N
! Error Flag.
      INTEGER, intent(inout) :: IER

! Local Variables.
!-----------------------------------------------------------------------
      DOUBLE PRECISION :: xT
!-----------------------------------------------------------------------

! Initialize the bounded temperature and error flag.
      xT = T
      IER = 0

! Verify that the temperature is in a valid range.
      if(T > Thigh(M,N)) THEN
        xT = Thigh(M,N)
        IER = 101
      elseif(T < Tlow(M,N)) THEN
        xT = Tlow(M,N)
        IER = 102
      endif

! Integrate the polynomial from 0.0 to T.
      if (xT < Tcom(M,N)) then
        calc_ICpoR = calc_ICpoR0(xT, Alow(1:5,M,N),  ICpoR_l(M,N))
      else
        calc_ICpoR = calc_ICpoR0(xT, Ahigh(1:5,M,N), ICpoR_h(M,N))
      endif

      RETURN
      END FUNCTION calc_ICpoR


!**********************************************************************!
! Function: calc_ICpoR                                                 !
! Purpose: Integrate the polynomial form of the specific heat.         !
!                                                                      !
!**********************************************************************!
      PURE DOUBLE PRECISION FUNCTION calc_ICpoR0(T, A, REF_ICpoR)

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Evaluation temperaure (K)
      DOUBLE PRECISION, intent(in) :: T
! Polynomial coefficients.
      DOUBLE PRECISION, intent(in) :: A(1:5)
! Referece Integral
      DOUBLE PRECISION, intent(in) :: REF_ICpoR

! Local Variables.
!-----------------------------------------------------------------------
! Integral of specific heat polynomial (from 0 to T) over T
      DOUBLE PRECISION ICpoRT

!-----------------------------------------------------------------------

      ICpoRT = (((A(5)*T/5.0d0 + A(4)/4.0d0)*T + A(3)/3.0d0)*T +       &
         A(2)/2.0d0)*T + A(1)

      calc_ICpoR0 = T*ICpoRT - REF_ICpoR

      RETURN
      END FUNCTION calc_ICpoR0



!**********************************************************************!
! Function: calc_H0oR                                                  !
! Purpose: Calculate the heat of formation from the first six poly-    !
!          nomial coefficients.                                        !
!                                                                      !
! >>> This function is currently unused.                               !
!                                                                      !
!**********************************************************************!
      PURE DOUBLE PRECISION FUNCTION calc_H0oR(T, Th, Tl, Tc, Ah, Al)

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Polynomial coefficients. (High/Low)
      DOUBLE PRECISION, intent(in) :: Ah(7), Al(7)
! Temperature ranges of polynomials.
      DOUBLE PRECISION, intent(in) :: Th   ! Max temp (for Ahigh)
      DOUBLE PRECISION, intent(in) :: Tl   ! Min temp (for Alow)
      DOUBLE PRECISION, intent(in) :: Tc   ! switch from low to high
! Evaluation temperaure (K)
      DOUBLE PRECISION, intent(in) :: T

! Local Variables.
!-----------------------------------------------------------------------

      !ICp = calc_ICpoR(T, Th, Tl, Tc, Ah, Al)
      !If (T < Tc) then
      !  calc_H0oR = ICp + Al(6)
      !else
      !  calc_H0oR = ICp + Ah(6)
      !endif

      calc_H0oR = 0.0

      return
      END FUNCTION calc_H0oR


!**********************************************************************!
! SUBROUTINE: replaceTab                                               !
! Purpose: Replace all instances of a tab with a single space.         !
!**********************************************************************!
      SUBROUTINE replaceTab(C)

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Incoming string that will have tabs removed.
      CHARACTER(len=*) :: C

! Local Variables:
!-----------------------------------------------------------------------
! Loop counter
      INTEGER :: I

      DO I = 1, len(C)
        IF(C(I:I) == CHAR(9)) C(I:I)=' '
      ENDDO

      RETURN
      END SUBROUTINE replaceTab


!**********************************************************************!
! SUBROUTINE: trimTab                                                  !
! Purpose: Search a string for the first instance of a tab. The        !
!          location of the tab and all remaing string entries are      !
!          replaced with blank spaces.                                 !
!**********************************************************************!
      SUBROUTINE trimTab(C)

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Incoming string that will have tabs removed.
      CHARACTER(len=*) :: C

! Local Variables:
!-----------------------------------------------------------------------
! Loop counter
      INTEGER :: I
! Logical indicating that a tab was located.
      LOGICAL :: tabFound

! Initialize flag
      tabFound = .FALSE.

! Look at each entry of the string. Once a tab is located, the rest of
! the string is replaced by blank spaces.
      DO I = 1, len(C)
        IF(C(I:I) == CHAR(9) ) tabFound = .TRUE.
        if(tabFound) C(I:I)=' '
      ENDDO

      RETURN
      END SUBROUTINE trimTab

    END MODULE read_thermochemical
