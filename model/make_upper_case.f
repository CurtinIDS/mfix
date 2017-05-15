!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MAKE_UPPER_CASE (LINE_STRING,MAXCOL)                   C
!  Purpose: change lowercase characters to uppercase in input line     C
!                                                                      C
!  Author: P.Nicoletti                                Date: 26-NOV-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: A_UP, A_LO, Z_LO, A_DIFF, INT_C, L                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE MAKE_UPPER_CASE(LINE_STRING, MAXCOL)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                   input line to change to uppercase
      CHARACTER(len=*) LINE_STRING
!
!                   number of characters to look at in LINE_STRING
      INTEGER       MAXCOL
!
! local variables:
!
!                   ICHAR value for UPPERCASE A
      INTEGER       A_UP
!
!                   ICHAR value for lowercase a
      INTEGER       A_LO
!
!                   ICHAR value for lowercase z
      INTEGER       Z_LO
!
!                   ICHAR differnce between lower and uppercase letters
      INTEGER       A_DIFF
!
!                   holds ICHAR value of current character
      INTEGER       INT_C
!
!                   loop index
      INTEGER       L
!-----------------------------------------------
!
!
      A_UP = ICHAR('A')
      A_LO = ICHAR('a')
      Z_LO = ICHAR('z')
      A_DIFF = A_LO - A_UP
!
      DO L = 1, MAXCOL
         INT_C = ICHAR(LINE_STRING(L:L))
         IF (A_LO<=INT_C .AND. INT_C<=Z_LO) THEN
            INT_C = INT_C - A_DIFF
            LINE_STRING(L:L) = CHAR(INT_C)
         ENDIF
      END DO
      RETURN
      END SUBROUTINE MAKE_UPPER_CASE
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: REPLACE_TAB (LINE_STRING,MAXCOL)                       C
!  Purpose: replace tab characters with space                          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 10-JUL-03  C
!  Reviewer:                                          Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: A_UP, A_LO, Z_LO, A_DIFF, INT_C, L                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE REPLACE_TAB(LINE_STRING, MAXCOL)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
      CHARACTER, PARAMETER          :: TAB = CHAR(9)
      CHARACTER, PARAMETER          :: CRET = CHAR(13)

!                   input line to change to uppercase
      CHARACTER(len=*) LINE_STRING
!
!                   number of characters to look at in LINE_STRING
      INTEGER       MAXCOL
!
! local variables:
!
!                   loop index
      INTEGER       L
!-----------------------------------------------
!
!
!
      DO L = 1, MAXCOL
        if(LINE_STRING(L:L) .eq. TAB) LINE_STRING(L:L) = ' '
        if(LINE_STRING(L:L) .eq. CRET) LINE_STRING(L:L) = ' '
      END DO
      RETURN
      END SUBROUTINE REPLACE_TAB
