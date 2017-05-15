!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_SUBSTR                                             C
!  Purpose: Get a substring of data                                    C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_SUBSTR(STRING, L1, SUBSTR)
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) STRING, SUBSTR
      INTEGER       L, LMAX1, LMAX2, L1, Ls
      LOGICAL       FINISH
!
      SUBSTR = ' '
!
      LMAX1 = LEN(STRING)
      LMAX2 = LEN(SUBSTR)
      FINISH = .FALSE.
!
      Ls = 1
      DO 100 L = L1,LMAX1
        IF(STRING(L:L) .EQ. ',' .OR.&
           STRING(L:L) .EQ. ';'     )THEN     !terminate at comma or semicolon
          L1 = L + 1
          RETURN
        ELSEIF(STRING(L:L) .EQ. ' ' .AND. Ls .EQ. 1) THEN !ignore leading sp.
          CONTINUE
        ELSEIF(STRING(L:L) .EQ. ' ') THEN ! trailing sp. => next field
          FINISH = .TRUE.
        ELSEIF(FINISH) THEN    !character for next field
          L1 = L
          RETURN
        ELSE
          SUBSTR(Ls:Ls) = STRING(L:L)
          Ls = Ls+1
          IF(Ls .GT. LMAX2)RETURN
        ENDIF
100   CONTINUE
      L1 = LMAX1
      RETURN
      END
