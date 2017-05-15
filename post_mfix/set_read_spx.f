!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_READ_SPX(READ_SPX,VAR_INDEX,MAP_FILE,COUNT)        C
!  Purpose: Set READ_SPX and MAP_FILE, based on VAR_INDEX              C
!                                                                      C
!  Author: P. Nicoletti                               Date: 05-FEB-92  C
!  Reviewer:                                                           C
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
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_READ_SPX(READ_SPX,VAR_INDEX,MAP_FILE,COUNT)
!
      IMPLICIT NONE
!
! Passed arguments
!
      LOGICAL READ_SPX(*)
      INTEGER MAP_FILE(*)
      INTEGER VAR_INDEX
      INTEGER COUNT
!
      IF (VAR_INDEX.EQ.01) READ_SPX(1) = .TRUE.
      IF (VAR_INDEX.EQ.02) READ_SPX(2) = .TRUE.
      IF (VAR_INDEX.EQ.03) READ_SPX(2) = .TRUE.
      IF (VAR_INDEX.EQ.04) READ_SPX(3) = .TRUE.
      IF (VAR_INDEX.EQ.05) READ_SPX(3) = .TRUE.
      IF (VAR_INDEX.EQ.06) READ_SPX(3) = .TRUE.
      IF (VAR_INDEX.EQ.07) READ_SPX(4) = .TRUE.
      IF (VAR_INDEX.EQ.08) READ_SPX(4) = .TRUE.
      IF (VAR_INDEX.EQ.09) READ_SPX(4) = .TRUE.
      IF (VAR_INDEX.EQ.10) READ_SPX(5) = .TRUE.
      IF (VAR_INDEX.EQ.11) READ_SPX(6) = .TRUE.
      IF (VAR_INDEX.EQ.12) READ_SPX(6) = .TRUE.
      IF (VAR_INDEX.EQ.13) READ_SPX(6) = .TRUE.
      IF (VAR_INDEX.EQ.01) MAP_FILE(COUNT) = 1
      IF (VAR_INDEX.EQ.02) MAP_FILE(COUNT) = 2
      IF (VAR_INDEX.EQ.03) MAP_FILE(COUNT) = 2
      IF (VAR_INDEX.EQ.04) MAP_FILE(COUNT) = 3
      IF (VAR_INDEX.EQ.05) MAP_FILE(COUNT) = 3
      IF (VAR_INDEX.EQ.06) MAP_FILE(COUNT) = 3
      IF (VAR_INDEX.EQ.07) MAP_FILE(COUNT) = 4
      IF (VAR_INDEX.EQ.08) MAP_FILE(COUNT) = 4
      IF (VAR_INDEX.EQ.09) MAP_FILE(COUNT) = 4
      IF (VAR_INDEX.EQ.10) MAP_FILE(COUNT) = 5
      IF (VAR_INDEX.EQ.11) MAP_FILE(COUNT) = 6
      IF (VAR_INDEX.EQ.12) MAP_FILE(COUNT) = 6
      IF (VAR_INDEX.EQ.13) MAP_FILE(COUNT) = 6
!
      RETURN
      END
