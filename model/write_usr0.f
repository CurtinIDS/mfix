!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR0                                             C
!  Purpose: Write initial part of user-defined output                  C
!                                                                      C
!  Author: Shiv Meka                                  Date: dd-mmm-yy  C
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
      SUBROUTINE WRITE_USR0
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!
      RETURN
      END SUBROUTINE WRITE_USR0

      SUBROUTINE WRITE_CUST(TIMESTEP,DT_PRINT,field,particle_info,&
                                &mydat,mydat_size)
              use run,only: RUN_TYPE
              use discretelement, only: is_restart
              IMPLICIT NONE
              DOUBLE PRECISION, DIMENSION(mydat_size/3,3) :: mydat
              INTEGER, DIMENSION(mydat_size/3,5) :: particle_info
              CHARACTER(LEN=10) :: field
              CHARACTER(LEN=100) :: FILENAME
              INTEGER :: TIMESTEP,mydat_size,I,DT_PRINT
              LOGICAL :: is_file=0
              WRITE (FILENAME,'(I5.5)') &
                      &(int((TIMESTEP-1)/DT_PRINT))

!Timestep resets to 0 after every run
              INQUIRE(FILE=trim("Z_"//field(1:3))//'.'//&
                                 &adjustl(trim(FILENAME)),&
                                 &EXIST=is_file)
!Check if whether to create a new file 
              IF ((MOD(TIMESTEP-1,DT_PRINT).EQ.0&
                                &).OR.(TIMESTEP.EQ.1)) THEN
!The logic behind formatting integer to 5.5 is lifted from:
!http://computer-programming-forum.com/49-fortran/c5f91cf2ac9dab68.htm
                    IF (is_file) THEN
                        OPEN(unit=1, file=&
                                   &trim("Z_"//field(1:3))//'.'//&
                                   &adjustl(trim(FILENAME)//"_r"))
                    ELSE
                        OPEN(unit=1, file=&
                                   &trim("Z_"//field(1:3))//'.'//&
                                   &adjustl(trim(FILENAME)))
                    ENDIF
                    WRITE (1,*)'Timestep: ', TIMESTEP
                    WRITE (1,*)'PARTICLE #  X,Y,Z,',&
                                    &'CELL_X,CELL_Y,CELL_Z,PHASE'
                    DO I=1,mydat_size/3
                       IF (.NOT.particle_info(I,1).EQ.0) THEN
                       WRITE (1,*),particle_info(I,1),mydat(I,:),&
                                & particle_info(I,2:4), &
                                & particle_info(I,5)
                       END IF
                    ENDDO
                    WRITE (1,*),''
                    CLOSE(1)
              ENDIF
      END

