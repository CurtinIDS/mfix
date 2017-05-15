!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: IN_BIN_512                                             C
!  Purpose: read an array in chunks of 512 bytes    (DP WORDS)         C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
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
!  Local variables: NWORDS, DS, L, NSEG, NREM, LC, N1, N2              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE IN_BIN_512(IUNIT,ARRAY,N,NEXT_REC)
!
!
      Use machine

      IMPLICIT NONE
!
! passed arguments
!
!                      array to write out
      DOUBLE PRECISION ARRAY(*)
!
!                      output unit number
      INTEGER          IUNIT
!
!                      number of elements in ARRAY
      INTEGER          N
!
!                      next record number in direct access output file
      INTEGER          NEXT_REC
!
! local variables
!
!                      number of words for 512 bytes
      INTEGER          NWORDS
!
!                      loop counter
      INTEGER          L
!
!                      number of full 512 byte segments need to write N
!                      double precision words
      INTEGER          NSEG
!
!                      number of double precision words in the partially
!                      filled last record
      INTEGER          NREM
!
!                      loop counter
      INTEGER          LC
!
!                      write out array elements N1 to N2
      INTEGER          N1 , N2
!
      NWORDS = NWORDS_DP
      IF (N.LE.NWORDS) THEN
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      END IF
!
      NSEG = N / NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
!
! read the full 512 byte segments
!
      DO 100 LC = 1,NSEG
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
100   CONTINUE
!
! read the partially filled last record
!
      IF (NREM.NE.0) THEN
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      END IF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: IN_BIN_512I                                            C
!  Purpose: read an array in chunks of 512 bytes (INTEGER WORDS)       C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
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
!  Local variables: NWORDS, L, NSEG, NREM, LC, N1, N2                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE IN_BIN_512I(IUNIT,ARRAY,N,NEXT_REC)
!
!
      Use machine
      IMPLICIT NONE
!
! passed arguments
!
!                      array to write out
      INTEGER          ARRAY(*)
!
!                      output unit number
      INTEGER          IUNIT
!
!                      number of elements in ARRAY
      INTEGER          N
!
!                      next record number in direct access output file
      INTEGER          NEXT_REC
!
! local variables
!
!                      number of words for 512 bytes (nwords * 4 = 512)
      INTEGER          NWORDS
!
!                      loop counter
      INTEGER          L
!
!                      number of full 512 byte segments need to write N
!                      double precision words
      INTEGER          NSEG
!
!                      number of double precision words in the partially
!                      filled last record
      INTEGER          NREM
!
!                      loop counter
      INTEGER          LC
!
!                      write out array elements N1 to N2
      INTEGER          N1 , N2
!
      NWORDS = NWORDS_I
      IF (N.LE.NWORDS) THEN
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      END IF
!
      NSEG = N / NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
!
! read the full 512 byte segments
!
      DO 100 LC = 1,NSEG
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
100   CONTINUE
!
! read the partially filled last record
!
      IF (NREM.NE.0) THEN
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      END IF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: IN_BIN_512R                                            C
!  Purpose: read in an array in chunks of 512 bytes (REAL    WORDS)    C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
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
!  Local variables: NWORDS, L, NSEG, NREM, LC, N1, N2                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE IN_BIN_512R(IUNIT,ARRAY,N,NEXT_REC)
!
!
      Use machine
      IMPLICIT NONE
!
! passed arguments
!
!                      array to write out
      REAL             ARRAY(*)
!
!                      output unit number
      INTEGER          IUNIT
!
!                      number of elements in ARRAY
      INTEGER          N
!
!                      next record number in direct access output file
      INTEGER          NEXT_REC
!
! local variables
!
!                      number of words for 512 bytes (nwords * 4 = 512)
      INTEGER          NWORDS
!
!                      loop counter
      INTEGER          L
!
!                      number of full 512 byte segments need to write N
!                      double precision words
      INTEGER          NSEG
!
!                      number of double precision words in the partially
!                      filled last record
      INTEGER          NREM
!
!                      loop counter
      INTEGER          LC
!
!                      write out array elements N1 to N2
      INTEGER          N1 , N2
!
      NWORDS = NWORDS_R
      IF (N.LE.NWORDS) THEN
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      END IF
!
      NSEG = N / NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
!
! write out the full 512 byte segments
!
      DO 100 LC = 1,NSEG
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
100   CONTINUE
!
! write out the partially filled last record
!
      IF (NREM.NE.0) THEN
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      END IF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: IN_BIN_R(IUNIT,ARRAY,IJKMAX2,NEXT_REC)                 C
!  Purpose: read in a time-dependent restart variable (REAL)           C
!                                                                      C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
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
!  Local variables: LC, ARRAY_REAL                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE IN_BIN_R(IUNIT,ARRAY,IJKMAX2,NEXT_REC)
!
      Use param
      Use param1
      IMPLICIT NONE
!
! passed arguments
!
!                      double precision array to write out
      DOUBLE PRECISION ARRAY(*)
!
!                      unit number to write to
      INTEGER          IUNIT
!
!                      record pointer into file IUNIT
      INTEGER          NEXT_REC
!
!                      number of indices in ARRAY to write out
      INTEGER          IJKMAX2
!
! local variables
!
!                      single precision version of ARRAY
      REAL             ARRAY_REAL(DIMENSION_3)
!
!                      loop counter
      INTEGER          LC
!
      CALL IN_BIN_512R (IUNIT,ARRAY_REAL,IJKMAX2,NEXT_REC)
!
      DO 100 LC = 1,IJKMAX2
         ARRAY(LC) = ARRAY_REAL(LC)
100   CONTINUE
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_FILEP(RUN_NAME, RUN_TYPE, NO_FILES)               C
!  Purpose: open all the files for this run (modified for POST_MFIX)   C
!                                                                      C
!  Author: P. Nicoletti                               Date: 12-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
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
!  Local variables: EXT, FILE_NAME, LC, NB                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      LOGICAL FUNCTION OPEN_FILEP(RUN_NAME, RUN_TYPE, NO_FILES)
!
      Use param
      Use param1
      Use post3d
      Use machine
      Use funits
      Use compar
!
      IMPLICIT NONE
!
! passed arguments
!
!                   run_name (as specified in input file)
      CHARACTER(LEN=*) RUN_NAME
!
!                   run type
      CHARACTER(LEN=*) RUN_TYPE
!
! local variables
!
!                   extension to filename
      CHARACTER(LEN=4) :: EXT
!
!                   run_name + extension
      CHARACTER(LEN=64) :: FILE_NAME
!
!                   Loop counter
      INTEGER       LC
!
!                   index to first blank character in run_name
      INTEGER       NB
!
!                   Number of files to be opened
      INTEGER       NO_FILES
      CHARACTER(LEN=35) :: EXT_END
!-----------------------------------------------

      ext_end = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!
      OPEN_FILEP = .FALSE.
!
! DETERMINE THE FIRST BLANK CHARCATER IN RUN_NAME
!
      DO 100 LC = 1,LEN(RUN_NAME)
         IF (RUN_NAME(LC:LC).EQ.' ') THEN
            NB = LC
            GOTO 125
         END IF
100   CONTINUE
      WRITE (*,*) 'RUN_NAME TOOOOOOO LOOOONG'
      RETURN
!
125   IF (NB+4.GT.LEN(FILE_NAME)) THEN
         WRITE (*,*) 'RUN_NAME TOOOOOOO LOOOONG'
         RETURN
      END IF
!
!  Open RES file
!
      EXT = '.RES'
      FILE_NAME          = ' '
      FILE_NAME(1:NB-1)  = RUN_NAME(1:NB-1)
      FILE_NAME(NB:NB+3) = EXT(1:4)
      IF(RUN_TYPE .EQ. 'NEW')THEN
        OPEN (UNIT=UNIT_RES,FILE=FILE_NAME,STATUS='NEW',RECL=OPEN_N1,&
            ACCESS='DIRECT',FORM='UNFORMATTED',ERR=300,CONVERT='BIG_ENDIAN')
      ELSE
        OPEN (UNIT=UNIT_RES,FILE=FILE_NAME,STATUS='OLD',RECL=OPEN_N1,&
            ACCESS='DIRECT',FORM='UNFORMATTED',ERR=300,CONVERT='BIG_ENDIAN')
      ENDIF
      IF(NO_FILES .EQ. 0) THEN
        OPEN_FILEP = .TRUE.
        RETURN
      ENDIF
!
! Open SPx files
!
      EXT = '.SPx'
      DO 200 LC = 1,N_SPX
        ext(4:4) = ext_end(LC:LC)
        FILE_NAME          = ' '
        FILE_NAME(1:NB-1)  = RUN_NAME(1:NB-1)
        FILE_NAME(NB:NB+3) = EXT(1:4)
        IF(RUN_TYPE .EQ. 'NEW') THEN
          OPEN (UNIT=UNIT_SPX+LC,FILE=FILE_NAME,STATUS='NEW',&
                RECL=OPEN_N1,ACCESS='DIRECT',FORM='UNFORMATTED',ERR=150,CONVERT='BIG_ENDIAN')
        ELSE
          OPEN (UNIT=UNIT_SPX+LC,FILE=FILE_NAME,STATUS='OLD',&
                RECL=OPEN_N1,ACCESS='DIRECT',FORM='UNFORMATTED',ERR=150,CONVERT='BIG_ENDIAN')
        ENDIF
        SPX_OPEN(LC) = .TRUE.
        GOTO 200
150     WRITE (*,*) 'ERROR OPENING FILE: ', FILE_NAME
        SPX_OPEN(LC) = .FALSE.
200   CONTINUE
      OPEN_FILEP = .TRUE.
      RETURN
300   WRITE (*,*) 'ERROR OPENING FILE: ', FILE_NAME
      RETURN
1000  FORMAT(I1)
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_ARRAY (ARRAY, MESSAGE)                             C
!  Purpose: print out a 3D array to standard output                    C
!                                                                      C
!  Author: P.Nicoletti                                Date: 02-DEC-91  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: KMAX2                                         C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: POINTER, LK                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE OUT_ARRAY (ARRAY,MESSAGE)
!
      Use param
      Use param1
      Use geometry
      Use fldvar
      Use physprop
      Use indices
      Use funits
      Use compar
      Use functions

      IMPLICIT NONE
!
! passed arguments
!
!                       array to print out
      DOUBLE PRECISION  ARRAY(*)
!
!                       message to print out
      CHARACTER(LEN=*) :: MESSAGE
!
! local variables
!
!                       pointer into array (points to start of a k-plane)
      INTEGER           IJK
!
!                       loop counter
      INTEGER           L

      DO 100 L = 1,KMAX2
         IJK = FUNIJK (1,1,L)
         WRITE (UNIT_OUT,1100) MESSAGE , L
         CALL OUT_ARRAY_K (ARRAY(IJK))
100   CONTINUE
!
1100  FORMAT(/,1X,A,' at K = ' ,I4,/)
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_ARRAY_C (ARRAY,MESSAGE)                            C
!  Purpose: print out a 3D array to standard output (character)        C
!                                                                      C
!  Author: P.Nicoletti                                Date: 10-JAN-92  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: KMAX2                                         C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: POINTER, LK                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE OUT_ARRAY_C(ARRAY,MESSAGE)
!
      Use param
      Use param1
      Use geometry
      Use fldvar
      Use physprop
      Use indices
      Use funits
      Use compar
      Use functions

      IMPLICIT NONE
!
! passed arguments
!
!                       array to print out
      CHARACTER(LEN=3) ::       ARRAY(*)
!
!                       message to print out
      CHARACTER(LEN=*)     MESSAGE
!
! local variables
!
!                       pointer into array (points to start of a k-plane)
      INTEGER           IJK
!
!                       loop counter
      INTEGER           L
!
      DO 100 L = 1,KMAX2
         IJK = FUNIJK (1,1,L)
         WRITE (UNIT_OUT,1100) MESSAGE , L
         CALL OUT_ARRAY_KC (ARRAY(IJK))
100   CONTINUE
!
1100  FORMAT(/,1X,A,' at K = ' ,I4,/)
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_ARRAY_K(ARRAY)                                     C
!  Purpose: print out a 2D (constant k-plane) array to standard output C
!                                                                      C
!  Author: P.Nicoletti                                Date: 02-DEC-91  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX2, JMAX2                                  C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables: NCOL, NTAB, L1, L2, L3, IFORM1, IFORM2, IJ1, IJ2   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OUT_ARRAY_K (ARRAY)
!
      Use param
      Use param1
      Use geometry
      Use fldvar
      Use physprop
      Use indices
      Use funits
      Use compar
      Use functions

      IMPLICIT NONE
!
! passed arguments
!
!                      2D array to print out
      DOUBLE PRECISION ARRAY(*)
!
! local variables
!
!                      number of columns to print out across the page
      INTEGER          NCOL
!
!                      number of tables the 2D array must be split into
!                      for printing
      INTEGER          NTAB
!
!                      loop indices
      INTEGER          L1x, L2x, L3, IJK
!
!                      start and end 'I' for current table
      INTEGER          IFORM1 , IFORM2
!
!                      start 'IJ' and end 'IJ' for a given 'J' to print out
      INTEGER          IJ1 , IJ2
!
! NOTE:  IF NCOL IS CHANGED TO A NUMBER GREATER THAN 30, THEN THE "30"
!        IN FORMATS 5050 AND 5100 MUST BE CHANGED TO THAT NUMBER.
!
      NCOL = 10
      NTAB = IMAX2 / NCOL   +  1
      IF ( MOD(IMAX2,NCOL).EQ.0 ) NTAB = NTAB - 1
!
      DO 100 L1x = 1,NTAB
         IFORM1 = 1 + NCOL*(L1x-1)
         IFORM2 = NCOL * L1x
         IFORM2 = MIN(IFORM2,IMAX2)
         WRITE (UNIT_OUT,5050) (L3,L3=IFORM1,IFORM2)
         DO 50 L2x = JMAX2,1,-1
            IJ1 = FUNIJK(IFORM1,L2x,1)
            IJ2 = FUNIJK(IFORM2,L2x,1)
            WRITE (UNIT_OUT,5100) L2x , (ARRAY(L3),L3=IJ1,IJ2)
50       CONTINUE
100   CONTINUE
!
5050  FORMAT (3X,'J',3X,'I=',3X,10(I3,9X))
5100  FORMAT (1X,I3,3X,10(1PE12.4))
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_ARRAY_KC (ARRAY)                                   C
!  Purpose: print out a 2D (constant k-plane) array to standard output C
!           (character)                                                C
!                                                                      C
!  Author: P.Nicoletti                                Date: 02-DEC-91  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX2, JMAX2                                  C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: NCOL, NTAB, L1, L2, L3, IFORM1, IFORM2, IJ1, IJ2   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OUT_ARRAY_KC (ARRAY)
!
      Use param
      Use param1
      Use geometry
      Use fldvar
      Use physprop
      Use indices
      Use funits
      Use compar
      Use functions

      IMPLICIT NONE
!
! passed arguments
!
!                      2D array to print out
      CHARACTER(LEN=3) :: ARRAY(*)
!
! local variables
!
!                      number of columns to print out across the page
      INTEGER          NCOL
!
!                      number of tables the 2D array must be split into
!                      for printing
      INTEGER          NTAB
!
!                      loop indices
      INTEGER          L1x, L2x, L3, IJK
!
!                      start and end 'I' for current table
      INTEGER          IFORM1 , IFORM2
!
!                      start 'IJ' and end 'IJ' for a given 'J' to print out
      INTEGER          IJ1 , IJ2
!
! NOTE:  IF NCOL IS CHANGED TO A NUMBER GREATER THAN 30, THEN THE "30"
!        IN FORMATS 5050 AND 5100 MUST BE CHANGED TO THAT NUMBER.
!
      NCOL = 30
      NTAB = IMAX2 / NCOL   +  1
      IF ( MOD(IMAX2,NCOL).EQ.0 ) NTAB = NTAB - 1
!
      DO 100 L1x = 1,NTAB
         IFORM1 = 1 + NCOL*(L1x-1)
         IFORM2 = NCOL * L1x
         IFORM2 = MIN(IFORM2,IMAX2)
         WRITE (UNIT_OUT,5050) (L3,L3=IFORM1,IFORM2)
         DO 50 L2x = JMAX2,1,-1
            IJ1 = FUNIJK(IFORM1,L2x,1)
            IJ2 = FUNIJK(IFORM2,L2x,1)
            WRITE (UNIT_OUT,5100) L2x , (ARRAY(L3),L3=IJ1,IJ2)
50       CONTINUE
100   CONTINUE
!
5050  FORMAT (3X,'J',3X,'I=',3X,30(I3,1X))
5100  FORMAT (1X,I3,8X,30(A3,1X))
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_BIN_512                                            C
!  Purpose: write an array in chunks of 512 bytes    (DP WORDS)        C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
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
!  Local variables: NWORDS, DS, L, NSEG, NREM, LC, N1, N2              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE OUT_BIN_512(IUNIT,ARRAY,N,NEXT_REC)
!
      Use machine
!
      IMPLICIT NONE
!
! passed arguments
!
!                      array to write out
      DOUBLE PRECISION ARRAY(*)
!
!                      output unit number
      INTEGER          IUNIT
!
!                      number of elements in ARRAY
      INTEGER          N
!
!                      next record number in direct access output file
      INTEGER          NEXT_REC
!
! local variables
!
!                      number of words for 512 bytes
      INTEGER          NWORDS
!
!                      loop counter
      INTEGER          L
!
!                      number of full 512 byte segments need to write N
!                      double precision words
      INTEGER          NSEG
!
!                      number of double precision words in the partially
!                      filled last record
      INTEGER          NREM
!
!                      loop counter
      INTEGER          LC
!
!                      write out array elements N1 to N2
      INTEGER          N1 , N2
!
      NWORDS = NWORDS_DP
      IF (N.LE.NWORDS) THEN
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      END IF
!
      NSEG = N / NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
!
! read the full 512 byte segments
!
      DO 100 LC = 1,NSEG
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
100   CONTINUE
!
! read the partially filled last record
!
      IF (NREM.NE.0) THEN
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      END IF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_BIN_512I                                           C
!  Purpose: write out an array in chunks of 512 bytes (INTEGER WORDS)  C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
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
!  Local variables: NWORDS, L, NSEG, NREM, LC, N1, N2                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE OUT_BIN_512I(IUNIT,ARRAY,N,NEXT_REC)
!
      Use machine
!
      IMPLICIT NONE
!
! passed arguments
!
!                      array to write out
      INTEGER          ARRAY(*)
!
!                      output unit number
      INTEGER          IUNIT
!
!                      number of elements in ARRAY
      INTEGER          N
!
!                      next record number in direct access output file
      INTEGER          NEXT_REC
!
! local variables
!
!                      number of words for 512 bytes (nwords * 4 = 512)
      INTEGER          NWORDS
!
!                      loop counter
      INTEGER          L
!
!                      number of full 512 byte segments need to write N
!                      double precision words
      INTEGER          NSEG
!
!                      number of double precision words in the partially
!                      filled last record
      INTEGER          NREM
!
!                      loop counter
      INTEGER          LC
!
!                      write out array elements N1 to N2
      INTEGER          N1 , N2
!
      NWORDS = NWORDS_I
      IF (N.LE.NWORDS) THEN
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      END IF
!
      NSEG = N / NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
!
! write out the full 512 byte segments
!
      DO 100 LC = 1,NSEG
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
100   CONTINUE
!
! write out the partially filled last record
!
      IF (NREM.NE.0) THEN
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      END IF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_BIN_512R                                           C
!  Purpose: write out an array in chunks of 512 bytes (REAL    WORDS)  C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
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
!  Local variables: NWORDS, L, NSEG, NREM, LC, N1, N2                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE OUT_BIN_512R(IUNIT,ARRAY,N,NEXT_REC)
!
      Use machine
!
      IMPLICIT NONE
!
! passed arguments
!
!                      array to write out
      REAL             ARRAY(*)
!
!                      output unit number
      INTEGER          IUNIT
!
!                      number of elements in ARRAY
      INTEGER          N
!
!                      next record number in direct access output file
      INTEGER          NEXT_REC
!
! local variables
!
!                      number of words for 512 bytes (nwords * 4 = 512)
      INTEGER          NWORDS
!
!                      loop counter
      INTEGER          L
!
!                      number of full 512 byte segments need to write N
!                      double precision words
      INTEGER          NSEG
!
!                      number of double precision words in the partially
!                      filled last record
      INTEGER          NREM
!
!                      loop counter
      INTEGER          LC
!
!                      write out array elements N1 to N2
      INTEGER          N1 , N2
!
      NWORDS = NWORDS_R
      IF (N.LE.NWORDS) THEN
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      END IF
!
      NSEG = N / NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
!
! write out the full 512 byte segments
!
      DO 100 LC = 1,NSEG
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
100   CONTINUE
!
! write out the partially filled last record
!
      IF (NREM.NE.0) THEN
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      END IF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_BIN_R                                              C
!  Purpose: write out a time-dependent restart variable (REAL)         C
!                                                                      C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
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
!  Local variables: LC, ARRAY_REAL                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE OUT_BIN_R(IUNIT,ARRAY,IJKMAX2,NEXT_REC)
      Use param
      Use param1
!
      IMPLICIT NONE
!
! passed arguments
!
!                      double precision array to write out
      DOUBLE PRECISION ARRAY(*)
!
!                      unit number to write to
      INTEGER          IUNIT
!
!                      record pointer into file IUNIT
      INTEGER          NEXT_REC
!
!                      number of indices in ARRAY to write out
      INTEGER          IJKMAX2
!
! local variables
!
!                      single precision version of ARRAY
      REAL             ARRAY_REAL(DIMENSION_3)
!
!                      loop counter
      INTEGER          LC
!
      DO 100 LC = 1,IJKMAX2
         ARRAY_REAL(LC) = SNGL(ARRAY(LC))
100   CONTINUE
!
      CALL OUT_BIN_512R (IUNIT,ARRAY_REAL,IJKMAX2,NEXT_REC)
!
      RETURN
      END
!
!
      subroutine convert_from_io_dp(arr_io,arr_internal,n)
!
      use geometry
      use indices
      use compar
      use functions
!
      implicit none
!
      double precision   arr_io(*) , arr_internal(*)
      integer            n,i,j,k,ijk,ijk_io
!
      do k = 1,kmax2
         do j = 1,jmax2
            do i = 1,imax2
               ijk    = funijk(i,j,k)
               ijk_io = funijk_io(i,j,k)
               arr_internal(ijk) = arr_io(ijk_io)
            end do
         end do
      end do
!
      return
      end
!
!
      subroutine convert_to_io_dp(arr_internal,arr_io,n)
!
      use geometry
      use indices
      use compar
      use functions
!
      implicit none
!
      double precision   arr_io(*) , arr_internal(*)
      integer            n,i,j,k,ijk,ijk_io
!
      do k = 1,kmax2
         do j = 1,jmax2
            do i = 1,imax2
               ijk    = funijk(i,j,k)
               ijk_io = funijk_io(i,j,k)
               arr_io(ijk_io) = arr_internal(ijk)
            end do
         end do
      end do
!
      return
      end

!
!
      subroutine convert_from_io_i(arr_io,arr_internal,n)
!
      use geometry
      use indices
      use compar
      use functions
!
      implicit none
!
      integer   arr_io(*) , arr_internal(*)
      integer   n,i,j,k,ijk,ijk_io
!
      do k = 1,kmax2
         do j = 1,jmax2
            do i = 1,imax2
               ijk    = funijk(i,j,k)
               ijk_io = funijk_io(i,j,k)
               arr_internal(ijk) = arr_io(ijk_io)
            end do
         end do
      end do
!
      return
      end
!
!
      subroutine convert_to_io_i(arr_internal,arr_io,n)
!
      use geometry
      use indices
      use compar
      use functions
!
      implicit none
!
      integer   arr_io(*) , arr_internal(*)
      integer   n,i,j,k,ijk,ijk_io
!
      do k = 1,kmax2
         do j = 1,jmax2
            do i = 1,imax2
               ijk    = funijk(i,j,k)
               ijk_io = funijk_io(i,j,k)
               arr_io(ijk_io) = arr_internal(ijk)
            end do
         end do
      end do
!
      return
      end
!
!
      subroutine convert_to_io_c(arr_internal,arr_io,n)
!
      use geometry
      use indices
      use compar
      use functions
!
      implicit none
!
      character(len=3)   arr_io(*) , arr_internal(*)
      integer       n,i,j,k,ijk,ijk_io
!
      do k = 1,kmax2
         do j = 1,jmax2
            do i = 1,imax2
               ijk    = funijk(i,j,k)
               ijk_io = funijk_io(i,j,k)
               arr_io(ijk_io) = arr_internal(ijk)
            end do
         end do
      end do
!
      return
      end
