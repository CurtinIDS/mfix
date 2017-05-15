!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: REPORT_STATS_PIC                                        !
!                                                                      !
!  Purpose: Output stats about PIC simulation.                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE REPORT_STATS_PIC

! Global Variables:
!---------------------------------------------------------------------//
! Flag to report minimum EP_G
      use mfix_pic, only: PIC_REPORT_MIN_EPG
! Gas phase volume fraction
      use fldvar, only: EP_G
! Location of cell faces (East, North, Top)
      use discretelement, only: XE, YN, ZT

      use param1, only: large_number
      use mpi_utility
      USE error_manager
      USE functions

      IMPLICIT NONE

! Local Variables:
!----------------------------------------------------------------------!
! Loop counters
      INTEGER I, J, K, IJK, IPROC

      INTEGER :: EPg_MIN_loc(0:numpes-1, 4), EPg_MIN_loc2(1)
      DOUBLE PRECISION :: EPg_MIN(0:numpes-1), EPg_min2

!-----------------------------------------------

      CALL INIT_ERR_MSG("REPORT_STATS_PIC")


      IF(PIC_REPORT_MIN_EPG) THEN

         EPG_MIN(:) = 0
         EPG_MIN(mype) = LARGE_NUMBER

         EPG_MIN_LOC(:,:) = 0
         EPG_MIN_LOC(mype,:) = -1

         DO K = KSTART1, KEND1
            DO J = JSTART1, JEND1
               DO I = ISTART1, IEND1
                  IJK = funijk(I,J,K)

                  IF(EP_G(IJK) < EPG_MIN(mype)) THEN
                     EPG_MIN_LOC(mype,1) = I
                     EPG_MIN_LOC(mype,2) = J
                     EPG_MIN_LOC(mype,3) = K
                     EPG_MIN_LOC(mype,4) = IJK
                     EPG_MIN(mype) = EP_G(IJK)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         call GLOBAL_ALL_SUM(EPg_MIN)
         CALL GLOBAL_ALL_SUM(EPg_MIN_loc)

         epg_min2     = MINVAL(epg_min(0:numpes-1))
         epg_min_loc2 = MINLOC(epg_min(0:numpes-1)) - 1
         !-1, since minloc goes from 1:size of the array.
         !If not corrected by -1, then the proc id will be off by 1

         iproc = epg_min_loc2(1)

         I     = epg_min_loc(iproc, 1)
         J     = epg_min_loc(iproc, 2)
         K     = epg_min_loc(iproc, 3)
         IJK   = epg_min_loc(iproc, 4)
         WRITE(ERR_MSG,1014) EPG_MIN2, Iproc, I, J, K, IJK, &
            XE(I) - 0.5*DX(I), YN(J)-0.5*DY(J), ZT(K) - 0.5*DZ(K)

 1014       FORMAT( /, &
            &      5x,'EPGMIN                    = ', 2x,g17.8,/ &
            &      5x,'EPGMIN PROC RANK          = ', 2x, I10, / &
            &      5x,'EPGMIN (I, J, K, IJK)     = ', 3(2x,i5),2x,i10,/ &
            &      5x,'XMID, YMID, ZMID FOR CELL = ', 3(2x,g17.8))

            call flush_err_msg(header = .false., footer = .false.)

      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE REPORT_STATS_PIC


