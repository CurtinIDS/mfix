      PROGRAM POST_MFIX

      INCLUDE 'xforms.inc'

      DO_XFORMS = .FALSE.

      CALL F_INIT

      STOP
    END PROGRAM POST_MFIX

      SUBROUTINE  CHECK_INTER(inter)
        integer          :: inter
        RETURN
      END SUBROUTINE CHECK_INTER

      SUBROUTINE       ADD_TO_RESULTS_BROWSER(line)
        character(len=*) :: line
        RETURN
      END SUBROUTINE ADD_TO_RESULTS_BROWSER

      SUBROUTINE       ADD_TO_SPX_BR(sel,line)
        character(len=*) :: line
        logical          :: sel
        RETURN
      END SUBROUTINE ADD_TO_SPX_BR

      SUBROUTINE       SPX_TIME_SELECTED(i1,i2)
        integer          :: i1 , i2
        RETURN
      END SUBROUTINE SPX_TIME_SELECTED

      SUBROUTINE       SPX_DESELECT_TIME(i1)
        integer          :: i1
        RETURN
      END SUBROUTINE SPX_DESELECT_TIME

      SUBROUTINE       GET_PTX_G
        RETURN
      END SUBROUTINE GET_PTX_G

