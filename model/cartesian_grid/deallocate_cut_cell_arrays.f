
      SUBROUTINE Deallocate_CUT_CELL_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: Deallocate_CUT_CELL_ARRAYS
!  Purpose: Deallocate arrays
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      Use indices

      USE cutcell

      IMPLICIT NONE


      Deallocate(  INTERSECT_X  )
      Deallocate(  INTERSECT_Y  )
      Deallocate(  INTERSECT_Z  )

      Deallocate(  X_int  )
      Deallocate(  Y_int  )
      Deallocate(  Z_int  )

      Deallocate(  Xn_int )
      Deallocate(  Xn_U_int )
      Deallocate(  Xn_V_int )
      Deallocate(  Xn_W_int )

      Deallocate(  Ye_int )
      Deallocate(  Ye_U_int )
      Deallocate(  Ye_V_int )
      Deallocate(  Ye_W_int )

      Deallocate(  Zt_int )
      Deallocate(  Zt_U_int )
      Deallocate(  Zt_V_int )
      Deallocate(  Zt_W_int )

      Deallocate( DELX_Ue   )
      Deallocate( DELX_Uw   )
!      Deallocate( DELY_Un   )
!      Deallocate( DELY_Us   )
!      Deallocate( DELZ_Ut   )
!      Deallocate( DELZ_Ub   )

!      Deallocate( DELX_Ve  )
!      Deallocate( DELX_Vw  )
      Deallocate( DELY_Vn  )
      Deallocate( DELY_Vs  )
!      Deallocate( DELZ_Vt  )
!      Deallocate( DELZ_Vb  )

!      Deallocate( DELX_We  )
!      Deallocate( DELX_Ww  )
      Deallocate( DELY_Wn  )
      Deallocate( DELY_Ws  )
      Deallocate( DELZ_Wt  )
      Deallocate( DELZ_Wb  )

      Deallocate( X_U_ec  )
      Deallocate( Y_U_ec  )
      Deallocate( Z_U_ec  )
      Deallocate( X_U_nc  )
      Deallocate( Y_U_nc  )
      Deallocate( Z_U_nc  )
      Deallocate( X_U_tc  )
      Deallocate( Y_U_tc  )
      Deallocate( Z_U_tc  )

      Deallocate( X_V_ec  )
      Deallocate( Y_V_ec  )
      Deallocate( Z_V_ec  )
      Deallocate( X_V_nc  )
      Deallocate( Y_V_nc  )
      Deallocate( Z_V_nc  )
      Deallocate( X_V_tc  )
      Deallocate( Y_V_tc  )
      Deallocate( Z_V_tc  )

      Deallocate( X_W_ec  )
      Deallocate( Y_W_ec  )
      Deallocate( Z_W_ec  )
      Deallocate( X_W_nc  )
      Deallocate( Y_W_nc  )
      Deallocate( Z_W_nc  )
      Deallocate( X_W_tc  )
      Deallocate( Y_W_tc  )
      Deallocate( Z_W_tc  )

      Deallocate(  SNAP )


      RETURN
      END SUBROUTINE Deallocate_CUT_CELL_ARRAYS


