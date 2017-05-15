!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_ROg                                            !
!  Purpose: User hook for calculating the gas phase density.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_ROg(IJK)

! Modules
!---------------------------------------------------------------------//
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use fldvar, only: t_g, x_g, p_g, ro_g
      use param1, only: undefined_i, zero, one, half
      use physprop, only: mw_g, mw_avg, nmax
      use run, only: units
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!
! if using this quantity then remove definition of ier
      ier = 1

! Assign the fluid density
      RO_G(IJK) = ZERO

      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("gas density")')
         CALL INIT_ERR_MSG('USR_PROP_ROg')
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF
      RETURN
      END SUBROUTINE USR_PROP_ROg


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_CPg                                            !
!  Purpose: User hook for calculating the gas phase constant pressure  !
!  specific heat.                                                      !
!                                                                      !
!  Comments:                                                           !
!  - The specific heat assigned in this routine only applies to the    !
!    mixture average specific heat invoked in the general energy       !
!    equations. MFIX has no global variable representing species       !
!    specific heats. The reason for this limitation is explained.      !
!                                                                      !
!    Species specific heat values are locally evaluated based on a     !
!    specific polynominal format set by the Burcat database. Values    !
!    for the polynominal coefficients are read from either the         !
!    database or the mfix.dat. Species specific heats are needed by    !
!    reacting flow simulations.                                        !
!                                                                      !
!    Inconsistencies may arise in reacting flow systems if a user      !
!    specifies a phase average specific heat that is not consistent    !
!    with the values of the species specific heats that comprise that  !
!    phase.                                                            !
!    - IT may be possible to circumvent the matter long as the formula !
!      for the specific heat follows the burcat database form and      !
!      the quantites Thigh(M,N), Tlow(M,N), Tcom, Alow(M,N) and        !
!      Ahigh(M,N) are all appropriately assigned. However, the same    !
!      can be achieved by simply entering the data into the mfix.dat   !
!      as indicated by the user guide.                                 !
!    - To permit different forms of the polynominal for species        !
!      specific heats would require more code development to ensure    !
!      reacting flow simulations evaluate/reference the species        !
!      specific heats accordingly.                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_CPg(IJK)

! Modules
!---------------------------------------------------------------------//
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use constant, only: RGAS => GAS_CONST_cal
      use fldvar, only: t_g, x_g
      use param1, only: undefined_i, zero, one, half
      use physprop, only: mw_g, c_pg, nmax
      use read_thermochemical, only: calc_CpoR
      use run, only: units
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!
! if using this quantity then remove definition of ier
      ier = 1

! Assign the fluid specific heat
      C_PG(IJK) = ZERO

      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("gas specific heat")')
         CALL INIT_ERR_MSG('USR_PROP_CPg')
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF
      RETURN
      END SUBROUTINE USR_PROP_CPg


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_Mug                                            !
!  Purpose: User hook for calculating the gas phase viscosity.         !
!                                                                      !
!  Comments:                                                           !
!  - Assign a value to gas phase viscosity (mu_g). MFIX uses the       !
!    concept of second viscosity (lambda_g) for its internal           !
!    calculations of stress. Second viscosity is defined as follows:   !
!       lambda_g = mu_gbulk - 2/3mu_g                                  !
!       where mu_gbulk is the gas phase bulk viscosity.                !
!    MFIX automatically calculates a value for the gas phase second    !
!    viscosity by assuming the gas phase bulk viscosity is zero. As a  !
!    result, gas phase second  viscosity is evaluated as follows:      !
!       lambda_g = -2/3mu_g                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_Mug(IJK)

! Modules
!---------------------------------------------------------------------//
      use constant, only: to_si
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use fldvar, only: t_g, x_g, p_g, ro_g
      use param1, only: undefined_i, zero, one, half
      use physprop, only: mu_g, mw_g, mw_avg, nmax
      use run, only: units
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!
! if using this quantity then remove definition of ier
      ier = 1

! Assign the fluid viscosity
      Mu_g(IJK) = ZERO

      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("gas viscosity")')
         CALL INIT_ERR_MSG('USR_PROP_MUG')
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF
      RETURN
      END SUBROUTINE USR_PROP_Mug


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_Kg                                             !
!  Purpose: User hook for calculating the gas phase conductivity.      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_Kg(IJK)

! Modules
!---------------------------------------------------------------------//
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use fldvar, only: t_g, x_g, p_g, ro_g
      use param1, only: undefined_i, zero, one, half
      use physprop, only: k_g, mw_g, mw_avg, nmax
      use run, only: units
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!
! if using this quantity then remove definition of ier
      ier = 1

! Assign the fluid conductivity
      K_g(IJK) = ZERO

      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("gas conductivity")')
         CALL INIT_ERR_MSG('USR_PROP_KG')
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF
      RETURN
      END SUBROUTINE USR_PROP_Kg


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_Difg                                           !
!  Purpose: User hook for calculating the diffusivity of the gas phase !
!  species.                                                            !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_Difg(IJK,N)

! Modules
!---------------------------------------------------------------------//
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use fldvar, only: t_g, x_g, p_g, ro_g, rop_g
      use param1, only: undefined_i, zero, one, half
      use physprop, only: dif_g, mw_g, mw_avg, nmax
      use run, only: units
      use scales, only: unscale_pressure
      use toleranc, only: zero_x_gs
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK
! species index
      INTEGER, INTENT(IN) :: N

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!
! if using this quantity then remove definition of ier
      ier = 1

! Assign the fluid phase species diffusivity
      Dif_g(IJK,N) = ZERO

      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("gas phase species ",I2," diffusivity")') N
         CALL INIT_ERR_MSG('USR_PROP_DIFG')
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF
      RETURN
      END SUBROUTINE USR_PROP_Difg


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_ROs                                            !
!  Purpose: User hook for calculating solids phase density.            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_ROs(IJK,M)

! Modules
!---------------------------------------------------------------------//
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use fldvar, only: ro_s, T_s, X_s
      use param1, only: undefined_i, zero, one, half
      use physprop, only: mw_s, nmax
      use run, only: units
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!
! if using this quantity then remove definition of ier
      ier = 1

! Assign the solids phase density
      RO_S(IJK,M) = ZERO

      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("solids phase ",I2," density")') M
         CALL INIT_ERR_MSG('USR_PROP_ROS')
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF
      RETURN
      END SUBROUTINE USR_PROP_ROs


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_CPs                                            !
!  Purpose: User hook for calculating solids phase constant pressure   !
!  specific heat.                                                      !
!                                                                      !
!  Comments:                                                           !
!  - See comments under USER_PROP_CPg                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_CPs(IJK, M)
! Modules
!---------------------------------------------------------------------//
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use constant, only: RGAS => GAS_CONST_cal
      use fldvar, only: t_s, x_s
      use param1, only: undefined_i, zero, one, half
      use physprop, only: mw_s, c_ps, nmax
      use read_thermochemical, only: calc_CpoR
      use run, only: units
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!
! if using this quantity then remove definition of ier
      ier = 1

! Assign the solids phase specific heat
      C_PS(IJK,M) = ZERO

      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("solids phase ",I2," specific heat")') M
         CALL INIT_ERR_MSG('USR_PROP_CPs')
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF
      RETURN
      END SUBROUTINE USR_PROP_CPs


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_Mus                                            !
!  Purpose: User hook for calculating the solids phase viscosity.      !
!                                                                      !
!  Comments:                                                           !
!  - Assign a value to solids phase viscosity (mu_g), second viscosity !
!    (lambda_s) and solids pressure (p_s). MFIX uses the concept of    !
!    second viscosity for its internal calculations of stress. Second  !
!    viscosity is defined as follows:                                  !
!       lambda_s = mu_sbulk - 2/3 mu_s                                 !
!       where mu_sbulk is the solids phase bulk viscosity.             !
!                                                                      !
!  - In assigning a value to solids pressure consider how MFIX's       !
!    governing equations have been posed and what physics at are       !
!    involved. Nominally the governing equations are written assuming  !
!    a gas-solids system wherein the gradient in gas phase pressure    !
!    is present in the solids phase momentum balances.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_Mus(IJK,M)

! Modules
!---------------------------------------------------------------------//
      use constant, only: to_si
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use fldvar, only: t_s, x_s, rop_s, ro_s, ep_s, p_s
      use param1, only: undefined_i, zero, one, half
      use physprop, only: mw_s, nmax
      use run, only: units
      use visc_s, only: mu_s, lambda_s
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!
! if using this quantity then remove definition of ier
      ier = 1

! Assign the solids phase viscosity, second viscosity and pressure
      Mu_s(IJK,M) = ZERO
      LAMBDA_S(IJK,M) = ZERO
      P_S(IJK,M) = ZERO

      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("solids phase ",I2," viscosity")') M
         CALL INIT_ERR_MSG('USR_PROP_Mus')
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF
      RETURN
      END SUBROUTINE USR_PROP_Mus


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_Ks                                             !
!  Purpose: User hook for calculating solids phase conductivity.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_Ks(IJK,M)

! Modules
!---------------------------------------------------------------------//
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use fldvar, only: ro_s, T_s, X_s, ep_s, ep_g, rop_s
      use param1, only: undefined_i, zero, one, half
      use physprop, only: mw_s, nmax, k_s, k_g
      use run, only: units
      use toleranc, only: dil_ep_s
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!
! if using this quantity then remove definition of ier
      ier = 1

! Assign the solids phase conductivity
      K_S(IJK,M) = ZERO

      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("solids phase ",I2," conductivity")') M
         CALL INIT_ERR_MSG('USR_PROP_KS')
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF
      RETURN
      END SUBROUTINE USR_PROP_Ks


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_Difs                                           !
!  Purpose: User hook for calculating diffusivity of 'solids' phase    !
!  species.                                                            !
!                                                                      !
!  Comments:                                                           !
!  - As always consider whether such a term is meaningful in the       !
!    system being modeled.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_Difs(IJK,M,N)

! Modules
!---------------------------------------------------------------------//
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use fldvar, only: ro_s, T_s, X_s, ep_s, ep_g, rop_s
      use param1, only: undefined_i, zero, one, half
      use physprop, only: mw_s, nmax, dif_s
      use run, only: units
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK
! solids phase index
      INTEGER, INTENT(IN) :: M
! species index
      INTEGER, INTENT(IN) :: N

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!
! if using this quantity then remove definition of ier
      ier = 1

! Assign the solids phase species diffusivity
      DIF_S(IJK,M,N) = ZERO

      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("solids phase ",I2," species ", I2, &
                           " diffusivity")') M, N
         CALL INIT_ERR_MSG('USR_PROP_DIFS')
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF
      RETURN
      END SUBROUTINE USR_PROP_Difs


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_Gama                                           !
!  Purpose: User hook for calculating the gas-solids heat transfer     !
!  coefficient due to relative temperature differences between the     !
!  gas phase (M=0) and each 'solids' phase (M=1 to MMAX)               !
!                                                                      !
!  Comments:                                                           !
!  - No solids-solids heat transfer is allowed at this time. To        !
!    account for this term would require appropriate closure and       !
!    additional code development.                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_Gama(IJK,M)

! Modules
!---------------------------------------------------------------------//
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use energy, only: gama_gs
      use fldvar, only: u_g, v_g, w_g, u_s, v_s, w_s
      use fldvar, only: ep_s, ep_g, rop_s, rop_g, ro_s, ro_g
      use fldvar, only: T_s, T_g, d_p
      Use fun_avg, only: avg_x_e, avg_y_n, avg_z_t
      use functions, only: im_of, jm_of, km_of
      use indices, only: i_of
      use param1, only: undefined_i, zero, one, half
      use param1, only: small_number, large_number
      use physprop, only: mu_g, k_g, c_pg
      use run, only: units
      use rxns, only: r_phase
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=50) :: err_prop
!......................................................................!
! if using this quantity then remove definition of ier
      ier = 1

! Assign the gas-solids phase heat transfer coefficient
      GAMA_GS(IJK,M) = ZERO

      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("gas-solids phase ",I2," heat transfer ", &
                           "coefficient")') M
         CALL INIT_ERR_MSG('USR_PROP_Gama')
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF
      RETURN
      END SUBROUTINE USR_PROP_GAMA


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_FSS                                            !
!  Purpose: User hook for calculating the solids-solids drag           !
!  coefficient due to relative velocity differences between each of    !
!  solids phases (M = 1 to MMAX).                                      !
!                                                                      !
!  Comments:                                                           !
!  - solids-solids drag is momentum transfer due to relative velocity  !
!    differences between solids phases M and L, where M and L range    !
!    from 1 to MMAX and M!=L.                                          !
!  - this implementation is currently restricted to solids-solids      !
!    drag between continuum solids phases. Further development is      !
!    needed to interact with discrete phases.                          !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_Fss(IJK,L,M)

! Modules
!---------------------------------------------------------------------//
      use constant, only: segregation_slope_coefficient
      use drag, only: f_ss
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use fldvar, only: u_s, v_s, w_s
      use fldvar, only: ep_s, ep_g, rop_s, ro_s
      use fldvar, only: d_p, p_star
      use fun_avg, only: avg_x_e, avg_y_n, avg_z_t
      use functions, only: funlm
      use functions, only: im_of, jm_of, km_of
      use indices, only: i_of
      use param1, only: undefined_i, zero, one, half
      use param1, only: small_number, large_number
      use physprop, only: smax, close_packed
      use run, only: units
      use rdf, only: g_0
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK
! solids phase indices
      INTEGER, INTENT(IN) :: M, L

! Local Variables:
!---------------------------------------------------------------------//
! index for storing solids-solids drag coefficients in the upper
! triangle of the matrix
      INTEGER :: LM
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=60) :: err_prop
!......................................................................!
! if using this quantity then remove definition of ier
      ier = 1

! Assign the solids-solids drag coefficient
      LM = FUNLM(L,M)
      F_SS(IJK,LM) = ZERO

      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("solids phase ",I2," solids phase ", I2, &
                           " drag coefficient")') M, L
         CALL INIT_ERR_MSG('USR_PROP_FSS')
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF
      RETURN
      END SUBROUTINE USR_PROP_FSS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: USR_PROPERTIES                                          C
!  Purpose: Hook for user defined physical properties and transport    C
!  coefficients including interphase exchange coefficients for all     C
!  phases. The exact quantities are listed here for clarity.           C
!                                                                      C
!         Quantity      GAS (M=0)     SOLIDS (M=1toMMAX)               C
!         density       RO_G(:)        RO_S(:,M)                       C
!      specific heat    C_PG(:)        C_PS(:,M)                       C
!       conductivity    K_G(:)         K_S(:,M)                        C
!       diffusivity     DIF_G(:,N)     DIF_S(:,M,N)                    C
!        viscosity      MU_G(:)        MU_S(:,M)                       C
!     gas-solids drag                  F_GS(:,M)                       C
!   solids-solids drag                 F_SS(:,LM)                      C
!   gas-solids heat tr.                GAMA(:,M)                       C
!                                                                      C
!  Comments:                                                           C
!  - gas-solids drag is momentum transfer due to relative velocity     C
!    differences (skin friction and form drag) between the gas phase   C
!    (M=0) and each solids phase (M=1 to MMAX).                        C
!  - solids-solids drag is momentum transfer due to relative velocity  C
!    differences between solids phases M and L, where M and L range    C
!    from 1 to MMAX and M!=L.                                          C
!  - gas-solids heat transfer is heat transfer due to relative         C
!    temperature differences between the gas phase phase (M=0) and     C
!    each solids phase (M=1 to MMAX).                                  C
!  - No solids-solids heat transfer is allowed. To account for such    C
!    would require appropriate closures and additional code            C
!    development.                                                      C
!  - the specific heat assigned in this routine only applies to the    C
!    mixture average specific heat invoked in the general energy       C
!    equations. reacting flow simulations require values for species   C
!    specific heat and are taken from the Burcat database or read      C
!    in that format from the mfix.dat file. therefore inconsitencies   C
!    may arise in calculations involving species specific heat. this   C
!    requires further development to fully address.                    C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR_PROPERTIES(lprop, IJK, M, N)

! Modules
!-----------------------------------------------
      use constant, only: pi, gas_const, gravity
      use error_manager

      use fldvar, only: u_g, v_g, w_g, ep_g
      use fldvar, only: u_s, v_s, w_s, ep_s
      use fldvar, only: p_g, rop_g, ro_g, T_g, X_g
      use fldvar, only: p_s, rop_s, ro_s, T_s, X_s
      use fldvar, only: d_p, theta_m

      use fldvar, only: scalar
      use functions
      use geometry

      use indices, only: i_of, j_of, k_of
      use indices, only: im1, ip1, jm1, jp1, km1, kp1

      use param1, only: zero, one, half, undefined, undefined_i
      use physprop, only: k_g, c_pg, dif_g, mu_g
      use physprop, only: k_s, c_ps, dif_s
      use scalars, only: phase4scalar
      use visc_s, only: mu_s, lambda_s
      use usr_prop
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------
! reference equation
      INTEGER, INTENT(IN) :: lprop
! index
      INTEGER, INTENT(IN) :: IJK
! Phase index
      INTEGER, INTENT(IN) :: M
! Species index or second solids phase index for solids-solids drag
! (if applicable otherwise undefined_i)
      INTEGER, INTENT(IN) :: N

! Local variables
!-----------------------------------------------
! Error flag
      INTEGER :: IER
      CHARACTER(len=40):: err_prop

!-----------------------------------------------
! initialize
      ier = undefined_i

! in each case the ier flag is set to ensure that if a user defined
! quantity is invoked that the associated user defined quantity is
! actually specified. this flag must be removed by the user or the
! code will fail.

      SELECT CASE(lprop)

!
      CASE (GAS_DENSITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("gas density")')
!        RO_G(IJK) =


      CASE (GAS_SPECIFICHEAT)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("gas specific heat")')
!        C_PG(IJK,M) =


! assign gas viscosity value mu_g. bulk viscosity is taken as zero.
! second viscosity (lambda_g) is automatically defined as
! lambda_g = -2/3mu_g
      CASE (GAS_VISCOSITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("gas viscosity")')
!        MU_G(IJK) =


      CASE (GAS_CONDUCTIVITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("gas conductivity")')
!        K_G(IJK) =


      CASE (GAS_DIFFUSIVITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("gas diffusivity")')
!        DIF_G(IJK,N) =


      CASE (SOLIDS_DENSITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("solids phase ",I2," density")') M
!        RO_S(IJK,M) =


      CASE (SOLIDS_SPECIFICHEAT)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("solids phase ",I2," specific heat")') M
!        C_PS(IJK,M) =


! assign solids phase M viscosity mu_s.
! assign second viscosity (lambda_s): lambda_s = mu_sbulk - 2/3mu_s
! assign solids pressure p_s.
      CASE (SOLIDS_VISCOSITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("solids phase ", I2," viscosity")') M
!        MU_S(IJK,M) =
!        LAMBDA_S(IJK,M) =
!        P_S(IJK,M) =


      CASE (SOLIDS_CONDUCTIVITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("solids phase ",I2," conductivity")') M
!        K_S(IJK,M) =


      CASE (SOLIDS_DIFFUSIVITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("solids phase",I2," diffusivity")') M
!        DIF_S(IJK,M,N) =


      CASE (GASSOLIDS_HEATTRANSFER)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, &
               '("gas-solids phase ",I2," heat transfer coefficient")') M
!       GAMA_GS(IJK,M) =


      CASE (GASSOLIDS_DRAG)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, &
               '("gas-solids phase ",I2," drag coefficient")') M
!        F_GS(IJK,M) =


      CASE (SOLIDSSOLIDS_DRAG)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, &
            '("solids-solids phases ",I2," & ",I2," drag coefficient")') M, N
!        F_SS(IJK,LM) =


      END SELECT

      IF (IER /= UNDEFINED_I) THEN
         CALL INIT_ERR_MSG('USR_PROPERTIES')
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exits.',&
         'Either choose a different model or correct',/,'mfix/model/',&
         'usr_properties.f')
      ENDIF

      RETURN
      END SUBROUTINE USR_PROPERTIES

