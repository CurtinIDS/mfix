!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: INIT_NAMELIST                                           !
!  Purpose: initialize the NAMELIST variables                          !
!                                                                      !
!  Author: P. Nicoletti                               Date: 26-NOV-91  !
!                                                                      !
!  Keyword Documentation Format:                                       !
!                                                                      !
!<keyword category="category name" required="true"/FALSE               !
!                                    legacy=TRUE/FALSE>                !
!  <description></description>                                         !
!  <arg index="" id="" max="" min=""/>                                 !
!  <dependent keyword="" value="DEFINED"/>                             !
!  <conflict keyword="" value="DEFINED"/>                              !
!  <valid value="" note="" alias=""/>                                  !
!  <range min="" max="" />                                             !
!  MFIX_KEYWORD=INIT_VALUE                                             !
!</keyword>                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE INIT_NAMELIST

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE cdist
      USE cg_init_namelist, ONLY: CARTESIAN_GRID_INIT_NAMELIST
      USE compar
      USE constant
      USE fldvar
      USE geometry
      USE ic
      USE indices
      USE is
      USE iterate, only: max_nit
      USE leqsol
      USE output
      USE parallel
      USE param
      USE param1
      USE physprop
      USE ps
      USE residual
      USE run
      USE rxns
      USE scalars
      USE scales
      USE stiff_chem
      USE toleranc
      USE ur_facs
      use usr_src, only: call_usr_source
! user defined flags
      use usr_prop, only: usr_rog, usr_cpg, usr_kg, usr_mug, usr_difg
      use usr_prop, only: usr_ros, usr_cps, usr_ks, usr_mus, usr_difs
      use usr_prop, only: usr_gama, usr_fgs, usr_fss

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop counters
      INTEGER :: LC

!#####################################################################!
!                             Run Control                             !
!#####################################################################!

!<keyword category="Run Control" required="true">
!  <description> Name used to create output files. The name should
!    generate legal file names after appending extensions.
!    Ex: Given the input, RUN_NAME = "bub01", MFIX will generate
!    the output files: BUB01.LOG, BUB01.OUT, BUB01.RES, etc.
!  </description>
      RUN_NAME = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Problem description. Limited to 60 characters.</description>
      DESCRIPTION = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="true">
!  <description> Simulation input/output units.</description>
!  <valid value="cgs" note="All input and output in CGS units (g, cm, s, cal)."/>
!  <valid value="si" note="All input and output in SI units (kg, m, s, J)."/>
      UNITS = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="true">
!  <description>Type of run.</description>
!  <valid value="new" note="A new run. There should be no .RES, .SPx,
!    .OUT, or .LOG files in the run directory."/>
!  <valid value="RESTART_1" note="Traditional restart. The run continues
!    from the last time the .RES file was updated and new data is added
!    to the SPx files."/>
!  <valid value="RESTART_2"
!    note="Start a new run with initial conditions from a .RES file
!      created from another run. No other data files (SPx) should be
!      in the run directory."/>
      RUN_TYPE = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Simulation start time. This is typically zero.
!  </description>
!  <range min="0.0" max="+Inf" />
      TIME = UNDEFINED
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Simulation stop time.
!  </description>
!  <range min="0.0" max="+Inf" />
      TSTOP = UNDEFINED
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Initial time step size. If left undefined, a steady-state
!    calculation is performed.
!  </description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT = UNDEFINED
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Maximum time step size.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT_MAX = ONE
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Minimum time step size.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT_MIN = 1.0D-6
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Factor for adjusting time step.
!    * The value must be less than or equal to 1.0.
!    * A value of 1.0 keeps the time step constant which may help overcome
!      initial non-convergence.
!  </description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="1" />
      DT_FAC = 0.9D0
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Force a forward time-step if the maximum number of iterations,
!    MAX_NIT, is reached. The forward time-step is only forced after
!    reaching the minimum time-step, DT_MIN, for adjustable time-step
!    simulations (DT_FAC /= 1). This option should be used with caution
!    as unconverged time-steps may lead to poor simulation results and/or
!    additional convergence issues.
!  </description>
!  <valid value=".TRUE." note="Force forward time-step when DT=DT_MIN and
!    the maximum number of iterations are reached."/>
!  <valid value=".FALSE." note="Abort run when DT < DT_MIN."/>
      PERSISTENT_MODE = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Flag to restart the code when DT < DT_MIN.
!  </description>
      AUTO_RESTART = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Flag to enable/disable solving the X-momentum equations.
!  </description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".TRUE." note="Solve X-momentum equations."/>
!  <valid value=".FALSE." note="The X velocity initial conditions
!   persist throughout the simulation."/>
      MOMENTUM_X_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Flag to enable/disable solving the Y-momentum equations.
! </description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".TRUE." note="Solve Y-momentum equations."/>
!  <valid value=".FALSE." note="The Y velocity initial conditions
!   persist throughout the simulation."/>
      MOMENTUM_Y_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Flag to enable/disable solving the Z-momentum equations.
!  </description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".TRUE." note="Solve Z-momentum equations."/>
!  <valid value=".FALSE." note="The Z velocity initial conditions
!   persist throughout the simulation."/>
      MOMENTUM_Z_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Flag to enable Jackson form of momentum equations.
!    See Anderson and Jackson, (1967), IECF, 6(4), p.527.
!  </description>
!  <valid value=".TRUE." note="Solve Jackson form of momentum equations."/>
!  <valid value=".FALSE." note="Default form."/>
      JACKSON = .FALSE.
!</keyword>
!<keyword category="Run Control" required="false">
!  <description>
!    Flag to enable Ishii form of momentum equations.
!    See Ishii, (1975), Thermo-fluid dynamic theory of two-phase flow.
!  </description>
!  <valid value=".TRUE." note="Solve Ishii form of momentum equations."/>
!  <valid value=".FALSE." note="Default form."/>
      ISHII = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Solve energy equations.</description>
!  <valid value=".TRUE." note="Solve energy equations."/>
!  <valid value=".FALSE." note="Do not solve energy equations."/>
      ENERGY_EQ = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Solve species transport equations.</description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".TRUE." note="Solve species equations."/>
!  <valid value=".FALSE." note="Do not solve species equations."/>
      SPECIES_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false" tfm="true">
!  <description>Granular energy formulation selection.</description>
!  <valid value=".FALSE."
!    note="Use algebraic granular energy equation formulation."/>
!  <valid value=".TRUE."
!    note="Use granular energy transport equation (PDE) formulation."/>
      GRANULAR_ENERGY = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    The K-Epsilon turbulence model (for single-phase flow).
!    o Numerical parameters (like under-relaxation) are the same as the
!      ones for SCALAR (index = 9).
!    o All walls must be defined (NSW, FSW or PSW) in order to use
!      standard wall functions. If a user does not specify a wall type,
!      the simulation will not contain the typical turbulent profile in
!      wall-bounded flows.
!  </description>
!  <dependent keyword="MU_GMAX" value="DEFINED"/>
!  <conflict keyword="L_SCALE0" value="DEFINED"/>
!  <valid value=".TRUE."  note="Enable the K-epsilon turbulence model
!    (for single-phase flow) using standard wall functions."/>
!  <valid value=".FALSE." note="Do not use K-epsilon turbulence model"/>
      K_EPSILON = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Value of turbulent length initialized. This may be overwritten
!    in specific regions with the keyword IC_L_SCALE.
!</description>
!  <dependent keyword="MU_GMAX" value="DEFINED"/>
!  <conflict keyword="K_EPSILON" value=".TRUE."/>
      L_SCALE0 = ZERO
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Maximum value of the turbulent viscosity of the fluid, which
!    must be defined if any turbulence model is used.
!    A value MU_GMAX =1.E+03 is recommended. (see calc_mu_g.f)
!  </description>
      MU_GMAX = UNDEFINED
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!     Available gas-solids drag models.
!     Note: The extension _PCF following the specified drag model
!     indicates that the polydisperse correction factor is available.
!     For PCF details see:
!     o Van der Hoef MA, Beetstra R, Kuipers JAM. (2005)
!       Journal of Fluid Mechanics.528:233-254.
!     o Beetstra, R., van der Hoef, M. A., Kuipers, J.A.M. (2007).
!       AIChE Journal, 53:489-501.
!     o Erratum (2007), AIChE Journal, Volume 53:3020
!  </description>
!
!  <valid value="SYAM_OBRIEN" note="Syamlal M, OBrien TJ (1988).
!   International Journal of Multiphase Flow 14:473-481.
!   Two additional parameters may be specified: DRAG_C1, DRAG_D1"/>
!
!  <valid value="GIDASPOW" note="Ding J, Gidaspow D (1990).
!   AIChE Journal 36:523-538"/>
!
!  <valid value="GIDASPOW_BLEND" note="Lathouwers D, Bellan J (2000).
!    Proceedings of the 2000 U.S. DOE
!        Hydrogen Program Review NREL/CP-570-28890."/>
!
!  <valid value="WEN_YU" note="Wen CY, Yu YH (1966).
!   Chemical Engineering Progress Symposium Series 62:100-111."/>
!
!  <valid value="KOCH_HILL" note="Hill RJ, Koch DL, Ladd JC (2001).
!   Journal of Fluid Mechanics, 448: 213-241. and 448:243-278."/>
!
!  <valid value="BVK" note="Beetstra, van der Hoef, Kuipers (2007).
!   Chemical Engineering Science 62:246-255"/>
!
!  <valid value="HYS" note="Yin, X, Sundaresan, S. (2009).
!   AIChE Journal 55:1352-1368
!   This model has a lubrication cutoff distance, LAM_HYS, that can be
!   specified."/>
!
!  <valid value="USER_DRAG" note="Invoke user-defined drag law. (usr_drag.f)"/>
!
!  <valid value="GIDASPOW_PCF" note="see GIDASPOW"/>
!  <valid value="GIDASPOW_BLEND_PCF" note="see GIDASPOW_BLEND"/>
!  <valid value="WEN_YU_PCF" note="see WEN_YU"/>
!  <valid value="KOCH_HILL_PCF" note="see KOCH_HILL"/>
!
      DRAG_TYPE = 'SYAM_OBRIEN'
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Quantity for calibrating Syamlal-O'Brien drag correlation using Umf
!    data.  This is determined using the Umf spreadsheet.
!  </description>
      DRAG_C1 = 0.8d0
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Quantity for calibrating Syamlal-O'Brien drag correlation using Umf
!    data.  This is determined using the Umf spreadsheet.
!  </description>
      DRAG_D1 = 2.65d0
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    The lubrication cutoff distance for HYS drag model.  In practice
!    this number should be on the order of the mean free path of the
!    gas for smooth particles, or the RMS roughness of a particle if
!    they are rough (if particle roughness is larger than the mean
!   free path).
!  </description>
!  <dependent keyword="DRAG_TYPE" value="HYS"/>
      LAM_HYS = UNDEFINED
!</keyword>

!<keyword category="Run Control" required="false" tfm="true">
!  <description>
!    Subgrid models.
!  </description>
!
!  <valid value="Igci" note="
!   Igci, Y., Pannala, S., Benyahia, S., and Sundaresan S. (2012).
!   Industrial & Engineering Chemistry Research, 2012, 51(4):2094-2103"/>
!
!  <valid value="Milioli" note="
!   Milioli, C.C., Milioli, F. E., Holloway, W., Agrawal, K. and
!   Sundaresan, S. (2013). AIChE Journal, 59:3265-3275."/>
!
      SUBGRID_TYPE = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="false" tfm="true">
!  <description>
!    Ratio of filter size to computational cell size.
!  </description>
      FILTER_SIZE_RATIO = 2.0D0
!</keyword>

!<keyword category="Run Control" required="false" tfm="true">
!  <description>Flag for subgrid wall correction.</description>
!  <valid value=".FALSE." note="Do not include wall correction."/>
!  <valid value=".TRUE." note="Include subgrid wall correction."/>
      SUBGRID_Wall = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Shared gas-pressure formulation. See Syamlal, M. and Pannala, S.
!    "Multiphase continuum formulation for gas-solids reacting flows,"
!    chapter in Computational Gas-Solids Flows and Reacting Systems:
!    Theory, Methods and Practice, S. Pannala, M. Syamlal and T.J.
!    O'Brien (editors), IGI Global, Hershey, PA, 2011.
!  </description>
!  <valid value=".FALSE." note="Use Model A"/>
!  <valid value=".TRUE."  note="Use Model B. Bouillard, J.X.,
!    Lyczkowski, R.W., Folga, S., Gidaspow, D., Berry, G.F. (1989).
!    Canadian Journal of Chemical Engineering 67:218-229."/>
      MODEL_B = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description> The number of user-defined scalar transport equations
!    to solve.
!  </description>
!  <range min="0" max="DIM_SCALAR" />
      NScalar = 0
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    The phase convecting the indexed scalar transport equation.
!  </description>
!  <arg index="1" id="Scalar Equation" min="0" max="DIM_SCALAR"/>
!  <range min="0" max="DIM_M" />
      Phase4Scalar(:DIM_SCALAR) = UNDEFINED_I
!</keyword>

!#####################################################################!
!                           Physical Parameters                       !
!#####################################################################!


!<keyword category="Physical Parameters" required="false">
!  <description>Reference pressure. [0.0]</description>
      P_REF = ZERO
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>Scale factor for pressure. [1.0]</description>
      P_SCALE = ONE
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>Gravitational acceleration. [980.7 in CGS]</description>
      GRAVITY = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    X-component of gravitational acceleration vector. By default, the
!    gravity force acts in the negative y-direction.
!  </description>
      GRAVITY_X = ZERO
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    Y-component of gravitational acceleration vector. By default, the
!    gravity force acts in the negative y-direction.
!  </description>
      GRAVITY_Y = ZERO
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    Z-component of gravitational acceleration vector. By default, the
!    gravity force acts in the negative y-direction.
!  </description>
      GRAVITY_Z = ZERO
!</keyword>





!#####################################################################!
!                          Numerical Parameters                       !
!#####################################################################!



!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Maximum number of iterations [500].
!  </description>
      MAX_NIT = 500
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Factor to normalize the gas continuity equation residual. The
!    residual from the first iteration is used if NORM_G is left
!    undefined. NORM_G=0 invokes a normalization method based on the
!    dominant term in the continuity equation. This setting may speed up
!    calculations, especially near a steady state and incompressible
!    fluids. But, the number of iterations for the gas phase pressure
!    should be increased, LEQ_IT(1), to ensure mass balance
!  </description>
      NORM_G = UNDEFINED
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Factor to normalize the solids continuity equation residual. The
!    residual from the first iteration is used if NORM_S is left
!    undefined. NORM_S = 0 invokes a normalization method based on the
!    dominant term in the continuity equation. This setting may speed up
!    calculations, especially near a steady state and incompressible
!    fluids. But, the number of iterations for the solids volume
!    fraction should be increased, LEQ_IT(2), to ensure mass balance.
!  </description>
      NORM_S = UNDEFINED
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Maximum residual at convergence (Continuity + Momentum) [1.0d-3].
!  </description>
      TOL_RESID = 1.0D-3
!</keyword>


!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Maximum residual at convergence (Energy) [1.0d-4].
!  </description>
      TOL_RESID_T = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Maximum residual at convergence (Species Balance) [1.0d-4].
!  </description>
      TOL_RESID_X = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Maximum residual at convergence (Granular Energy) [1.0d-4].
!  </description>
      TOL_RESID_Th = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Maximum residual at convergence (Scalar Equations) [1.0d-4].
!  </description>
      TOL_RESID_Scalar = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Maximum residual at convergence (K_Epsilon Model) [1.0d-4].
!  </description>
      TOL_RESID_K_Epsilon = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Minimum residual for declaring divergence [1.0d+4].
!    This parameter is useful for incompressible fluid simulations
!    because velocity residuals can take large values for the second
!    iteration (e.g., 1e+8) before dropping down to smaller values for
!    the third iteration.
!  </description>
      TOL_DIVERGE = 1.0D+4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Reduce the time step if the residuals stop decreasing. Disabling this
!    feature may help overcome initial non-convergence.
!  </description>
!  <valid value=".FALSE." note="Continue iterating if residuals stall."/>
!  <valid value=".TRUE."  note="Reduce time step if residuals stall."/>
      DETECT_STALL = .TRUE.
!</keyword>


!<keyword category="Numerical Parameters" required="false">
!  <description>
!    LEQ Solver selection. BiCGSTAB is the default method for all
!    equation types.
!  </description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
!  <valid value="1" note="SOR - Successive over-relaxation"/>
!  <valid value="2" note="BiCGSTAB - Biconjugate gradient stabilized."/>
!  <valid value="3" note="GMRES - Generalized minimal residual method"/>
!  <valid value="5" note="CG - Conjugate gradient"/>
      LEQ_METHOD(:) = 2
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Linear Equation tolerance [1.0d-4].
!  </description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
!  <dependent keyword="LEQ_METHOD" value="2"/>
!  <dependent keyword="LEQ_METHOD" value="3"/>
      LEQ_TOL(:) = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Number of iterations in the linear equation solver.
!    o 20 iterations for equation types 1-2
!    o  5 iterations for equation types 3-5,10
!    o 15 iterations for equation types 6-9
!  </description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
      LEQ_IT(1) =  20
      LEQ_IT(2) =  20
      LEQ_IT(3) =   5
      LEQ_IT(4) =   5
      LEQ_IT(5) =   5
      LEQ_IT(6) =  15
      LEQ_IT(7) =  15
      LEQ_IT(8) =  15
      LEQ_IT(9) =  15
      LEQ_IT(10) =  5
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Linear equation sweep direction. This applies when using GMRES or
!    when using the LINE preconditioner with BiCGSTAB or CG methods.
!    'RSRS' is the default for all equation types.
!  </description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
!  <valid value="RSRS" note="(Red/Black Sweep, Send Receive) repeated twice"/>
!  <valid value="ISIS" note="(Sweep in I, Send Receive) repeated twice"/>
!  <valid value="JSJS" note="(Sweep in J, Send Receive) repeated twice"/>
!  <valid value="KSKS" note="(Sweep in K, Send Receive) repeated twice"/>
!  <valid value="ASAS" note="(All Sweep, Send Receive) repeated twice"/>
      LEQ_SWEEP(:) = 'RSRS'
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Linear precondition used by the BiCGSTAB and CG LEQ solvers. 'LINE'
!    is the default for all equation types.
!  </description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
!  <valid value="NONE" note="No preconditioner"/>
!  <valid value="LINE" note="Line relaxation"/>
!  <valid value="DIAG" note="Diagonal Scaling"/>
      LEQ_PC(:) = 'LINE'
!</keyword>


!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Under relaxation factors.
!    o 0.8 for equation types 1,9
!    o 0.5 for equation types 2,3,4,5,8
!    o 1.0 for equation types 6,7,10
!  </description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
      UR_FAC(1)  = 0.8D0     ! pressure
      UR_FAC(2)  = 0.5D0     ! rho, ep
      UR_FAC(3)  = 0.5D0     ! U
      UR_FAC(4)  = 0.5D0     ! V
      UR_FAC(5)  = 0.5D0     ! W
      UR_FAC(6)  = 1.0D0     ! T
      UR_FAC(7)  = 1.0D0     ! X
      UR_FAC(8)  = 0.5D0     ! Th
      UR_FAC(9)  = 0.8D0     ! Scalar
      UR_FAC(10) = 1.0D0     ! DES Diffusion
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    The implicitness calculation of the gas-solids drag coefficient
!    may be underrelaxed by changing ur_f_gs, which takes values
!    between 0 to 1.
!    o  0 updates F_GS every time step
!    o  1 updates F_GS every iteration
!  </description>
!  <range min="0" max="1" />
      UR_F_gs = 1.0D0
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Under relaxation factor for conductivity coefficient associated
!    with other solids phases for IA Theory [1.0].
!  </description>
      UR_Kth_sml = 1.0D0
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>Discretization scheme of equations.</description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
!  <valid value="0" note="First-order upwinding."/>
!  <valid value="1" note="First-order upwinding (using down-wind factors)."/>
!  <valid value="3" note="Smart."/>
!  <valid value="2" note="Superbee (recommended method)."/>
!  <valid value="5" note="QUICKEST (does not work)."/>
!  <valid value="4" note="ULTRA-QUICK."/>
!  <valid value="7" note="van Leer."/>
!  <valid value="6" note="MUSCL."/>
!  <valid value="8" note="minmod."/>
!  <valid value="9" note="Central (often unstable; useful for testing)."/>
      DISCRETIZE(:) = 0
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Use deferred correction method for implementing higher order
!    discretization.
!  </description>
!  <valid value=".FALSE." note="Use down-wind factor method (default)."/>
!  <valid value=".TRUE."  note="Use deferred correction method."/>
      DEF_COR  =  .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    This scheme guarantees that the set of differenced species mass
!    balance equations maintain the property that the sum of species
!    mass fractions sum to one. This property is not guaranteed when
!    a flux limiter is used with higher order spatial discretization
!    schemes. Note: The chi-scheme is implemented for SMART and MUSCL
!    discretization schemes.
!    Darwish, M.S., Moukalled, F. (2003). Computer Methods in Applied
!    Mech. Eng., 192(13):1711-1730.
!  </description>
!  <valid value=".FALSE." note="Do not use the chi-scheme."/>
!  <valid value=".TRUE."  note="Use the chi-scheme correction."/>
      Chi_scheme = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Four point fourth order interpolation and is upstream biased.
!    Notes:
!    o DISCRETIZE(*) defaults to Superbee if this scheme is chosen
!      and DISCRETIZE(*) < 2.
!    o Set C_FAC between 0 and 1 when using this scheme.
!  </description>
!  <dependent keyword="C_FAC" value="DEFINED"/>
!  <valid value=".FALSE." note="Do not use fourth order interpolation."/>
!  <valid value=".TRUE."  note="Use fourth order interpolation."/>
      FPFOI = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Factor between zero and one used in the universal limiter when
!    using four point, fourth order interpolation (FPFOI).
!    o Choosing one gives (diffusive) first order upwinding.
!    o The scheme becomes more compressive as values near zero.
!  </description>
!  <range min="0.0" max="1.0" />
!  <dependent keyword="fpfoi" value=".TRUE."/>
      C_FAC = UNDEFINED
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>Temporal discretization scheme.</description>
!  <valid value=".FALSE."
!    note="Implicit Euler based temporal discretization scheme employed
!      (first order accurate in time)."/>
!  <valid value=".TRUE."
!    note="Two-step implicit Runge-Kutta method based temporal
!      discretization scheme employed. This method should be second
!      order accurate in time excluding pressure terms and restart
!      time step which are first order accurate. However, recent testing
!      shows that second order accuracy in time is not observed."/>
      CN_ON = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    The code declares divergence if the velocity anywhere in the domain
!    exceeds a maximum value.  This maximum value is automatically
!    determined from the boundary values. The user may scale the maximum
!    value by adjusting this scale factor [1.0d0].
!  </description>
      MAX_INLET_VEL_FAC = ONE
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Solve transpose of linear system. (BICGSTAB ONLY).
!  </description>
!  <dependent keyword="LEQ_METHOD" value="2"/>
      DO_TRANSPOSE = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Frequency to check for convergence. (BICGSTAB ONLY)
!  </description>
!  <dependent keyword="LEQ_METHOD" value="2"/>
      icheck_bicgs = 1
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Sets optimal LEQ flags for parallel runs.
!  </description>
      OPT_PARALLEL = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Use do-loop assignment over direct vector assignment.
!  </description>
      USE_DOLOOP = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Calculate dot-products more efficiently (Serial runs only.)
!  </description>
      IS_SERIAL = .TRUE.
!</keyword>


!#####################################################################!
!                      Geometry and Discretization                    !
!#####################################################################!


!<keyword category="Geometry and Discretization" required="false">
!  <description>Coordinates used in the simulation.</description>
!  <valid value="cartesian" note="Cartesian coordinates."/>
!  <valid value="cylindrical" note="Cylindrical coordinates."/>
      COORDINATES = UNDEFINED_C
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>(Do not use.)</description>
!  <valid value=".FALSE." note="x (r) direction is considered."/>
!  <valid value=".TRUE." note="x (r) direction is not considered."/>
!     NO_I = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Number of cells in the x (r) direction.</description>
      IMAX = UNDEFINED_I
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Cell sizes in the x (r) direction. Enter values from DX(0) to
!    DX(IMAX-1).
!    o Use uniform mesh size with higher-order discretization methods.
!    o DX should be kept uniform in cylindrical coordinates
!      for strict momentum conservation.
!  </description>
!  <arg index="1" id="Cell" min="0" max="DIM_I"/>
      DX(:DIM_I) = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    The inner radius in the simulation of an annular cylindrical region.
!  </description>
      XMIN = ZERO
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Reactor length in the x (r) direction.</description>
      XLENGTH = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>(Do not use.)</description>
!  <valid value=".FALSE. note="y-direction is considered."/>
!  <valid value=".TRUE." note="y-direction is not considered."/>
!     NO_J = .FALSE.
!</keyword>


!<keyword category="Geometry and Discretization" required="false">
!  <description>Number of cells in the y-direction.</description>
      JMAX = UNDEFINED_I
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Cell sizes in the y-direction. Enter values from DY(0) to
!    DY(IMAX-1). Use uniform mesh size with second-order
!    discretization methods.
!  </description>
!  <arg index="1" id="Cell" min="0" max="DIM_J"/>
      DY(:DIM_J) = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Reactor length in the y-direction.</description>
      YLENGTH = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag to disable the third dimension (i.e., 2D simulation).
!      o Z axis in Cartesian coordinate system
!      o Theta in Cylindrical coordinate system
!  </description>
!  <valid value=".FALSE." note="3D simulation."/>
!  <valid value=".TRUE."  note="2D simulation."/>
      NO_K = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Number of cells in the z-direction.</description>
      KMAX = UNDEFINED_I
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Cell sizes in the z (theta) direction. Enter values from DZ(0) to
!    DZ(IMAX-1). Use uniform mesh size with second-order discretization
!    methods.
!  </description>
!  <arg index="1" id="Cell" min="0" max="DIM_K"/>
      DZ(:DIM_K) = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Reactor length in the z (theta) direction.</description>
      ZLENGTH = UNDEFINED
!</keyword>


!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag for making the x-direction cyclic without pressure drop. No other
!    boundary conditions for the x-direction should be specified.
!</description>
!  <valid value=".FALSE." note="No cyclic condition at x-boundary."/>
!  <valid value=".TRUE." note="Cyclic condition at x-boundary."/>
      CYCLIC_X = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag for making the x-direction cyclic with pressure drop. If the
!    keyword FLUX_G is given a value this becomes a cyclic boundary
!    condition with specified mass flux. No other boundary conditions
!    for the x-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="No cyclic condition at x-boundary."/>
!  <valid value=".TRUE." note="Cyclic condition with pressure drop at x-boundary."/>
      CYCLIC_X_PD = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Fluid pressure drop across XLENGTH when a cyclic boundary condition
!    with pressure drop is imposed in the x-direction.
!  </description>
      DELP_X = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag for making the y-direction cyclic without pressure drop. No
!    other boundary conditions for the y-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="No cyclic condition at y-boundary."/>
!  <valid value=".TRUE." note="Cyclic condition at x-boundary."/>
      CYCLIC_Y = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag for making the y-direction cyclic with pressure drop. If the
!    keyword FLUX_G is given a value this becomes a cyclic boundary
!    condition with specified mass flux. No other boundary conditions
!    for the y-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="No cyclic condition at y-boundary."/>
!  <valid value=".TRUE." note="Cyclic condition with pressure drop at y-boundary."/>
      CYCLIC_Y_PD = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Fluid pressure drop across YLENGTH when a cyclic boundary condition
!    with pressure drop is imposed in the y-direction.
!  </description>
      DELP_Y = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag for making the z-direction cyclic without pressure drop. No
!    other boundary conditions for the z-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="No cyclic condition at z-boundary."/>
!  <valid value=".TRUE." note="Cyclic condition at z-boundary."/>
      CYCLIC_Z = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag for making the z-direction cyclic with pressure drop. If the
!    keyword FLUX_G is given a value this becomes a cyclic boundary
!    condition with specified mass flux. No other boundary conditions
!    for the z-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="No cyclic condition at z-boundary."/>
!  <valid value=".TRUE." note="Cyclic condition with pressure drop at
!    z-boundary."/>
      CYCLIC_Z_PD = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Fluid pressure drop across ZLENGTH when a cyclic boundary condition
!    with pressure drop is imposed in the z-direction.
!  </description>
      DELP_Z = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Imposes a mean shear on the flow field as a linear function of the
!    x coordinate. This feature should only be used when CYCLIC_X is
!    .TRUE. and the keyword V_SH is set.
!  </description>
!  <dependent keyword="CYCLIC_X" value=".TRUE."/>
!  <dependent keyword="V_SH" value="DEFINED"/>
      SHEAR = .FALSE.
!</keyword>


!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Specifies the mean y velocity component at the eastern boundary
!    of the domain (V_SH), and the mean Y velocity (-V_SH) at the
!    western boundary of the domain.
!  </description>
      V_sh = 0.0d0
!</keyword>


!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    If a value is specified (in units of g/cm^2.s), the domain-averaged gas
!    flux is held constant at that value in simulations over a periodic
!    domain.  A pair of boundaries specified as periodic with fixed
!    pressure drop is then treated as periodic with fixed mass flux.
!    Even for this case a pressure drop must also be specified, which
!    is used as the initial guess in the simulations.
!  </description>
      Flux_g = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Applies the 2.5D model for cylindrical column by combining 2D assumption
!    and axi-symmetric assumption.
!    Li et al. (2015). A 2.5D computational method to simulate
!    cylindrical fluidized beds, Chemical Engineering Science,
!    123:236-246.
!  </description>
      CYLINDRICAL_2D = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Parameter to control the plate half width and the wedge radius
!    in the 2.5D cylindrical model. This value should be less than
!    half the grid cells in the radial direction (IMAX/2).  [1]
!  </description>
!  <dependent keyword="CYLINDRICAL_2D" value=".TRUE."/>
      I_CYL_NUM = 1
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Parameter to smooth the transition from cylindrical to 2D in
!    the 2.5D cylindrical model. [2]
!  </description>
!  <valid value="2" note="Two cell smoothing transition."/>
!  <valid value="1" note="One cell smoothing transition."/>
!  <valid value="0" note="No smoothing."/>
!  <dependent keyword="CYLINDRICAL_2D" value=".TRUE."/>
      I_CYL_TRANSITION = 2
!</keyword>

!#####################################################################!
!                               Gas Phase                             !
!#####################################################################!

!<keyword category="Gas Phase" required="false">
!  <description>
!    Specified constant gas density [g/cm^3 in CGS]. An equation of
!    state -the ideal gas law by default- is used to calculate the gas
!    density if this parameter is undefined. The value may be set to
!    zero to make the drag zero and to simulate granular flow in a
!    vacuum. For this case, users may turn off solving for gas momentum
!    equations to accelerate convergence.
!  </description>
      RO_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>
!    Specified constant gas viscosity [g/(cm.s) in CGS].
!  </description>
      MU_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>
!    Specified constant gas conductivity [cal/(s.cm.K) in CGS].
!  </description>
      K_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>
!    Specified constant gas specific heat [cal/(g.s.K) in CGS].
!  </description>
      C_PG0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>
!    Specified constant gas diffusivity [(cm^2/s) in CGS].
!  </description>
      DIF_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>
!    Average molecular weight of gas [(g/mol) in CGS]. Used in
!    calculating the gas density for non-reacting flows when the gas
!    composition is not defined.
!  </description>
      MW_AVG = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>
!    Molecular weight of gas species [(g/mol) in GCS].
!  </description>
!  <arg index="1" id="Species" min="1" max="DIM_N_G"/>
      MW_G(:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>Number of species comprising the gas phase.</description>
      NMAX_g = UNDEFINED_I
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>
!    Name of gas phase species as it appears in the materials database.
!  </description>
!  <arg index="1" id="Species" min="1" max="DIM_N_G"/>
      SPECIES_g = UNDEFINED_C
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>
!    User defined name for gas phase species. Aliases are used in
!    specifying chemical equations and must be unique.
!  </description>
!  <arg index="1" id="Species" min="1" max="DIM_N_G"/>
      SPECIES_ALIAS_g = UNDEFINED_C
!</keyword>



!#####################################################################!
!                            Solids Phase                             !
!#####################################################################!

!<keyword category="Solids Phase" required="false">
!  <description>
!    Defines the model used for the solids phase. For TFM/DEM
!    hybrid simulations, first define all TFM solids, then
!    define the DEM solids phases.
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <valid value='TFM' note='Two-fluid Model (continuum)' />
!  <valid value='DEM' note='Discrete Element Model' />
!  <valid value='PIC' note='Multiphase-Particle in Cell' />
      SOLIDS_MODEL(:DIM_M) = 'TFM'
!</keyword>

!<keyword category="Solids Phase" required="false"
!  tfm="true" dem="true" pic="true">
!  <description>Number of solids phases.</description>
      MMAX = 1
!</keyword>

!<keyword category="Solids Phase" required="false"
!  tfm="true" dem="true" pic="true">
!  <description>
!    Initial particle diameters [cm in CGS].
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      D_P0(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false"
!  tfm="true" dem="true" pic="true">
!  <description>
!    Specified constant solids density [g/cm^3 in CGS]. Reacting flows
!    may use variable solids density by leaving this parameter
!    undefined and specifying X_S0 and RO_XS0 as well as the index
!    of the inert species.
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      RO_S0(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false" tfm="true" dem="true">
!  <description>
!    Baseline species mass fraction. Specifically, the mass fraction
!    of an unreacted sample (e.g., proximate analysis).
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_s"/>
!  <dependent keyword="SPECIES_EQ" value=".TRUE."/>
!  <dependent keyword="RO_Xs0" value="DEFINED"/>
!  <dependent keyword="INERT_SPECIES" value="DEFINED"/>
!  <conflict keyword="RO_s0" value="DEFINED"/>
      X_s0(:DIM_M,:DIM_N_s) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false" tfm="true" dem="true">
!  <description>
!    Specified constant solids species density [g/cm^3 in CGS].
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_s"/>
!  <dependent keyword="SPECIES_EQ" value=".TRUE."/>
!  <dependent keyword="X_s0" value="DEFINED"/>
!  <dependent keyword="INERT_SPECIES" value="DEFINED"/>
!  <conflict keyword="RO_s0" value="DEFINED"/>
      RO_Xs0(:DIM_M,:DIM_N_s) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false" tfm="true" dem="true">
!  <description>
!    Index of inert solids phase species. This species should not be a
!    product or reactant of any chemical reaction.
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_s"/>
!  <dependent keyword="SPECIES_EQ" value=".TRUE."/>
!  <dependent keyword="X_s0" value="DEFINED"/>
!  <dependent keyword="RO_Xs0" value="DEFINED"/>
!  <conflict keyword="RO_s0" value="DEFINED"/>
      INERT_SPECIES(:DIM_M) = UNDEFINED_I
!</keyword>

!<keyword category="Solids Phase" required="false" tfm="true" dem="true">
!  <description>
!    Mass fraction of inert solids phase species in the dilute region.
!    In dilute region (see DIL_FACTOR_VSD), the solids density is computed based
!    on this inert species mass fraction, rather than the current inert species mass fraction.
!    This may help convergence when the Variable Solids Density model is invoked.
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <dependent keyword="SPECIES_EQ" value=".TRUE."/>
!  <dependent keyword="X_s0" value="DEFINED"/>
!  <dependent keyword="RO_Xs0" value="DEFINED"/>
!  <conflict keyword="RO_s0" value="DEFINED"/>
      DIL_INERT_X_VSD(:DIM_M) = ONE
!</keyword>

!<keyword category="Solids Phase" required="false" tfm="true" dem="true">
!  <description>
!    Factor to define the dilute region where the solids density is set using DIL_INERT_X_VSD.
!    Cells where the solids volume fraction is between DIL_EP_S and DIL_EP_S x DIL_FACTOR_VSD
!    will automatically set the solids density using DIL_INERT_X_VSD instead of the current
!    inerts species mass fraction. Set this factor to zero to always use the current inert
!    species mass fraction.
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <dependent keyword="SPECIES_EQ" value=".TRUE."/>
!  <dependent keyword="X_s0" value="DEFINED"/>
!  <dependent keyword="RO_Xs0" value="DEFINED"/>
!  <conflict keyword="RO_s0" value="DEFINED"/>
      DIL_FACTOR_VSD = 10.0D0
!</keyword>

!<keyword category="Solids Phase" required="false" tfm="true" dem="true">
!  <description>
!    Specified constant solids conductivity [cal/(s.cm.K) in CGS].
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      K_S0(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false" tfm="true" dem="true">
!  <description>
!    Specified constant solids specific heat [cal/(g.s.K) in CGS].
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      C_PS0(:DIM_M) = UNDEFINED
!</keyword>


!<keyword category="Solids Phase" required="false" tfm="true" dem="true">
!  <description>
!    Molecular weight of solids phase species [(g/mol) in CGS].
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_s"/>
      MW_S(:DIM_M,:DIM_N_s) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false" tfm="true" dem="true">
!  <description>
!    Number of species comprising the solids phase.
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      NMAX_s(:DIM_M) = UNDEFINED_I
!</keyword>

!<keyword category="Solids Phase" required="false" tfm="true" dem="true">
!  <description>
!    Name of solids phase M, species N as it appears in the materials
!    database.
!</description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_s"/>
      SPECIES_s(:DIM_M,:DIM_N_s) = UNDEFINED_C
!</keyword>

!<keyword category="Solids Phase" required="false" tfm="true" dem="true">
!  <description>
!    User defined name for solids phase species. Aliases are used in
!    specifying chemical equations and must be unique.
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_s"/>
      SPECIES_ALIAS_s(:DIM_M,:DIM_N_s) = UNDEFINED_C
!</keyword>

!#####################################################################!
!                           Two Fluid Model                           !
!#####################################################################!


!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Solids phase stress model [LUN_1984]. This is only needed when
!    solving the granular energy PDE (GRANULAR_ENERGY = .TRUE.).
!  </description>
!  <dependent keyword="GRANULAR_ENERGY" value=".TRUE."/>
!  <valid value="AHMADI"
!    note="Cao and Ahmadi (1995). Int. J. Multiphase Flow 21(6), 1203."/>
!  <valid value="GD_99"
!     note="Garzo and Dufty (1999). Phys. Rev. E 59(5), 5895."/>
!  <valid value="GHD"
!    note="Garzo, Hrenya and Dufty (2007). Phys. Rev. E 76(3), 31304"/>
!  <valid value="GTSH"
!    note="Garzo, Tenneti, Subramaniam, Hrenya (2012). J.Fluid Mech. 712, 129."/>
!  <valid value="IA_NONEP"
!     note="Iddir & Arastoopour (2005). AIChE J. 51(6), 1620"/>
!  <valid value="LUN_1984"
!    note="Lun et al (1984). J. Fluid Mech., 140, 223."/>
!  <valid value="SIMONIN"
!    note="Simonin (1996). VKI Lecture Series, 1996-2"/>
      KT_TYPE = "LUN_1984"
!</keyword>

! Retired keyword for specifying Ahmadi KT Theory.
! Use: KT_TYPE = "AHMADI"
      AHMADI = .FALSE.

! Retired keyword for specifying Simonin KT Theory.
! Use: KT_TYPE = "SIMONIN"
      SIMONIN = .FALSE.

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Solids stress model selection.
!  </description>
!  <valid value=".FALSE." note="Do not use the Princeton solids stress model."/>
!  <valid value=".TRUE."  note="Use the Princeton solids stress model"/>
!  <dependent keyword="GRANULAR_ENERGY" value=".TRUE."/>
!  <dependent keyword="PHI" value="DEFINED"/>
!  <dependent keyword="PHI_W" value="DEFINED"/>
      FRICTION = .FALSE.
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    For a term appearing in the frictional stress model
!    invoked with FRICTION keyword.
!  </description>
!  <valid value="0" note="Use S:S in the frictional stress model."/>
!  <valid value="1" note="Use an alternate form suggested by Savage."/>
!  <valid value="2" note="An appropriate combination of above."/>
!  <dependent keyword="friction" value=".TRUE."/>
      SAVAGE = 1
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Schaeffer frictional stress tensor formulation. </description>
!  <dependent keyword="PHI" value="DEFINED"/>
!  <valid value=".TRUE." note="Use the Schaeffer model."/>
!  <valid value=".FALSE." note="Do not use the Schaeffer model."/>
      SCHAEFFER = .TRUE.
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Blend the Schaeffer stresses with the stresses resulting from
!    algebraic kinetic theory around the value of EP_STAR.
!  </description>
      BLENDING_STRESS = .FALSE.
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="ture">
!  <description>
!    Hyperbolic tangent function for blending frictional stress models.
!  </description>
!  <dependent keyword="BLENDING_STRESS" value=".TRUE."/>
!  <conflict keyword="SIGM_BLEND" value=".TRUE."/>
      TANH_BLEND = .TRUE.
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    A scaled and truncated sigmoidal function for blending
!    frictional stress models.
!  </description>
!  <dependent keyword="BLENDING_STRESS" value=".TRUE."/>
!  <conflict keyword="TANH_BLEND" value=".TRUE."/>
      SIGM_BLEND = .FALSE.
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Correlation to compute maximum packing for polydisperse systems.
!  </description>
!  <valid value=".TRUE."
!    note="Use the Yu and Standish correlation."/>
!  <valid value=".FALSE."
!    note="Do not use the Yu and Standish correlation."/>
      YU_STANDISH = .FALSE.
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Correlation to compute maximum packing for binary (only)
!    mixtures of powders.
!  </description>
!  <valid value=".TRUE."
!    note="Use the Fedors and Landel correlation."/>
!  <valid value=".FALSE."
!    note="Do not use the Fedors and Landel correlation."/>
      FEDORS_LANDEL = .FALSE.
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Radial distribution function at contact for polydisperse systems.
!    Do not specify any RDF for monodisperse systems because Carnahan-
!    Starling is the model only available.
!
!    Carnahan, N.F. and Starling K.E., (1969).
!    The Journal of Chemical Physics, Vol. 51(2):635-636.
!  </description>
!
!  <valid value="LEBOWITZ" note="Lebowitz, J.L. (1964)
!   The Physical Review, A133, 895-899"/>
!
!  <valid value="MODIFIED_LEBOWITZ" note="
!    Iddir, H. Y., Modeling of the multiphase mixture of particles
!    using the kinetic theory approach. Doctoral Dissertation,
!    Illinois Institute of Technology, Chicago, Illinois, 2004,
!    (chapter 2, equations 2-49 through 2-52.)"/>
!
!  <valid value="MANSOORI" note="
!   Mansoori, GA, Carnahan N.F., Starling, K.E. Leland, T.W. (1971).
!    The Journal of Chemical Physics, Vol. 54:1523-1525."/>
!
!  <valid value="MODIFIED_MANSOORI" note="van Wachem, B.G.M., Schouten, J.C.,
!    van den Bleek, C.M., Krishna, R. and Sinclair, J. L. (2001)
!    AIChE Journal 47:1035–1051."/>
      RDF_TYPE = 'LEBOWITZ'
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Flag to include the added (or virtual) mass force. This force
!    acts to increase the inertia of the dispersed phase, which
!    tends to stabilize simulations of bubbly gas-liquid flows.
!  </description>
!  <dependent keyword="M_AM" value="DEFINED"/>
      Added_Mass = .FALSE.
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    The disperse phase number to which the added mass is applied.
!  </description>
      M_AM = UNDEFINED_I
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Coefficient of restitution for particle-particle collisions.
!  </description>
!  <range min="0.0" max="1.0" />
      C_E = UNDEFINED
!</keyword>

!<keyword category="Two Fluid Model" required="false">
!  <description>
!    Coefficient of restitution for particle-particle collisions specific
!    to GHD theory implementation.
!  </description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <arg index="2" id="Phase" min="0" max="DIM_M"/>
      r_p(:DIM_M, :DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Two Fluid Model" required="false">
!  <description>
!    Coefficient of restitution for particle-wall collisions when using
!    Johnson and Jackson partial slip BC (BC_JJ_PS).</description>
!  <range min="0.0" max="1.0" />
      E_W = 1.D0
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Specularity coefficient associated with particle-wall collisions
!    when using Johnson and Jackson partial slip BC (BC_JJ_PS). If
!    Jenkins small frictional BC are invoked (JENKINS) then phip is
!    not used.
!  </description>
!  <range min="0.0" max="1.0" />
      PHIP = 0.6D0
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Specify the value of specularity coefficient when the normalized
!    slip velocity goes to zero when BC_JJ_M is .TRUE.. This variable
!    is calculated internally in the code. Do not modify unless an
!    accurate number is known.
!  </description>
!  <dependent keyword="BC_JJ_M" value=".TRUE."/>
      phip0 = undefined
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Coefficient of friction between the particles of two solids phases.
!  </description>
      C_F = UNDEFINED
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!     Angle of internal friction (in degrees). Set this value
!     to zero to turn off plastic regime stress calculations.
!  </description>
      PHI = UNDEFINED
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Angle of internal friction (in degrees) at walls. Set this
!    value to non-zero (PHI_W = 11.31 means TAN_PHI_W = MU = 0.2)
!    when using Johnson and Jackson partial slip BC (BC_JJ_PS) with
!    Friction model or Jenkins small frictional boundary condition.
!  </description>
      PHI_W = UNDEFINED
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Minimum solids fraction above which friction sets in. [0.5] (when
!    FRICTION = .TRUE.)
!  </description>
!  <dependent keyword="FRICTION" value=".TRUE."/>
      EPS_F_MIN = 0.5D0
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Maximum solids volume fraction at packing for polydisperse
!    systems (more than one solids phase used). The value of
!    EP_STAR may change during the computation if solids phases
!    with different particle diameters are specified and
!    Yu_Standish or Fedors_Landel correlations are used.
!  </description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <range min="0" max="1-EP_STAR" />
      EP_S_MAX(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Used in calculating the initial slope of segregation: see
!    Gera et al. (2004) - recommended value 0.3. Increasing this
!    coefficient results in decrease in segregation of particles
!    in binary mixtures.
!  </description>
      SEGREGATION_SLOPE_COEFFICIENT=0.D0
!</keyword>


!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>Excluded volume in Boyle-Massoudi stress.</description>
!  <valid value="0.0" note="b-m stress is turned off."/>
      V_EX = ZERO
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Specified constant viscosity. If any value is specified then:
!    1) kinetic theory calculations (granular_energy) are off, which
!       means zero granular pressure contribution (P_S = 0),
!    2) frictional/plastic calculations are off, which means zero
!       frictional viscosity contributions, however, a plastic pressure
!       term is still invoked (P_STAR), and
!    3) LAMBDA_S = -2/3 MU_S0.
!  </description>
!  <conflict keyword="GRANULAR_ENERGY" value=".TRUE."/>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      MU_S0(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Specified constant solids diffusivity [(cm^2)/s in CGS].
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      DIF_S0(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Packed bed void fraction. Used to calculate plastic stresses (for
!    contribution to viscosity) and when to implement plastic pressure,
!    P_STAR. Specifically, if EP_G < EP_STAR, then plastic pressure is
!    employed in the momentum equations.
!  </description>
!  <range min="0.0" max="1.0" />
      EP_STAR = UNDEFINED
!</keyword>

!<keyword category="Two Fluid Model" required="false" tfm="true">
!  <description>
!    Flag to enable/disable a phase from forming a packed bed.
!    Effectively removes plastic pressure term from the solids phase
!    momentum equation.
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <valid value=".TRUE." note="The phase forms a packed bed with void
!    fraction EP_STAR."/>
!  <valid value=".FALSE." note="The phase can exceed close pack conditions
!    so that it maybe behave like a liquid."/>
      CLOSE_PACKED(:DIM_M) = .TRUE.
!</keyword>


!#####################################################################!
!                   Initial Conditions Section                        !
!#####################################################################!


      DO LC = 1, DIMENSION_IC

!<keyword category="Initial Condition" required="false">
!  <description>X coordinate of the west face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>X coordinate of the east face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Y coordinate of the south face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Y coordinate of the north face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Z coordinate of the bottom face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Z coordinate of the top face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>I index of the west-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>I index of the east-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>J index of the south-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>J index of the north-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>K index of the bottom-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>K index of the top-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Type of initial condition. Mainly used in restart runs to overwrite
!    values read from the .RES file by specifying it as _PATCH_. The
!    user needs to be careful when using the _PATCH_ option, since the
!    values from the .RES file are overwritten and no error checking is
!    done for the patched values.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial void fraction in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_EP_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial gas pressure in the IC region. If this quantity is not
!    specified, MFIX will set up a hydrostatic pressure profile,
!    which varies only in the y-direction.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_P_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial solids pressure in the IC region. Usually, this value is
!    specified as zero.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_P_STAR(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Turbulence length scale in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_L_SCALE(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial bulk density (rop_s = ro_s x ep_s) of solids phase-m in the
!    IC region. Users need to specify this IC only for polydisperse flow
!    (MMAX > 1). Users must make sure that summation of ( IC_ROP_s(ic,m)
!    / RO_s(m) ) over all solids phases is equal to ( 1.0 - IC_EP_g(ic)).
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_ROP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial solids volume fraction of solids phase-m in the IC region.
!    This may be specified in place of IC_ROP_s.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_EP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial gas phase temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial solids phase-m temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial solids phase-m granular temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Gas phase radiation coefficient in the IC region. Modify file
!    rdtn2.inc to change the source term.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_GAMA_RG(LC) = ZERO
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Gas phase radiation temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_T_RG(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Solids phase-m radiation coefficient in the IC region. Modify file
!    energy_mod.f to change the source term.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_GAMA_RS(LC,:DIM_M) = ZERO
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Solids phase-m radiation temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_T_RS(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial x-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial x-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial y-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial y-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial z-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_W_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial z-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial mass fraction of gas species.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         IC_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial mass fraction of solids species.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         IC_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial value of Scalar n.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
        IC_SCALAR(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial value of K in K-Epsilon.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_K_Turb_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial value of Epsilon in K-Epsilon.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_E_Turb_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Flag for inflating initial lattice distribution
! to the entire IC region. </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
          IC_DES_FIT_TO_REGION(LC) = .FALSE.
!</keyword>


!<keyword category="Initial Condition" required="false">
!  <description>Flag to specify the initial constant number
! of particles per cell for the PIC method initialization.
!Statistical weight of parcels will be calculated by the code.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <dependent keyword="SOLIDS_MODEL" value="PIC"/>
!  <conflict keyword="IC_PIC_CONST_STATWT" value="DEFINED"/>
          IC_PIC_CONST_NPC(LC, :DIM_M) = 0
!</keyword>


!<keyword category="Initial Condition" required="false">
!  <description>Flag to specify the initial constant statistical
! weight for computational particles/parcels. Actual number of
! parcels will be automatically computed. </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <dependent keyword="SOLIDS_MODEL" value="PIC"/>
!  <conflict keyword="IC_PIC_CONST_NPC" value="DEFINED"/>
          IC_PIC_CONST_STATWT(LC, :DIM_M) = ZERO
!</keyword>
      ENDDO




!#####################################################################!
!                        Boundary Conditions                          !
!#####################################################################!
      DO LC = 1, DIMENSION_BC


!<keyword category="Boundary Condition" required="false">
!  <description>X coordinate of the west face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>X coordinate of the east face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y coordinate of the south face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y coordinate of the north face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z coordinate of the top face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>I index of the west-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>I index of the east-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>J index of the south-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>J index of the north-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>K index of the bottom-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>K index of the top-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Type of boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!
!  <valid value='DUMMY'
!    note='The specified boundary condition is ignored. This is
!      useful for turning off some boundary conditions without having
!      to delete them from the file.' />
!
!  <valid value='MASS_INFLOW' alias='MI'
!    note='Mass inflow rates for gas and solids phases are
!      specified at the boundary.'/>
!
!  <valid value='MASS_OUTFLOW' alias='MO'
!    note='The specified values of gas and solids mass outflow
!      rates at the boundary are maintained, approximately. This
!      condition should be used sparingly for minor outflows, when
!      the bulk of the outflow is occurring through other constant
!      pressure outflow boundaries.' />
!
!  <valid value='P_INFLOW' alias='PI'
!    note='Inflow from a boundary at a specified constant
!      pressure. To specify as the west, south, or bottom end of
!      the computational region, add a layer of wall cells to the
!      west, south, or bottom of the PI cells. Users need to specify
!      all scalar quantities and velocity components. The specified
!      values of fluid and solids velocities are only used initially
!      as MFIX computes these values at this inlet boundary.' />
!
!  <valid value='P_OUTFLOW' alias='PO'
!    note='Outflow to a boundary at a specified constant pressure.
!      To specify as the west, south, or bottom end of the computational
!      region, add a layer of wall cells to the west, south, or bottom of
!      the PO cells.' />
!
!  <valid value='FREE_SLIP_WALL' alias='FSW'
!    note='Velocity gradients at the wall vanish. If BC_JJ_PS is
!      equal to 1, the Johnson-Jackson boundary condition is used for
!      solids.  A FSW is equivalent to using a PSW with hw=0.' />
!
!  <valid value='NO_SLIP_WALL' alias='NSW'
!    note='All components of the velocity vanish at the wall. If
!      BC_JJ_PS is equal to 1, the Johnson-Jackson boundary condition is
!      used for solids.  A NSW is equivalent to using a PSW with vw=0
!      and hw undefined.' />
!
!  <valid value='PAR_SLIP_WALL' alias='PSW'
!    note='Partial slip at the wall implemented as
!      dv/dn + hw (v - vw) = 0, where n is the normal pointing from the
!      fluid into the wall. The coefficients hw and vw should be
!      specified. For free slip set hw = 0. For no slip leave hw
!      undefined (hw=+inf) and set vw = 0. To set hw = +inf, leave it
!      unspecified. If BC_JJ_PS is equal to 1, the Johnson-Jackson
!      boundary condition is used for solids.' />
         BC_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase hw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_HW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase hw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_HW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Uw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_UW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase Uw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_UW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Vw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_VW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase Vw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_VW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Ww for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_WW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase Ww for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_WW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!   Johnson and Jackson partial slip BC.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <valid value='0'
!    note='Do not use Johnson and Jackson partial slip bc. Default
!      if granular energy transport equation is not solved.'/>
!  <valid value='1'
!    note='Use Johnson and Jackson partial slip bc. Default if
!      granular energy transport equation is solved.'/>
!  <dependent keyword="GRANULAR_ENERGY" value=".TRUE."/>
         BC_JJ_PS(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Use a modified version of Johnson and Jackson
!   partial slip BC (BC_JJ_PS BC) with a variable specularity
!   coefficient.
!  </description>
!  <dependent keyword="E_w" value="DEFINED"/>
!  <dependent keyword="PHI_w" value="DEFINED"/>
!  <conflict keyword="JENKINS" value=".TRUE."/>
         BC_JJ_M = .FALSE.
!</keyword>

!<keyword category="Two Fluid Model" required="false">
!  <description>
!    This flag effects how the momentum and granular energy boundary
!    conditions are implemented when using BC_JJ_PS BC.
!  </description>
!  <dependent keyword="PHI_w" value="DEFINED"/>
!  <dependent keyword="E_w" value="DEFINED"/>
!  <conflict keyword="BC_JJ_M" value=".TRUE."/>
!  <valid value=".FALSE." note="Use standard boundary conditions."/>
!  <valid value=".TRUE."
!    note="Use Jenkins small frictional boundary condition."/>
         JENKINS = .FALSE.
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified wall value, THETAw_M, in diffusion boundary condition:
!    d(Theta_M)/dn + Hw (THETA_M - THETAw_M) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_THETAW_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Transfer coefficient, Hw, in diffusion boundary condition:
!    d(Theta_M)/dn + Hw (THETA_M - THETAw_M) = C, where n is the fluid-to-wall normal.
!  </description>
!  <description>Hw for granular energy bc.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_HW_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant flux, C, in diffusion boundary condition:
!    d(Theta_M)/dn + Hw (THETA_M - THETAw_M) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_C_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Gas phase heat transfer coefficient, Hw, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_HW_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified gas phase wall temperature, Tw_g, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_TW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant gas phase heat flux, C, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_C_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Solids phase heat transfer coefficient, Hw, in diffusion boundary condition:
!    d(T_s)/dn + Hw (T_s - Tw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <description>Solids phase hw for heat transfer.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_HW_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified solids phase wall temperature, Tw_s, in diffusion boundary condition:
!    d(T_s)/dn + Hw (T_s - Tw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_TW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant solids phase heat flux, C, in diffusion boundary condition:
!    d(T_s)/dn + Hw (T_s - Tw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_C_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Gas phase species mass transfer coefficient, Hw, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_HW_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified wall gas species mass fraction, Xw, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <description>Gas phase Xw for mass transfer.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_XW_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant gas species mass flux, C, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_C_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Solid phase species mass transfer coefficient, Hw, in diffusion boundary condition:
!    d(X_s)/dn + Hw (X_s - Xw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         BC_HW_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified solids species mass fraction at the wall, Xw, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         BC_XW_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant solids species mass flux, C, in diffusion boundary condition:
!    d(X_s)/dn + Hw (X_s - Xw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         BC_C_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Scalar transfer coefficient, Hw, in diffusion boundary condition:
!    d(Scalar)/dn + Hw (Scalar - ScalarW) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
         BC_HW_Scalar(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified scalar value at the wall, ScalarW, in diffusion boundary condition:
!    d(Scalar)/dn + Hw (Scalar - ScalarW) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
         BC_ScalarW(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant scalar flux, C, in diffusion boundary condition:
!    d(Scalar)/dn + Hw (Scalar - ScalarW) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
         BC_C_Scalar(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Void fraction at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_EP_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas pressure at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_P_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Bulk density of solids phase at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_ROP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids volume fraction at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_EP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase temperature at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase-m temperature at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase-m granular temperature at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Mass fraction of gas species at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Mass fraction of solids species at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         BC_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>X-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>X-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_W_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas volumetric flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_VOLFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids volumetric flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_VOLFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas mass flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_MASSFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids mass flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_MASSFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>The interval at the beginning when the normal
!    velocity at the boundary is equal to BC_Jet_g0. When restarting,
!    run this value and BC_Jet_g0 should be specified such that the
!    transient jet continues correctly. MFIX does not store the jet
!    conditions. For MASS_OUTFLOW boundary conditions, BC_DT_0 is
!    the time period to average and print the outflow rates. The
!    adjustment of velocities to get a specified mass or volumetric
!    flow rate is based on the average outflow rate.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_DT_0(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Value of normal velocity during the initial interval BC_DT_0.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_JET_G0(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>The interval when normal velocity is equal to BC_Jet_gh.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_DT_H(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Value of normal velocity during the interval BC_DT_h.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_JET_GH(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>The interval when normal velocity is equal to BC_JET_gL.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_DT_L(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Value of normal velocity during the interval BC_DT_L.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_JET_GL(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Boundary value for user-defined scalar equation.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
         BC_Scalar(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Boundary value of K for K-Epsilon Equation.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_K_Turb_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Boundary value of Epsilon for K-Epsilon Equation.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_E_Turb_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Magnitude of gas velocity in a specified boundary region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_VELMAG_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Magnitude of gas velocity in a specified boundary region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_VELMAG_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Flag to specify the constant number
! of computational particles per cell for the PIC solids inflow BC.
!Statistical weight of parcels will be calculated by the code.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <conflict keyword="BC_PIC_CONST_STATWT" value="DEFINED"/>
!  <dependent keyword="SOLIDS_MODEL" value="PIC"/>
          BC_PIC_MI_CONST_NPC(LC, :DIM_M) = 0
!</keyword>


!<keyword category="Boundary Condition" required="false">
!  <description>Flag to specify the constant statistical
! weight for inflowing computational particles/parcels. Actual number of
! parcels will be automatically computed. </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <conflict keyword="IC_PIC_CONST_NPC" value="DEFINED"/>
          BC_PIC_MI_CONST_STATWT(LC, :DIM_M) = ZERO
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Flag to make the PO BC invisible to discrete solids.
! Set this flag to.FALSE.to remove this BC for discrete solids. </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_IC"/>
         BC_PO_APPLY_TO_DES(LC) = .TRUE.
!</keyword>


         BC_ROP_G(LC) = UNDEFINED
      ENDDO




!#####################################################################!
!                         Internal Surfaces                           !
!#####################################################################!
      DO LC = 1, DIMENSION_IS


!<keyword category="Internal Surface" required="false">
!  <description>X coordinate of the west face or edge.</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>X coordinate of the east face or edge.</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>Y coordinate of the south face or edge</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>Y coordinate of the north face or edge</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>Z coordinate of the bottom face or edge</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>Z coordinate of the top face or edge</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>I index of the west-most cell.</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>I index of the east-most cell</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>J index of the south-most cell</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>J index of the north-most cell</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>K index of the bottom-most cell</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>K index of the top-most cell</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>Type of internal surface</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
!  <valid value="IMPERMEABLE"
!    note="No gas or solids flow through the surface." alias="IP"/>
!  <valid value="SEMIPERMEABLE" alias='SP'
!    note="Gas flows through the surface with an additional resistance.
!      Solids velocity through the surface is set to zero or to a user-
!      specified fixed value (i.e., solids momentum equation for this
!      direction is not solved)." />
         IS_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>
!    Parameters defining the internal surface. These values need to be
!    specified for semipermeable surfaces only. The thickness used for
!    pressure drop computation is that of the momentum cell (DX_e,
!    DY_n, or DZ_t). To turn off the resistance, use a large value
!    for permeability.
!    o IDX=1: Permeability [1.0E32]
!    o IDX=2: Inertial resistance coefficient [0.0]
!  </description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
!  <arg index="2" id="IDX" min="1" max="2"/>
         IS_PC(LC,1) = UNDEFINED
         IS_PC(LC,2) = ZERO
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>Value of fixed solids velocity through semipermeable surfaces.</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IS_VEL_S(LC,:DIM_M) = ZERO
!</keyword>
      ENDDO


!#####################################################################!
!                     Point Source Mass Inlets                        !
!#####################################################################!
      DO LC = 1, DIMENSION_PS

!<keyword category="Point Source" required="false">
!  <description>X coordinate of the west face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>X coordinate of the east face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Y coordinate of the south face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Y coordinate of the north face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Z coordinate of the top face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>I index of the west-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>I index of the east-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>J index of the south-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>J index of the north-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>K index of the bottom-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>K index of the top-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_K_T(LC) = UNDEFINED_I
!</keyword>


!<keyword category="Point Source" required="false">
!  <description>X-component of incoming gas velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Y-component of incoming gas velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Z-component of incoming gas velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_W_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Gas mass flow rate through the point source.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_MASSFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Temperature of incoming gas.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Gas phase incoming species n mass fraction.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         PS_X_G(LC,:DIM_N_g) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>X-component of incoming solids velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Y-component of incoming solids velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Z-component of incoming solids velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Solids mass flow rate through the point source.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_MASSFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Temperature of incoming solids.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Solids phase incoming species n mass fraction.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         PS_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

      ENDDO


!#####################################################################!
!                          Output Control                             !
!#####################################################################!

!<keyword category="Output Control" required="true">
!  <description>
!    Interval at which restart (.res) file is updated.
!  </description>
      RES_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Interval at which a backup copy of the restart file is created.
!  </description>
      RES_BACKUP_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    The number of backup restart files to retain.
!  </description>
      RES_BACKUPS = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Interval at which .SPX files are updated.
!    o SP1: void fraction (EP_G)
!    o SP2: Gas pressure (P_G) and Solids pressure (P_star)
!    o SP3: Gas velocity (U_G, V_G, W_G)
!    o SP4: Solids velocity (U_S, V_S, W_S)
!    o SP5: Solids bulk density (ROP_s)
!    o SP6: Gas and solids temperature (T_G, T_S)
!    o SP7: Gas and solids mass fractions (X_G, X_S)
!    o SP8: Granular temperature (THETA_M)
!    o SP9: User defined scalars. (SCALAR)
!    o SPA: Reaction Rates (ReactionRates)
!    o SPB: Turbulence quantities (K_TURB_G, E_TURB_G)
!  </description>
!  <arg index="1" id="SP Value" min="1" max="N_SPX"/>
      SPX_DT(:N_SPX) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    The number of user defined chemical reactions stored
!    in the *.SPA file.
!  </description>
      nRR = 0
!</keyword>

!<keyword category="Output Control" required="false">
!  <description> Interval at which standard output (.OUT) file is updated.
!    Only run configuration information is written if left undefined. Otherwise
!    all field variables for the entire domain are written in ASCII
!    format to the .OUT file at OUT_DT intervals.
!  </description>
      OUT_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Number of time steps between .LOG file updates.</description>
      NLOG = 25
!</keyword>

!<keyword category="Output Control" required="false">
!  <description> Display the residuals on the screen and provide
!    messages about convergence on the screen and in the .LOG file.
!  </description>
      FULL_LOG = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Specifies the residuals to display. </description>
!  <arg index="1" id="Residual Index" max="8" min="1"/>
!  <valid value="P0" note="Gas pressure"/>
!  <valid value="PM" note="Solids phase M pressure"/>
!  <valid value="R0" note="Gas density"/>
!  <valid value="RM" note="Solids phase M density"/>
!  <valid value="U0" note="Gas phase U-velocity"/>
!  <valid value="V0" note="Gas phase V-velocity"/>
!  <valid value="W0" note="Gas phase W-velocity"/>
!  <valid value="UM" note="Solids phase M U-velocity"/>
!  <valid value="VM" note="Solids phase M V-velocity"/>
!  <valid value="WM" note="Solids phase M W-velocity"/>
!  <valid value="T0" note="Gas temperature"/>
!  <valid value="TM" note="Solids phase M temperature"/>
!  <valid value="X0NN" note="Gas phase species NN mass fraction"/>
!  <valid value="XMNN" note="Solids phase M species NN mass fraction"/>
!  <valid value="K0" note="K-Epsilon model residuals"/>
      RESID_STRING(:8) = UNDEFINED_C
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Display residuals by equation.  </description>
      GROUP_RESID = .FALSE.
!</keyword>


!<keyword category="Output Control" required="false">
!  <description>
!    Provide detailed logging of negative density errors.
!  </description>
!  <valid value=".FALSE." note="Do not log negative density errors."/>
!  <valid value=".TRUE." note="Log negative density errors."/>
      REPORT_NEG_DENSITY = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Provide detailed logging of zero or negative specific heat errors.
!  </description>
!  <valid value=".FALSE." note="Do not log zero or negative specific heat errors."/>
!  <valid value=".TRUE." note="Log zero or negative specific heat errors."/>
      REPORT_NEG_SPECIFICHEAT = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Frequency to perform an overall species mass balance. Leaving
!    undefined suppresses the mass balance calculations which can
!    slightly extend run time.
!  </description>
      REPORT_MASS_BALANCE_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Output the variable specularity coefficient when BC_JJ_M is
!    .TRUE.. The specularity coefficient will be stored in ReactionRates
!    array for post-processing by post-mfix. User needs to set NRR to 1
!    for this purpose. Be careful with this setting when reacting flow
!    is simulated.
!  </description>
      PHIP_OUT_JJ=.FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Use distributed IO :: Each MPI process generates RES/SPx files.
!  </description>
      bDist_IO = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Restart a serial IO run (only one RES file was created) with
!    distributed IO.
!  </description>
!  <dependent keyword="RUN_TYPE" value="RESTART_2"/>
!  <dependent keyword="bDist_IO" value=".TRUE."/>
      bStart_with_one_RES = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Flag to write variable in NetCDF output file. NetCDF support is not
!    included in MFIX by default. The executable must be compiled and
!    linked with an appropriate NetCDF library to use this functionality.
!
!    Variable Index List:
!     1: void fraction (EP_G)
!     2: Gas pressure (P_G)
!     3: Solids pressure (P_star)
!     4: Gas velocity (U_G, V_G, W_G)
!     5: Solids velocity (U_S, V_S, W_S)
!     6: Solids bulk density (ROP_s)
!     7: Gas temperature (T_G)
!     8: Gas and solids temperature (T_S)
!     9: Gas mass fractions (X_G)
!    10: Solids mass fractions (X_S)
!    11: Granular temperature (THETA_M)
!    12: User defined scalars. (SCALAR)
!    13: Reaction Rates (ReactionRates)
!    14: Turbulence quantities (K_TURB_G, E_TURB_G)
!  </description>
!  <arg index="1" id="NetCDF Variable Reference" max="20" min="1"/>
!  <valid value=".TRUE." note="Write variable in NetCDF output."/>
!  <valid value=".FALSE." note="Do not include variable in NetCDF output."/>
      bWrite_netCDF(:20) = .FALSE.
!</keyword>


!#####################################################################!
!                           UDF  Control                              !
!#####################################################################!

!<keyword category="UDF Control" required="false">
!  <description>
!    Flag to enable user-defined subroutines: USR0, USR1, USR2, USR3,
!    USR0_DES, USR1_DES, USR2_DES, USR3_DES, USR4_DES.
!  </description>
!  <valid value=".TRUE." note="Call user-defined subroutines."/>
!  <valid value=".FALSE." note="Do NOT call user-defined subroutines."/>
      CALL_USR = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>
!    Flag to enable user_defined subroutine, usr_source, for
!    calculating source terms in the indicated equation.
!  </description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
!  <valid value=".TRUE." note="Call user-defined source."/>
!  <valid value=".FALSE." note="MFIX default: No additional source."/>
      CALL_USR_SOURCE(:) = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>
!    Flag to use the User Defined Function, USR_PROP_ROg,
!    in model/usr_prop.f for calculating the gas phase
!    density, RO_g.
!  </description>
!  <valid value=".TRUE." note="Call user-defined function."/>
!  <valid value=".FALSE." note="Use MFIX default calculation."/>
      USR_ROg = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>
!    Flag to use the User Defined Function, USR_PROP_CPg,
!    in model/usr_prop.f for calculating the gas phase
!    constant pressure specific heat, C_pg.
!  </description>
!  <valid value=".TRUE." note="Call user-defined function."/>
!  <valid value=".FALSE." note="Use MFIX default calculation."/>
      USR_CPg = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>
!    Flag to use the User Defined Function, USR_PROP_Kg,
!    in model/usr_prop.f for calculating the gas phase
!    conductivity, K_g.
!  </description>
!  <valid value=".TRUE." note="Call user-defined function."/>
!  <valid value=".FALSE." note="Use MFIX default calculation."/>
      USR_Kg = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>
!    Flag to use the User Defined Function, USR_PROP_Difg,
!    in model/usr_prop.f for calculating the gas phase
!    diffusivity, Dif_g.
!  </description>
!  <valid value=".TRUE." note="Call user-defined function."/>
!  <valid value=".FALSE." note="Use MFIX default calculation."/>
      USR_Difg = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>
!    Flag to use the User Defined Function, USR_PROP_Mug,
!    in model/usr_prop.f for calculating the gas phase
!    viscosity, Mu_g.
!  </description>
!  <valid value=".TRUE." note="Call user-defined function."/>
!  <valid value=".FALSE." note="Use MFIX default calculation."/>
      USR_Mug = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false" tfm="true">
!  <description>
!    Flag to use the User Defined Function, USR_PROP_ROs,
!    in model/usr_prop.f for calculating the solids phase
!    density, RO_s.
!  </description>
!  <valid value=".TRUE." note="Call user-defined function."/>
!  <valid value=".FALSE." note="Use MFIX default calculation."/>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      USR_ROs(:DIM_M) = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false" tfm="true">
!  <description>
!    Flag to use the User Defined Function, USR_PROP_CPs,
!    in model/usr_prop.f for calculating the solids phase
!    constant pressure specific heat, C_ps.
!  </description>
!  <valid value=".TRUE." note="Call user-defined function."/>
!  <valid value=".FALSE." note="Use MFIX default calculation."/>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      USR_CPs(:DIM_M) = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false" tfm="true">
!  <description>
!    Flag to use the User Defined Function, USR_PROP_Ks,
!    in model/usr_prop.f for calculating the solids phase
!    conductivity, K_s.
!  </description>
!  <valid value=".TRUE." note="Call user-defined function."/>
!  <valid value=".FALSE." note="Use MFIX default calculation."/>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      USR_Ks(:DIM_M) = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false" tfm="true">
!  <description>
!    Flag to use the User Defined Function, USR_PROP_Difs,
!    in model/usr_prop.f for calculating the solids phase
!    diffusivity, Dif_s.
!  </description>
!  <valid value=".TRUE." note="Call user-defined function."/>
!  <valid value=".FALSE." note="Use MFIX default calculation."/>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      USR_Difs(:DIM_M) = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false" tfm="true">
!  <description>
!    Flag to use the User Defined Function, USR_PROP_Mus,
!    in model/usr_prop.f for calculating the solids phase
!    viscosity, Mu_s; second viscosity, lambda_s; and pressure,
!    P_s.
!  </description>
!  <valid value=".TRUE." note="Call user-defined function."/>
!  <valid value=".FALSE." note="Use MFIX default calculation."/>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      USR_Mus(:DIM_M) = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false" tfm="true">
!  <description>
!    Flag to use the User Defined Function, USR_PROP_Gama,
!    in model/usr_prop.f for calculating the gas-solids phase
!    heat transfer coefficient, Gama_gs.
!  </description>
!  <valid value=".TRUE." note="Call user-defined function."/>
!  <valid value=".FALSE." note="Use MFIX default calculation."/>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      USR_Gama(:DIM_M) = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false" tfm="true">
!  <description>
!    Flag to use the User Defined Function, USR_PROP_Fgs, in
!    model/usr_prop.f for calculating the gas-solids phase drag
!    coefficient due to relative velocity differences, F_gs.
!    Currently unavailable.
!  </description>
!  <valid value=".TRUE." note="Call user-defined function."/>
!  <valid value=".FALSE." note="Use MFIX default calculation."/>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      USR_Fgs(:DIM_M) = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false" tfm="true">
!  <description>
!    Flag to use the User Defined Function, USR_PROP_Fss, in
!    model/usr_prop.f for calculating the solids-solids phase
!    drag coefficient due to relative velocity differences, F_ss.
!    Currently unavailable.
!  </description>
!  <valid value=".TRUE." note="Call user-defined function."/>
!  <valid value=".FALSE." note="Use MFIX default calculation."/>
!  <arg index="1" id="Phase" min="1" max="DIM_LM"/>
      USR_Fss( :((DIM_M*(DIM_M-1)/2)+1) ) = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>User defined constants.</description>
      C(:DIMENSION_C) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Name of user-defined constant. (20 character max)</description>
      C_NAME(:DIMENSION_C) = '....................'
!</keyword>

      DO LC=1, DIMENSION_USR
!<keyword category="UDF Control" required="false">
!  <description>
!    Intervals at which subroutine write_usr1 is called.
!  </description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_DT(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: x coordinate of the west face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: x coordinate of the east face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: y coordinate of the south face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: y coordinate of the north face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: z coordinate of the top face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: i index of the west-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: i index of the east-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: j index of the south-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: j index of the north-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: k index of the bottom-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: k index of the top-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: Type of user-defined output: Binary of ASCII.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook:
!    Variables to be written in the user-defined output files.
!  </description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_VAR(LC) = UNDEFINED_C
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook:
!    Format for writing user-defined (ASCII) output file.
!  </description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_FORMAT(LC) = UNDEFINED_C
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: File extension for the user-defined output.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_EXT(LC) = UNDEFINED_C
!</keyword>
      ENDDO


!#####################################################################!
!                        Chemical Reactions                           !
!#####################################################################!


!<keyword category="Chemical Reactions" required="false">
!  <description>Flag to use stiff chemistry solver (Direct Integration).</description>
!  <conflict keyword="USE_RRATES" value=".TRUE."/>
      STIFF_CHEMISTRY = .FALSE.
!</keyword>

!<keyword category="Chemical Reactions" required="false">
!  <description>
!    Maximum number of internal steps ODEPACK may use to integrate
!    over the time interval. Leaving this value unspecified permits
!    an unlimited number of steps. The stiff solver reports the
!    number of cells that exceed the number of steps as 'incomplete'.
!  </description>
!  <dependent keyword="STIFF_CHEMISTRY" value=".TRUE."/>
!  <conflict keyword="USE_RRATES" value=".TRUE."/>
      STIFF_CHEM_MAX_STEPS = UNDEFINED_I
!</keyword>

!<keyword category="Chemical Reactions" required="false">
!  <description>Flag to use legacy chemical reaction UDFs.</description>
      USE_RRATES = .FALSE.
!</keyword>

!<keyword category="Chemical Reactions" required="false" legacy=.TRUE.>
!  <description>
!    Names of gas and solids phase species as it appears in the
!    materials database. The first NMAX(0) are the names of gas
!    species. The next NMAX(1) are the names of solids phase-1
!    species, etc.
!  </description>
!  <dependent keyword="USE_RRATES" value=".TRUE."/>
      SPECIES_NAME(:DIM_N_ALL) = UNDEFINED_C
!</keyword>

!<keyword category="Chemical Reactions" required="false">
!  <description>
!    Number of species in phase m. Note that the gas phase is indicated
!    as m=0.
!  </description>
!  <dependent keyword="USE_RRATES" value=".TRUE."/>
      NMAX = UNDEFINED_I
!</keyword>


!#####################################################################!
!                    Parallelization Control                          !
!#####################################################################!


!<keyword category="Parallelization Control" required="false">
!  <description>Number of grid blocks in x-direction.</description>
      NODESI = UNDEFINED_I
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>Number of grid blocks in y-direction.</description>
      NODESJ = UNDEFINED_I
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>Number of grid blocks in z-direction.</description>
      NODESK = UNDEFINED_I
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>Print out additional statistics for parallel runs</description>
      solver_statistics = .FALSE.
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>Group residuals to reduce global collectives.</description>
      DEBUG_RESID = .TRUE.
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>All ranks write error messages.</description>
      ENABLE_DMP_LOG = .FALSE.
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>Print the index layout for debugging.</description>
      DBGPRN_LAYOUT = .FALSE.
!</keyword>


!#####################################################################!
!                       Batch Queue Environment                       !
!#####################################################################!


!<keyword category="Batch Queue Environment" required="false">
!  <description>
!    Enables controlled termination feature when running under batch
!    queue system to force MFIX to cleanly terminate before the end
!    of wall clock allocated in the batch session.
!  </description>
      CHK_BATCHQ_END = .FALSE.
!</keyword>

!<keyword category="Batch Queue Environment" required="false">
!  <description>Total wall-clock duration of the job, in seconds.</description>
      BATCH_WALLCLOCK = 9000.0    ! set to 2.5 hrs for jaguarcnl w/ nproc<=512
!</keyword>

!<keyword category="Batch Queue Environment" required="false">
!  <description>
!    Buffer time specified to allow MFIX to write out the files and
!    cleanly terminate before queue wall clock time limit is reached
!    such that (BATCH_WALLCLOCK-TERM_BUFFER) is less than then batch
!    queue wall clock time limit, in seconds.
!  </description>
      TERM_BUFFER = 180.0         ! set to 3 minutes prior to end of job
!</keyword>



!#####################################################################!
!          Direct Quadrature Method of Moments (DQMOM)                !
!#####################################################################!


!<keyword category="Direct Quadrature Method of Moments (DQMOM)" required="false">
!  <description>Variable to decide if the population balance equations are solved.</description>
      Call_DQMOM = .FALSE.
!</keyword>

!<keyword category="Direct Quadrature Method of Moments (DQMOM)" required="false">
!  <description>Success-factor for aggregation.</description>
      AGGREGATION_EFF=0.D0
!</keyword>

!<keyword category="Direct Quadrature Method of Moments (DQMOM)" required="false">
!  <description>Success-factor for breakage.</description>
      BREAKAGE_EFF=0.D0
!</keyword>








! ---------------------------------- questionable namelist entries below








!<keyword category="category name" required="false">
!  <description>Variable which triggers an automatic restart.</description>
      AUTOMATIC_RESTART = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>AUTO_RESTART counter.</description>
      ITER_RESTART = 1
!</keyword>



! NO_OF_RXNS is not a keyword. However, it is initialized here so that
! if there are no reactions, this value is assigned.
      NO_OF_RXNS = UNDEFINED_I


      U_G0 = UNDEFINED
      V_G0 = UNDEFINED
      W_G0 = UNDEFINED
      U_S0(:DIM_M) = UNDEFINED
      V_S0(:DIM_M) = UNDEFINED
      W_S0(:DIM_M) = UNDEFINED


      PHIP_OUT_ITER=0





      CALL DES_INIT_NAMELIST

      CALL QMOMK_INIT_NAMELIST

      CALL USR_INIT_NAMELIST

      CALL CARTESIAN_GRID_INIT_NAMELIST

      RETURN
      END SUBROUTINE INIT_NAMELIST
