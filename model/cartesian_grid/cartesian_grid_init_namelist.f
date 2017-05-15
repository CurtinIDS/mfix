MODULE CG_INIT_NAMELIST
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C
!     Module name: CARTESIAN_GRID_INIT_NAMELIST                           C
!     Purpose: initialize the cartesian_grid-namelist                     C
!                                                                         C
!                                                                         C
!     Author: Jeff Dietiker                              Date: 26-Aug-08  C
!     Reviewer:                                          Date:            C
!     Comments:                                                           C
!                                                                         C
!                                                                         C
!  Keyword Documentation Format:                                          C
!<keyword category="category name" required="true/false"                  C
!                                    legacy="true/false">                 C
!  <description></description>                                            C
!  <arg index="" id="" max="" min=""/>                                    C
!  <dependent keyword="" value="DEFINED"/>                                C
!  <conflict keyword="" value="DEFINED"/>                                 C
!  <valid value="" note="" alias=""/>                                     C
!  <range min="" max="" />                                                C
!  MFIX_KEYWORD=INIT_VALUE                                                C
!</keyword>                                                               C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CARTESIAN_GRID_INIT_NAMELIST

      USE param1
      USE quadric
      USE cutcell
      USE polygon
      USE vtk
      USE progress_bar
      USE dashboard
      Use stl

      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!-----------------------------------------------
!
!
      INCLUDE 'cartesian_grid_namelist.inc'

!<keyword category="cartesian grid" required="false">
!  <description>Activate Cartesian grid cut cell technique.</description>
!  <dependent keyword="COORDINATES" value="CARTESIAN"/>
!  <conflict keyword="COORDINATES" value="CYLINDRICAL"/>
!  <valid value=".false." note="Do not use Cartesian grid cut cell technique."/>
!  <valid value=".true." note="Use Cartesian grid cut cell
!  technique. one of the following methods must be used to define the
!  geometry:"/>
      CARTESIAN_GRID = .FALSE.
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Number of quadric surfaces defining the boundaries (<=100).</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <range min="0" max="100" />
      N_QUADRIC = 0
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Use STL file to describe geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <valid value=".false." note="Do not use STL file."/>
!  <valid value=".true." note="Read triangulated geometry (for 3d geometry only) from geometry.stl."/>
      USE_STL = .FALSE.
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Use .msh file to describe geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <valid value=".false." note="Do not use .msh file."/>
!  <valid value=".true." note="Read geometry (for 3d geometry only) from geometry.msh."/>
      USE_MSH = .FALSE.
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Use polygons to describe geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <valid value=".false." note="Do not use polygons."/>
!  <valid value=".true." note="Read polygon data (for 2d geometry only) from poly.dat."/>
      USE_POLYGON = .FALSE.
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Number of user-defined functions (currently limited to
!  0 or 1). If set to 1, the geometry is defined in the user
!  subroutine eval_usr_fct.f.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <valid value="0" note="Do not use user-defined function" alias=""/>
!  <valid value="1" note="Use one user-defined function" alias=""/>
!  <range min="0" max="1" />
      N_USR_DEF = 0
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Form of the quadric surface equation.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <valid value="normal" note="Use normal form, as defined in equation
!    (1). The LAMDBA\'s and D must be defined"/>
!  <valid value="plane" note="Plane. Needs to define N_X,N_Y,N_Z
!    (unit normal vector pointing away from fluid cells)."/>
!  <valid value="x_cyl_int" note="Cylinder aligned with x-axis,
!    internal flow. Needs to define RADIUS(QID)."/>
!  <valid value="x_cyl_ext" note="Cylinder aligned with x-axis,
!    external flow. Needs to define RADIUS(QID)."/>
!  <valid value="y_cyl_int" note="Cylinder aligned with y-axis,
!    internal flow. Needs to define RADIUS(QID)."/>
!  <valid value="y_cyl_ext" note="Cylinder aligned with y-axis,
!    external flow. Needs to define RADIUS(QID)."/>
!  <valid value="z_cyl_int" note="Cylinder aligned with z-axis,
!    internal flow. Needs to define RADIUS(QID)."/>
!  <valid value="z_cyl_ext" note="Cylinder aligned with z-axis,
!    external flow. Needs to define RADIUS(QID)."/>
!  <valid value="x_cone" note="Cone aligned with x-axis, internal flow.
!    Needs to define HALF_ANGLE(QID)."/>
!  <valid value="y_cone" note="Cone aligned with y-axis, internal flow.
!    Needs to define HALF_ANGLE(QID)."/>
!  <valid value="z_cone" note="Cone aligned with z-axis, internal flow.
!    Needs to define HALF_ANGLE(QID)."/>
!  <valid value="sphere_int" note="Sphere, internal flow.
!    Needs to define RADIUS(QID)."/>
!  <valid value="sphere_ext" note="Sphere, external flow.
!    Needs to define RADIUS(QID)."/>
!  <valid value="C2C" note="Cylinder-to-cylinder conical junction,
!    internal flow. Needs to be defined between two cylinders."/>
!  <valid value="Torus_int" note="Torus, internal flow. Needs to
!    define TORUS_R1(QID) and TORUS_R2(QID).A torus is not a quadric
!    surface but is defined as a basic shape."/>
!  <valid value="Torus_ext" note="Torus, external flow. Needs to define
!    TORUS_R1(QID) and TORUS_R2(QID)."/>
!  <valid value="Y_UCOIL_EXT" note="Pair of parallel cylinders (y-direction),
!     capped at both ends by half a torus, to create a U-shaped coil. Needs
!     UCOIL_R1, UCOIL_R2, UCOIL_Y1, UCOIL_Y2."/>
!  <valid value="Y_UCOIL_EXT" note="Pair of parallel cylinders (y-direction),
!     capped at both ends by a cylinder at 90 degree angle
!     to create a U-shaped coil. Needs
!     UCOIL_R1, UCOIL_R2, UCOIL_Y1, UCOIL_Y2."/>
!  <valid value="XY_BEND_INT" note="Bend between
!     two cylinders in the XY plane, Needs
!     BEND_R1,BEND_R2,BEND_THETA1,BEND_THETA2."/>
!  <valid value="Y_C2C_INT" note="connects two vertical cylinders by a
!     conical section. Needs C2C_R1,C2C_R2,C2C_Y1,C2C_Y2."/>
!  <valid value="REACTOR1" note="Reactor, made of two vertical cylinders,
!     connected by a conical section.Each cylinder is rounded and closed
!     by a conical cap. Needs REACTOR1_R1,REACTOR1_R2,REACTOR1_Y1,REACTOR1_Y2,
!     REACTOR1_YR1,REACTOR1_YR2,REACTOR1_RR1,REACTOR1_RR2,
!     REACTOR1_THETA1,REACTOR1_THETA2."/>
      quadric_form = 'NORMAL'
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Scaling factor, applied to all quadric geometry parameters. Must
!    be a positive number.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <range min="0.0" max="" />
      quadric_scale = ONE
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Coefficient LAMBDA_X in equation (1) ('NORMAL' form) or
!    x-component of normal vector defining plane in equation (5)
!    ('DEGENERATE' form).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      lambda_x = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Coefficient LAMBDA_Y in equation (1) ('NORMAL' form) or
!    y-component of normal vector defining plane in equation (5)
!    ('DEGENERATE' form).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      lambda_y = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Coefficient LAMBDA_Z in equation (1) ('NORMAL' form) or
!    z-component of normal vector defining plane in equation (5)
!    ('DEGENERATE' form).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      lambda_z = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Coefficient D in equation (1).</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      dquadric = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Rotation angle with respect to x-axis (degrees).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      theta_x = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Rotation angle with respect to y-axis (degrees).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      theta_y = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Rotation angle with respect to z-axis (degrees).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      theta_z = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Cylinder radius (used when QUADRIC_FORM = *_CYL_***)
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      Radius = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Cone half angle, expressed in degrees (used when
!    QUADRIC_FORM = *_CONE)
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      Half_angle = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Torus Radius 1 (used when QUADRIC_FORM = TORUS_*),
!    R1>R2 for a ring.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      Torus_R1 = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Torus Radius 2 (used when QUADRIC_FORM = TORUS_*),
!    R1>R2 for a ring.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      Torus_R2 = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    U-shaped coil Radius 1 (used when QUADRIC_FORM = UCOIL*),
!    UCOIL_R1>UCOIL_R2.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      UCOIL_R1 = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    U-shaped coil Radius 2 (used when QUADRIC_FORM = UCOIL*),
!    UCOIL_R1>UCOIL_R2.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      UCOIL_R2 = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    U-shaped coil ymax (used when QUADRIC_FORM = UCOIL*),
!    UCOIL_Y2>UCOIL_Y1.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      UCOIL_Y1 = -UNDEFINED
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    U-shaped coil ymin (used when QUADRIC_FORM = UCOIL*),
!    UCOIL_Y2>UCOIL_Y1.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      UCOIL_Y2 = UNDEFINED
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Bend Radius 1 (used when QUADRIC_FORM = BEND*),
!    BEND_R1>BEND_R2.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BEND_R1 = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Bend Radius 2 (used when QUADRIC_FORM = BEND*),
!    BEND_R1>BEND_R2.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BEND_R2 = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Bend start angle, in degrees (used when QUADRIC_FORM = BEND*).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BEND_THETA1 = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Bend end angle, in degrees (used when QUADRIC_FORM = BEND*).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BEND_THETA2 = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Cylinder-cone_cylinder Radius 1 (used when QUADRIC_FORM = C2C*).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      C2C_R1 = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Cylinder-cone_cylinder Radius 2 (used when QUADRIC_FORM = C2C*).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      C2C_R2 = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Cylinder-cone_cylinder Y1 (used when QUADRIC_FORM = C2C*).
!    If Y1=Y2, then R1=R2.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      C2C_Y1 = -UNDEFINED
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Cylinder-cone_cylinder Y2 (used when QUADRIC_FORM = C2C*).
!    If Y1=Y2, then R1=R2.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      C2C_Y2 = UNDEFINED
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Reactor 1, lower cylinder radius.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_R1 = UNDEFINED
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Reactor 1, upper cylinder radius.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_R2 = UNDEFINED
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Reactor 1, lower conical transition between cylinders.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_Y1 = UNDEFINED
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Reactor 1, upper conical transition between cylinders.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_Y2 = UNDEFINED
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Reactor 1, lower rounding below cylinder.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_YR1 = UNDEFINED
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Reactor 1, upper rounding above cylinder.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_YR2 = UNDEFINED
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Reactor 1, lower rounding radius.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_RR1 = UNDEFINED
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Reactor 1, upper rounding radius.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_RR2 = UNDEFINED
!</keyword>


!<keyword category="Cartesian grid" required="false">
!  <description>Reactor 1, lower rounding angle (degrees).</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_THETA1 = UNDEFINED
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Reactor 1, upper rounding angle (degrees).</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_THETA2 = UNDEFINED
!</keyword>


!<keyword category="Cartesian grid" required="false">
!  <description>
!    X-component of normal vector defining the plane (used when
!    QUADRIC_FORM = PLANE).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      n_x = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Y-component of normal vector defining the plane (used when
!    QUADRIC_FORM = PLANE).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      n_y = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Z-component of normal vector defining the plane (used when
!    QUADRIC_FORM = PLANE).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      n_z = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Translation in x-direction.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      t_x = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Translation in y-direction.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      t_y = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Translation in z-direction.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      t_z = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Lower x-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      clip_xmin = - LARGE_NUMBER
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Upper x-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      clip_xmax =   LARGE_NUMBER
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Lower y-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      clip_ymin = - LARGE_NUMBER
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Upper y-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      clip_ymax =   LARGE_NUMBER
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Lower z-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      clip_zmin = - LARGE_NUMBER
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Upper z-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      clip_zmax =   LARGE_NUMBER
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Lower x-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      piece_xmin = - LARGE_NUMBER
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Upper z-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      piece_xmax =   LARGE_NUMBER
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Lower y-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      piece_ymin = - LARGE_NUMBER
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Upper y-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      piece_ymax =   LARGE_NUMBER
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Lower z-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      piece_zmin = - LARGE_NUMBER
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Upper z-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      piece_zmax =   LARGE_NUMBER
!</keyword>


!<keyword category="Cartesian grid" required="false">
!  <description>
!    Flag defining the type of cells that are outside of the zone
!    defined by [CLIP_XMIN; CLIP_XMAX], [CLIP_YMIN; CLIP_YMAX],
!    [CLIP_ZMIN; CLIP_ZMAX].
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <valid value=".false." note="Remove cells from computational domain."/>
!  <valid value=".true." note="Treat cells as fluid cells."/>
      FLUID_IN_CLIPPED_REGION = .TRUE.
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Boundary condition flag.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BC_ID_Q = UNDEFINED_I
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Number of group(s) of quadrics (<=50).</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      N_GROUP = 1
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Number of quadrics in the group.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      GROUP_SIZE(1) = 1
!</keyword>
      GROUP_SIZE(2:DIM_GROUP) = 0

!<keyword category="Cartesian grid" required="false">
!  <description>Quadric ID assigned to a group.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      GROUP_Q = 0
!</keyword>

      GROUP_Q(1,1) = 1
!<keyword category="Cartesian grid" required="false">
!  <description>Relation among quadrics of a same group.</description>
!  <valid value="or" note="A point belongs to the computational domain
!    if at least one of f(x,y,z) among all quadrics is negative."/>
!  <valid value="and" note="A point belongs to the computational domain
!    if all of f(x,y,z) among all quadrics are negative."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      GROUP_RELATION = 'OR'
!</keyword>

!<keyword category="Cartesian grid" required="false">

!  <description>
!    Relation between current group and combination of all previous groups.
!  </description>
!  <valid value="or" note="A point belongs to the computational domain
!    if f-value for the current group or f-value for the combination of
!    previous groups is negative."/>
!  <valid value="and" note="A point belongs to the computational
!    domain if f-value for the current group and f-value for the
!    combination of previous groups is negative."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      RELATION_WITH_PREVIOUS = 'OR'
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Tolerance used to snap an intersection point onto an
!  existing cell corner (expressed as a fraction of edge length,
!  between 0.0 and 0.5). For stretched grids, three values can be
!  entered in the x, y and z directions.</description>
!  <range min="0.0" max="0.5" />
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_SNAP(1)    = 0.00D0  ! 0% of original edge length
!</keyword>
      TOL_SNAP(2)    = UNDEFINED
      TOL_SNAP(3)    = UNDEFINED

!<keyword category="Cartesian grid" required="false">
!  <description>Tolerance used to limit acceptable values of normal
!  distance to the wall (expressed as a fraction of cell diagonal,
!  between 0.0 and 1.0).</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_DELH       = 0.00D0  ! 0% of original Diagonal
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Tolerance used to detect small cells (expressed as a
!  fraction of cell volume, between 0.0 and 1.0).</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_SMALL_CELL = 0.01D0  ! 1% of original cell volume
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Tolerance used to remove duplicate nodes (expressed as
!  a fraction of cell diagonal, between 0.0 and 1.0).</description>
      TOL_MERGE      = 1.0D-12 ! fraction of original cell diagonal
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Tolerance used to detect small faces (expressed as a
!  fraction of original face area, between 0.0 and 1.0).</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_SMALL_AREA = 0.01D0  ! 1% of original face area
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Maximum acceptable value of interpolation correction factor.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      ALPHA_MAX      = ONE
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Tolerance used to find intersection of quadric surfaces or user-defined function with background grid.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_F     = 1.0D-9
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Tolerance used to find intersection of polygon with background grid.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_POLY  = 1.0D-9
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Maximum number of iterations used to find intersection points.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      ITERMAX_INT = 10000
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Tolerance used to find intersection of STL triangles with background grid.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_STL = 1.0D-6        ! Settings for STL file
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Smallest angle accepted for valid STL triangles (in
!  degrees). Triangles having an angle smaller that this value will be
!  ignored.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      STL_SMALL_ANGLE = 5.0   ! Degrees
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Dot product tolerance when determining if a point lies in a facet.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_STL_DP = 1.0D-3        ! Settings for STL file
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Maximum number of STL facets per cell.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      DIM_FACETS_PER_CELL  = 10        !
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Defines value of F_STL outside of the STL geometry. a
!  value of 1.0 means the domain outside of the STL geometry is
!  excluded from computation, i.e., an internal flow is
!  computed.</description>
!  <valid value="-1.0" note="model an external flow"/>
!  <valid value="1.0" note="model an internal flow"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      OUT_STL_VALUE = 1.0
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Boundary condition flag for the STL geometry</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      STL_BC_ID = UNDEFINED_I
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Translation in x-direction, applied to the STL geometry.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TX_STL = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Translation in y-direction, applied to the STL geometry.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TY_STL = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Translation in z-direction, applied to the STL geometry.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TZ_STL = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Scaling factor, applied to the STL geometry. Note that translation
!    occurs after scaling.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      SCALE_STL = ONE
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Tolerance used to find intersection of .msh file with background grid.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_MSH = 1.0D-6        ! Settings for MSH file
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Defines value of f outside of the .msh geometry. a
!  value of 1.0 means the domain outside of the .msh geometry is
!  excluded from computation, i.e., an internal flow is
!  computed.</description>
!  <valid value="-1.0" note="model an external flow"/>
!  <valid value="1.0" note="model an internal flow"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      OUT_MSH_VALUE = 1.0
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Translation in x-direction, applied to the .msh geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TX_MSH = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Translation in y-direction, applied to the .msh geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TY_MSH = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Translation in z-direction, applied to the .msh geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TZ_MSH = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Scaling factor, applied to the .msh geometry. Note that translation occurs after scaling.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      SCALE_MSH = ONE
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Ray propagation order used to determine whether any point is
!    located inside or outside of the STL surface.
!  </description>
!  <valid value="ijk" note="Propagation occurs in the
!    I, followed by J, and K directions"/>
!  <valid value="jki" note="Propagation occurs in the
!    J, followed by K, and I directions"/>
!  <valid value="kij" note="Propagation occurs in the
!    K, followed by I, and J directions"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      CAD_PROPAGATE_ORDER = '   '
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Ray direction when propagating CAD value</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      RAY_DIR = 'X-'
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Flag to detect and treat corner cells the same way as
!    in the original MFIX version (i.e. without cut cells).
!  </description>
!  <valid value=".true." note="Some cut cells may be treated as corner cells."/>
!  <valid value=".false." note="Do not treat cut cells as corner cells."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      SET_CORNER_CELLS = .FALSE.
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Factor used to allocate cut cell arrays (expressed as a fraction
!    of DIMENSION_3G).
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      FAC_DIM_MAX_CUT_CELL = 0.25
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write VTK files at regular intervals.</description>
!  <valid value=".false." note="Do not write VTK files. if there are
!  cut cells, they will not be displayed from the usual .res file"/>
!  <valid value=".true." note="Valid only if Cartesian_grid = .true."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      WRITE_VTK_FILES = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Use time-dependent VTK file names</description>
!  <valid value=".false." note="The VTK file overwrites the previous
!    file (recommended for steady-state computation)."/>
!  <valid value=".true." note="A sequential integer is appended to the
!    VTK filenames as they are written to create a series of files
!    (recommended for transient computation)."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TIME_DEPENDENT_FILENAME = .TRUE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Interval (expressed in seconds of simulation time) at which VTK
!    files are written.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_DT     = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>When VTK_DBG_FILE is .TRUE., the VTK region file is only written
!   when the subroutine WRITE_DBG_VTU_AND_VTP_FILES is called, typically in a UDF.
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_DBG_FILE = .FALSE.
!</keyword>



!<keyword category="Output Control" required="false">
!  <description>List of variables written in the VTK files.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <valid value="1" note="Void fraction (EP_g)"/>
!  <valid value="2" note="Gas pressure, solids pressure (P_g, P_star)"/>
!  <valid value="3" note="Gas velocity (U_g, V_g, W_g)"/>
!  <valid value="4" note="Solids velocity (U_s, V_s, W_s)"/>
!  <valid value="5" note="Solids density (ROP_s)"/>
!  <valid value="6" note="Gas and solids temperature (T_g, T_s)"/>
!  <valid value="7" note="Gas and solids mass fractions (X_g, X_s)"/>
!  <valid value="8" note="Granular temperature (Theta_m)"/>
!  <valid value="9" note="Scalar"/>
!  <valid value="10" note="Reaction rates"/>
!  <valid value="11" note="K and Epsilon"/>
!  <valid value="12" note="Vorticity magnitude and lambda_2"/>
!  <valid value="100" note="Grid Partition"/>
!  <valid value="101" note="Boundary Condition ID"/>
!  <valid value="102" note="Distance to wall"/>
!  <valid value="103" note="DEM facet count"/>
!  <valid value="104" note="DEM Neighboring facets"/>
!  <valid value="999" note="Cell IJK index"/>
!  <valid value="1000" note="Cut face normal vector"/>
      VTK_VAR(1) = UNDEFINED_I
!</keyword>
      VTK_VAR(2:20) = UNDEFINED_I

!<keyword category="Output Control" required="false">
!  <description>Write void fraction in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_EP_g = .FALSE.
!</keyword>


!<keyword category="Output Control" required="false">
!  <description>Write gas pressure in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_P_g = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write solids pressure in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_P_star = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write gas velocity vector in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_VEL_g = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write x-component of gas velocity vector in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_U_g = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write y-component of gas velocity vector in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_V_g = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write z-component of gas velocity vector in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_W_g = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write solids velocity vector in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_VEL_S = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write x-component of solids velocity vector in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_U_s = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write y-component of solids velocity vector in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_V_s = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write z-component of solids velocity vector in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_W_s = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write solids bulk density in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_ROP_s = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write gas temperature in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_T_g = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write solids temperature in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_T_s = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write gas phase species in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_X_g = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write solids phase species in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_X_s = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write granular temperature in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_Theta_m = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write scalar in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_Scalar =.FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write reaction rates in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_RRate = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write turbulent kinetic energy in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_K_Turb_G = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write turbulent dissipation rate in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_E_Turb_G = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write vorticity magnitude in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_VORTICITY = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write lambda_2 in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_LAMBDA_2  = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write void grid partition in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_PARTITION = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write boundary condition ID in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_BC_ID = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write wall distance in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_DWALL = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write STL facet count for DES in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_FACET_COUNT_DES = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write neighboring facets in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_NB_FACET_DES = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write cell IJK index in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_IJK = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write cut face normal vector in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_NORMAL = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write debug variable in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_DEBUG = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Type of data to write in the VTK file.</description>
!  <valid value="C" note="Cell data (VTU file)."/>
!  <valid value="P" note="Particle data (VTP file)."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_DATA = 'C'
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Tolerance to detect particle in a VTK region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_SLICE_TOL = ZERO
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Particle selection mode in a VTK region.</description>
!  <valid value="C" note="Select particles with centers inside
!    the VTK region."/>
!  <valid value="P" note="Select particles that are entirely inside
!    the VTK region."/>
!  <valid value="I" note="Select particles that are inside or
!    intersect the edges of the VTK region."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_SELECT_MODE = 'C'
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write particle radius in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_PART_DIAMETER= .FALSE.
!</keyword>


!<keyword category="Output Control" required="false">
!  <description>Write particle velocity in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_PART_VEL= .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write particle angular velocity in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_PART_ANGULAR_VEL= .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write particle angular velocity in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_PART_ORIENTATION= .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write particle user-defined variable in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_PART_USR_VAR= .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write particle temperature in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_PART_TEMP= .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write particle species mass fraction in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_PART_X_s= .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write particle cohesion in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_PART_COHESION= .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Write cut-cell data only in VTK file.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_CUTCELL_ONLY= .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Starting Index appended to VTU files</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      FRAME = -1
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Directory where vtk files are stored (default is run directory)</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTU_DIR = '.'
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>West location of VTK region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
         VTK_X_W = -UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>East location of VTK region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
         VTK_X_E = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>South location of VTK region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
         VTK_Y_S = -UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>North location of VTK region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
         VTK_Y_N = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Bottom location of VTK region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
         VTK_Z_B = -UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>West location of VTK region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
         VTK_Z_T = UNDEFINED
!</keyword>


!<keyword category="Cartesian grid" required="false">
!  <description>
!    Option for pressure gradient computation in cut cells.
!  </description>
!  <valid value="1" note="Use maximum of (east/west), (north/south),
!    and (top/bottom) pairs of velocity cells."/>
!  <valid value="2" note="Use both (east/west), (north/south), and
!    (top/bottom) areas of velocity cells."/>
!  <valid value="0" note="Use east, north and top areas of pressure cell (same as standard cells)."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      PG_OPTION = 0
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Run code in safe mode.</description>
!  <valid value="1" note="Performs initial preprocessing but use all
!    original MFIX subroutines during flow solution (using only cell
!    volumes and areas of cut cells)."/>
!  <valid value="0" note="Runs the code with modified subroutines for
!    cut cell treatment."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      CG_SAFE_MODE = 0
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Prints any warning message encountered during pre-processing
!    on the screen.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      PRINT_WARNINGS = .FALSE.
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Under-relaxation factor used in cut cells (only CG_UR_FAC(2) is used).
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      CG_UR_FAC = 1.0
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Print a progress bar during each major step of pre-processing stage.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      PRINT_PROGRESS_BAR = .FALSE.
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Width of the progress bar (complete status), expressed in number
!    of characters (between 10 and 80).
!  </description>
!  <range min="10" max="80" />
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BAR_WIDTH = 50
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Character used to create the progress bar.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BAR_CHAR  = '='
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Update frequency of progress bar, expressed in percent of total
!    length (between 1.0 and 100.0).
!  </description>
!  <range min="1.0" max="100.0" />
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BAR_RESOLUTION = 5.0
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Writes the file dashboard.txt at regular intervals. The file shows
!    a summary of the simulation progress.
!  </description>
      WRITE_DASHBOARD = .FALSE.
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Frequency, expressed in terms of iterations, at which the dashboard
!    is updated.
!  </description>
      F_DASHBOARD = 1
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Location of control points in x-direction.</description>
      CPX = ZERO
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Number of cells within a segment (x-direction).</description>
      NCX      = 0
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Expansion ratio (last DX/first DX) in a segment (x-direction).
!  </description>
      ERX      = ONE
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Value of first DX in a segment (x-direction). A negative
!      value will copy DX from previous segment (if available).
!  </description>
      FIRST_DX = ZERO
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Value of last DX in a segment (x-direction). A
!  negative value will copy DX from next segment (if
!  available).</description>
      LAST_DX  = ZERO
!</keyword>


!<keyword category="Geometry and Discretization" required="false">
!  <description>Location of control points in y-direction.</description>
      CPY = ZERO
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Number of cells within a segment (y-direction).
!  </description>
      NCY = 0
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Expansion ratio (last DY/first DY) in a segment (y-direction).
!  </description>
      ERY = ONE
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Value of first DY in a segment (y-direction). A negative value will
!    copy DY from previous segment (if available).
!  </description>
      FIRST_DY = ZERO
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Value of last DY in a segment (y-direction). A negative value will
!    copy DY from next segment (if available).
!  </description>
      LAST_DY  = ZERO
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Location of control points in z-direction.</description>
      CPZ = ZERO
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Number of cells within a segment (z-direction).
!  </description>
      NCZ = 0
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Expansion ratio (last DZ/first DZ) in a segment (z-direction).
!  </description>
      ERZ      = ONE
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Value of first DZ in a segment (z-direction). A negative value will
!    copy DZ from previous segment (if available).
!  </description>
      FIRST_DZ = ZERO
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Value of last DZ in a segment (z-direction). A negative value will
!    copy DZ from next segment (if available).
!  </description>
      LAST_DZ  = ZERO
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Turns on the re-indexing of cells. When true, inactive (dead)
!    cells are removed from computational domain.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      RE_INDEXING = .FALSE.
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Attempts to adjust grid partition. Each processor will be assigned
!    its own size to minimize load imbalance.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      ADJUST_PROC_DOMAIN_SIZE = .FALSE.
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Attempts to adjust grid partition. Each processor will
!    be assigned its own size to minimize load imbalance.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REPORT_BEST_DOMAIN_SIZE = .FALSE.
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Temporary setting used in serial run to report best domain size
!    for parallel run.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      NODESI_REPORT = 1
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Temporary setting used in serial run to report best domain size
!    for parallel run.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      NODESJ_REPORT = 1
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Temporary setting used in serial run to report best domain size
!    for parallel run.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      NODESK_REPORT = 1
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>
!    Attempts to minimize the size of the send/receive layers.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      MINIMIZE_SEND_RECV = .TRUE.
!</keyword>

!<keyword category="Cartesian grid" required="false">
!  <description>Brute force calculation of wall distance.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      DWALL_BRUTE_FORCE = .FALSE.
!</keyword>
      RETURN
      END SUBROUTINE CARTESIAN_GRID_INIT_NAMELIST
END MODULE CG_INIT_NAMELIST
