!
!  First, define the user-defined namelist variables as DOUBLE PRECISION,
!  INTEGER, LOGICAL, or CHARACTER.  Second, include the variable in the
!  appropriate common.  Third, add the variable name in the NAMELIST section.
!  A dummy variable is defined, so that the code can be compiled with the
!  template include file.
!
      DOUBLE PRECISION DUMMY_DPA
!
!     Double Precision variables
      COMMON / USR_DATA_DP / DUMMY_DPA
!
!     Integer variables
!      COMMON / USR_DATA_I /
!
!     Logical variables
!      COMMON / USR_DATA_L /
!
!     Character variables
!      COMMON / USR_DATA_C /
!
!
      NAMELIST / USR_INPUT_DATA / DUMMY_DPA
