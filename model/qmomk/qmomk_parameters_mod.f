MODULE qmomk_parameters

  IMPLICIT NONE


  !     QMOMK array dimensions

  !     Number of moments
  INTEGER, PARAMETER :: QMOMK_NMOM = 20
  !     Number of nodes and abscissas
  INTEGER, PARAMETER :: QMOMK_NN = 8
  !     Dimension of the covariance matrix
  INTEGER, PARAMETER :: QMOMK_NCOV = 3

  !     Eps values
  DOUBLE PRECISION, PARAMETER  :: eps = 1.D-6
  DOUBLE PRECISION, PARAMETER  :: eps2 = 1.D-4
  DOUBLE PRECISION, PARAMETER  :: epsu = 1.D0
  DOUBLE PRECISION, PARAMETER  :: epsv = 1.D0
  DOUBLE PRECISION, PARAMETER  :: epsw = 1.D0
  DOUBLE PRECISION, PARAMETER  :: epsn = 1.D-15

  !     Quadrature initialization
  DOUBLE PRECISION, PARAMETER  :: nleft = 1.D-4
  DOUBLE PRECISION, PARAMETER  :: nright = 1.D-4
  DOUBLE PRECISION, PARAMETER  :: uleft = 1.D-4
  DOUBLE PRECISION, PARAMETER  :: uright = 1.D-4
  DOUBLE PRECISION, PARAMETER  :: vleft = 1.D-4
  DOUBLE PRECISION, PARAMETER  :: vright = 1.D-4
  DOUBLE PRECISION, PARAMETER  :: wleft = 1.D-4
  DOUBLE PRECISION, PARAMETER  :: wright = 1.D-4

  !     Granular temperature lower limit
  DOUBLE PRECISION, PARAMETER  :: MINIMUM_THETA = 1.0D-6
  DOUBLE PRECISION, PARAMETER  :: MAXIMUM_SIGMA = 1.0D-6

END MODULE qmomk_parameters
