!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: solids_pressure                                        C
!  Purpose: To compute solids pressure and its inverse                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-FEB-93  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:  1                                                 C
!  Purpose:  allow SI units                                            C
!  Author: S. Dartevelle                              Date: 01-Jul-02  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Comments:                                                           C
!     See S_pr1.inc for the parameters                                 C
!     to_SI is a constant to change from CGS to SI (Ba --> Pa)         C
!        see set_constants.f and constant_mod.f                        C
!        if CGS: to_SI=ONE and                                         C
!        if SI : to_SI=0.1                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE solids_pressure

  USE constant
  USE param1

!     parameters for P_s = a_ps*(EP_star - EP_g)**b_ps
  DOUBLE PRECISION, PARAMETER          :: a_ps = 1D25
  INTEGER, PARAMETER                   :: b_ps = 10

!     coefficient in Jackson's model (dyne/cm^2)
!      DOUBLE PRECISION, PARAMETER          :: A_ps_jackson = 0.5D0*100D0*2D0*981D0*0.4D0
!     voidage of random close-packed bed
!      DOUBLE PRECISION, PARAMETER          :: EP_g_cp = 0.35D0

CONTAINS

! Solids pressure in plastic-flow stress formulation. MFIX default model

! function P_s(EP_g,Ep_star)
  DOUBLE PRECISION FUNCTION Neg_H(XXX,YYY)
    IMPLICIT NONE
    DOUBLE PRECISION :: XXX, YYY
    Neg_H = to_SI*a_ps * (MAX(ZERO, (YYY - XXX )))**b_ps
  END FUNCTION Neg_H

! inverse of Neg_H.  function EP_g(P_s)
  DOUBLE PRECISION FUNCTION INV_H(XXX,YYY)
    IMPLICIT NONE
    DOUBLE PRECISION :: XXX, YYY
    INV_H = YYY - ( XXX/(to_SI*a_ps) )**(ONE/dble(b_ps))
  END FUNCTION INV_H

! Differentiate P_s w.r.t. EP_s.  function dP_s/dEP_s (EP_s)
  DOUBLE PRECISION FUNCTION dPodEP_s(XXX,YYY)
    IMPLICIT NONE
    DOUBLE PRECISION :: XXX, YYY
    dPodEP_s = to_SI*a_ps * dble(b_ps)*(MAX(ZERO, (YYY - (ONE - XXX) )))**(b_ps-1)
  END FUNCTION dPodEP_s

! Solids pressure in plastic-flow stress formulation. Jackson's model
! function P_s(EP_g,Ep_star)
      !Neg_H(XXX,YYY) = to_SI*a_ps_jackson * (MAX(ZERO, (YYY - XXX )/(XXX - EP_g_cp)))

! inverse of Neg_H.  function EP_g(P_s)
      !INV_H(XXX,YYY) = (to_SI*a_ps_jackson * YYY + XXX*EP_g_cp)/(to_SI*a_ps_jackson + XXX)

! Differentiate P_s w.r.t. EP_s.  function dP_s/dEP_s (EP_s)
      !dPodEP_s(XXX,YYY) = to_SI*a_ps_jackson * (YYY - EP_g_cp)/(XXX - EP_g_cp)**2

END MODULE solids_pressure
