!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_RATES                                              !
!                                                                      !
!  Purpose: Silane Pyrolysis                                           !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jul-13  !
!                                                                      !
!  REF: B.Caussat, M.Hemati, and J.P.Couderc. Silicon deposition from  !
!       silane or disilane in a fluidized bed - Part 1: Experimental   !
!       study." Chemical Engineering Science,V ol. 50, No. 22,         !
!       pp. 3615-3624, 1995.                                           !
!                                                                      !
!  REF: T.Furusawa, T.Kojima, and H.Hiroha. Chemical vapor deposition  !
!       and homogeneous nucleation in monosilane pyrolysis within      !
!       interparticle spaces - Application of fines formation analysis !
!       to fluidized bed CVD. Chemical Engineering Science, Vol. 43,   !
!       No. 8, pp. 2037-2042, 1988.                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_RATES(IJK, RATES)

      USE param
      USE param1
      USE parallel
      USE fldvar
      USE rxns
      USE energy
      USE geometry
      USE run
      USE indices
      USE physprop
      USE constant
      USE funits
      USE compar
      USE sendrecv
      USE toleranc
      USE fun_avg
      USE usr
      USE functions

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJK

      DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(OUT) :: RATES

! Reaction specific variables:
!`````````````````````````````````````````````````````````````````````//
! Bounded phase temperatures (K)
      DOUBLE PRECISION xTg   ! gas phase
      DOUBLE PRECISION xTs   ! solids phase
! Gas phase pressure in KPa
      DOUBLE PRECISION Pg_KPa
! Gas phase concentrations (mol/cm^3)
      DOUBLE PRECISION c_SiH4, c_SiH2, c_Si2H6, c_H2
! Partial pressures (KPa)
      DOUBLE PRECISION p_SiH4, p_SiH2, p_H2
! Total solids surface area per volume (1/cm)
      DOUBLE PRECISION Sa
! Diffusion coefficient for SiH2 (cm^2/sec)
      DOUBLE PRECISION DIFF_SiH2
! SiH2 gas constant (KPa.cm^3)/(g-SiH2.K)
      DOUBLE PRECISION R_SiH2
! Calculate film diffusion resistance. (g-SiH2/(KPa.sec.cm^2)
      DOUBLE PRECISION K_f
! Apparent first order heterogeneous reation rate constant  (cm/sec)
      DOUBLE PRECISION k_so
! Constant of hydroggen inhibition of heterogeneous reaction (1/KPa)
      DOUBLE PRECISION K_H
! Constant of silane inhibition of heterogeneous reaction (1/KPa)
      DOUBLE PRECISION Ks
! Temporary values for reversible reactions
      DOUBLE PRECISION FWD, RVS, NET


! Prevent reactions when mass fraction falls below c_Limiter.
      DOUBLE PRECISION, parameter :: c_Limiter = 1.0d-7

!`````````````````````````````````````````````````````````````````````//
      INCLUDE 'species.inc'
      INCLUDE 'usrnlst.inc'
!`````````````````````````````````````````````````````````````````````//

! Calculate the bounded temperatures.
      xTg = max(min(T_g(IJK),   TMAX), 300.0d0)
      xTs = max(min(T_s(IJK,1), TMAX), 300.0d0)

! Convert gas pressure to KPa
      Pg_KPa = 1.0d-4 * P_g(IJK)  ! dyne/cm^2 --> KPa

! Initialize partial pressures and concentrations.
      c_SiH4  = ZERO;   p_SiH4 = ZERO;
      c_SiH2  = ZERO;   p_SiH2 = ZERO;
      c_H2    = ZERO;   p_H2   = ZERO;
      c_Si2H6 = ZERO;

! Calculate the partial pressure and molar concentration of SiH4
      IF(X_g(IJK,SiH4) > c_Limiter) then
         p_SiH4 = Pg_KPa * X_g(IJK,SiH4) * (MW_MIX_g(IJK) / MW_g(SiH4))
         c_SiH4 = RO_g(IJK) * X_g(IJK,SiH4) / Mw_g(SiH4)
      ENDIF

! Calculate the partial pressure and molar concentration of SiH2
      IF(X_g(IJK,SiH2) > c_Limiter) then
         p_SiH2 = Pg_KPa * X_g(IJK,SiH2) * (MW_MIX_g(IJK) / MW_g(SiH2))
         c_SiH2 = RO_g(IJK)*X_g(IJK,SiH2)/Mw_g(SiH2)
      ENDIF

! Calculate the molar concentration of Si2H6
      IF(X_g(IJK,Si2H6) > c_Limiter) then
         c_Si2H6 = RO_g(IJK)*X_g(IJK,Si2H6)/Mw_g(Si2H6)
      ENDIF

! Calculate the partial pressure and molar concentration of H2
      IF(X_g(IJK,H2) > c_Limiter) then
         p_H2 = Pg_KPa * X_g(IJK,H2) * (MW_MIX_g(IJK) / MW_g(H2))
         c_H2 = RO_g(IJK) * X_g(IJK,H2) / Mw_g(H2)
      ENDIF

!**********************************************************************!
!                         Homogeneous Reactions                        !
!**********************************************************************!

! SiH4 <--> SiH2 + H2  (g-mole/cm^3.s)
!---------------------------------------------------------------------//
! REF: Table 2, Caussat, Hemati, and Couderc (1995)
      FWD = 1.08d13 * exp(-2.56d4/xTg) * EP_g(IJK) * c_SiH4
      RVS = 1.26d12 * exp(-2.52d3/xTg) * EP_g(IJK) * c_SiH2 * c_H2

      NET = FWD - RVS
      IF(NET >= ZERO) THEN
         RATES(RX1F) =  NET
         RATES(RX1R) =  ZERO
      ELSE
         RATES(RX1F) =  ZERO
         RATES(RX1R) = -NET
      ENDIF

! Si2H6 <--> SiH4 + SiH2 (g-mole/cm^3.s)
!---------------------------------------------------------------------//
! REF: Table 2, Caussat, Hemati, and Couderc (1995)
      FWD = 5.62d15 * exp(-2.63d4/xTg) * EP_g(IJK) * c_Si2H6
      RVS = 1.00d14 * exp(-2.27d3/xTg) * EP_g(IJK) * c_SiH4 * c_SiH2

      NET = FWD - RVS
      IF(NET >= ZERO) THEN
         RATES(RX2F) =  NET
         RATES(RX2R) =  ZERO
      ELSE
         RATES(RX2F) =  ZERO
         RATES(RX2R) = -NET
      ENDIF


!**********************************************************************!
!                        Heterogeneous Reactions                       !
!**********************************************************************!
      IF(EP_s(IJK,1) > ZERO_EP_s) THEN

! Total surface area of solids in computational cell (cm^2/cm^3)
         Sa = 6.0d0 * EP_s(IJK,1) / D_p0(1)   ! (1/cm)

! SiH4 --> Si + 2H2  (g-mole/cm^3.s)
!---------------------------------------------------------------------//
! Ref. Furusawa, Kojima, and Hiroha (1988)
! Apparent first order heterogeneous reation rate constant  (cm/sec)
         k_so = 2.15d10 * exp(-2.3016d4/xTg)
! Constant of hydroggen inhibition of heterogeneous reaction (1/KPa)
         K_H  = 3.4d-2
! Constant of silane inhibition of heterogeneous reaction (1/KPa)
         KS  = 7.6d-3*exp(3.96d3/xTg)
! Reaction rate:
         RATES(RX3)  = (Sa*k_so*c_SiH4) / (ONE + K_H*p_H2 + KS*p_SiH4)


! SiH2 --> Si + H2  (g-mole/cm^3.s)
!---------------------------------------------------------------------//
! REF: Caussat, Hemati, and Couderc (1995)
! The reference states that 'all authors' consider their heterogeneous
! decomposition rate as infinite. This reaction rate is based on film
! layer resistance.

! Diffusion coefficient.
! This is based on Chapman-Enskog theory for diffusion in gases at low
! density. Approximately self-diffusivity of N2.
         DIFF_SiH2 = 3.45d-5*(xTg**1.5)  ! (cm^2/sec)
! SiH2 Gas constant. (Pa.m^3)/(kmol.K) == (KPa.cm^3)/(mol.K)
         R_SiH2 = (8.3147295d3/MW_g(SiH2)) ! (KPa.cm^3)/(g-SiH2.K)
! Calculate film diffusion resistance. (g-SiH2/(KPa.sec.cm^2)
         K_f = DIFF_SiH2 * N_sh(IJK,1) / (D_p0(1) * R_SiH2 * xTg)
! Calculate the reaction rate.
         RATES(RX4)  =  Sa * K_f * (p_SiH2 / MW_g(SiH2))
      ENDIF

      RETURN

      END SUBROUTINE USR_RATES
