!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_CONDUCTION                                         !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 25-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!  REF: Batchelor and O'Brien, "Thermal or electrical conduction       !
!       through a granular material," Proceedings of the Royal Society !
!       of London. Series A, Mathematical and Physical Sciences,       !
!       Vol. 355, no 1682, pp. 313-333, July 1977.                     !
!                                                                      !
!  REF: Rong and Horio, "DEM simulation of char combustion in a        !
!       fluidized bed," in Second International Conference on CFD in   !
!       the materials and Process Industries, Melbourne, 1999, pp.     !
!       65-70.                                                         !
!                                                                      !
!  REF: Xavier and Davidson, "Heat transfer to surfaces immersed in    !
!       fluidised beds, particularly tub arrays," in Fluidization:     !
!       Proceedings of the Second Engineering Foundation Conference,   !
!       Cambridge, England, 2-6 April, 1978, pp. 333-338.              !
!                                                                      !
!  REF: Zhou, Yu, and Zulli, "Particle scale study of heat transfer in !
!       packed and bubbling fluidized beds," AIChE Journal, Vol. 55,   !
!       no 4, pp 868-884, 2009.                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE DES_THERMO_COND
 

      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_Qw_cond
      LOGICAL :: DO_AREA_CORRECTION
      CONTAINS



      FUNCTION DES_CONDUCTION(I, J, CENTER_DIST, iM, iIJK)

      use constant
      use des_thermo
      use discretelement
      use funits
      use param1, only: zero, UNDEFINED
      use physprop
      use run, only: UNITS
      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
! Index of particle being looped over
      INTEGER, INTENT(IN) :: I,J
! Index of fluid cell containing particle I
      INTEGER, INTENT(IN) :: iIJK
! Solid phase indices for the given particles
      INTEGER, INTENT(IN) :: iM
! Distance between the centers of particle I and particle J (magnitude)
      DOUBLE PRECISION, INTENT(IN) :: CENTER_DIST
      DOUBLE PRECISION :: DES_CONDUCTION

! Local variables
!---------------------------------------------------------------------//
! Solid phase indices for the given particles
      INTEGER jM
! Radius of smaller particle
      DOUBLE PRECISION MIN_RAD
! Radius of larger particle
      DOUBLE PRECISION MAX_RAD
! Rate of particle-particle conduction
      DOUBLE PRECISION Q_pp
! Rate of particle-fluid-particle conduction
      DOUBLE PRECISION Q_pfp
! Outer radius of region delineating particle-fluid-particle conduction
      DOUBLE PRECISION RD_OUT
! Inner radius of region delineating particle-fluid-particle conduction
      DOUBLE PRECISION RD_IN
! The radius of the fluid lens containing the larger particle
      DOUBLE PRECISION LENS_RAD
! Temperature difference between two particles
      DOUBLE PRECISION DeltaTp
! Effective thermal conductivity
      DOUBLE PRECISION lK_eff
! Radius of contact area
      DOUBLE PRECISION lRadius
! Particle overlap (simulated and corrected)
      DOUBLE PRECISION OLAP_Sim, OLAP_actual
! Corrected center distance
      DOUBLE PRECISION CENTER_DIST_CORR
! Gas thermal conductivity
      DOUBLE PRECISION :: K_gas
! Overlap distance between each particle and contact plane
      DOUBLE PRECISION :: OLAP_1, OLAP_2
! Effective Lens for particles 1 and 2 (req'd for analytic calculation)
      DOUBLE PRECISION :: RLENS_1, RLENS_2
! Effective minimum conduction distance for particles 1 and 2
      DOUBLE PRECISION :: S_1, S_2
! Effective thermal conductance for particles 1 and 2
      DOUBLE PRECISION :: H_1, H_2
! Analytic thermal conductance between particles
      DOUBLE PRECISION :: H
! Dummy variable for calculations 
      DOUBLE PRECISION :: RATIO
! Functions
!---------------------------------------------------------------------//
 !     DOUBLE PRECISION :: EVAL_H_PFP

! Identify the solid phases of the neighbor particle
         jM = PIJK(J,5)

! Determine the radius of the larger and smaller particle
         MIN_RAD = MIN(DES_RADIUS(I), DES_RADIUS(J))
         MAX_RAD = MAX(DES_RADIUS(I), DES_RADIUS(J))

         DeltaTp = DES_T_s(J) - DES_T_s(I)

! Calculate the particle-particle conduction
! REF: Batchelor and O'Brien, 1977 (MODIFIED)
!---------------------------------------------------------------------//
         Q_pp = 0.0

! Initialize corrected center-center distance
         CENTER_DIST_CORR = CENTER_DIST
         IF(CENTER_DIST < (MAX_RAD + MIN_RAD)) THEN
! Correct particle overlap to account for artificial softening
            OLAP_Sim = MAX_RAD + MIN_RAD - CENTER_DIST
            IF(DO_AREA_CORRECTION)THEN
               IF(Im.ge.Jm)THEN
                  OLAP_Actual = CORRECT_OLAP(OLAP_Sim,Jm,Im)
               ELSE
                  OLAP_Actual = CORRECT_OLAP(OLAP_Sim,Im,Jm)
               ENDIF
               CENTER_DIST_CORR = MAX_RAD + MIN_RAD - OLAP_Actual
            ENDIF
! Effective thermal conductivity
            lK_eff = K_eff(K_s0(iM),K_s0(jM))
! Effective contact area's radius
            lRadius = RADIUS(MAX_RAD, MIN_RAD)
! Compute time-correction term (not done yet)
            ! CALL CALC_TIME_CORRECTION (.... )
! Inter-particle heat transfer
            Q_pp = 2.0d0 * lK_eff * lRadius * DeltaTp

! Assign the inter-particle heat transfer to both particles.
         ENDIF

!!         IF(.TRUE.)THEN ! NEW WAY
!---------------------------------------------------------------------//
! Calculate the particle-fluid-particle conduction
! REF: Rong and Horio, 1999 (MODIFIED)
!---------------------------------------------------------------------//
         LENS_RAD = MAX_RAD * (1.0D0 + FLPC)
         OLAP_Actual = MAX_RAD + MIN_RAD - CENTER_DIST_CORR
! check to see if particles are within lens thickness from eachother
! a negative overlap indicates particles are separated 
         IF(OLAP_actual > (-MAX_RAD*FLPC))THEN
            RD_OUT = RADIUS(LENS_RAD, MIN_RAD)
            ! get overlap between each particle and contact plane
            RATIO = OLAP_actual / MAX_RAD
            OLAP_1 = MAX_RAD*MAX_RAD*(RATIO-0.5D0*RATIO*RATIO)/CENTER_DIST_CORR
            OLAP_2 = OLAP_actual - OLAP_1
            ! get effective minimum conduction distance for each particle
            RATIO = (OLAP_actual + DES_MIN_COND_DIST) / MAX_RAD
            S_1 = (MAX_RAD*MAX_RAD*(RATIO-0.5D0*RATIO*RATIO)/&
            &(CENTER_DIST_CORR-DES_MIN_COND_DIST)) - OLAP_1
            S_2 = DES_MIN_COND_DIST - S_1
            ! get effective lens radius for particle-contact plane
            RLENS_1 = sqrt(RD_OUT*RD_OUT+(MIN_RAD-OLAP_1)**2.0)
            RLENS_2 = sqrt(RD_OUT*RD_OUT+(MAX_RAD-OLAP_2)**2.0)
            

! GET GAS THERMAL CONDUCTIVITY (default calculation is done in calc_k_g and is for air)
            if(k_g0.eq.UNDEFINED)then
                  ! Compute gas conductivity as is done in calc_k_g
                  ! But use average of particle and wall temperature to be gas temperature
               K_Gas = 6.02D-5*SQRT(0.5d0*(DES_T_S(I)+DES_T_S(J))/300.D0) ! cal/(s.cm.K)
                  ! 1 cal = 4.183925D0 J
               IF (UNITS == 'SI') K_Gas = 418.3925D0*K_Gas !J/s.m.K
            else
               K_Gas=k_g0
            endif
   
            H_1 = EVAL_H_PFP(RLENS_1, S_1, OLAP_1,MIN_RAD)*MIN_RAD*k_gas
            H_2 = EVAL_H_PFP(RLENS_2, S_2, OLAP_2,MAX_RAD)*MAX_RAD*k_gas
            IF(H_1.eq.ZERO.OR.H_2.eq.ZERO)THEN
               H = ZERO
            ELSE
               H = H_1*H_2/(H_1+H_2)
            ENDIF
     
            Q_pfp = H *DeltaTp
! Particle-fluid-particle is analytically computed using conductance for each
! particle to the contact plane.  The effective lens radius and minimum conduction
! distance are first calculated for each particle-fluid-contact_plane conduction. 
            
         ELSE
            Q_pfp = ZERO
         ENDIF
         DES_CONDUCTION = Q_pp + Q_pfp
               
            
 !!        ENDIF ! NEW WAY
!-------------------------------------------------------------------------------------
! OLD WAY
!-------------------------------------------------------------------------------------
!!         IF(.FALSE.)THEN
! Calculate the particle-fluid-particle conduction
! REF: Rong and Horio, 1999 (MODIFIED)
!---------------------------------------------------------------------//
! Calculate the radius of the fluid lens surrounding the larger particle
! Default FLPC = 0.2
!!            LENS_RAD = MAX_RAD * (1.0D0 + FLPC)

! Calculate the outer radial distance of the region for particle-fluid-
! particle heat conduction.
!!            RD_OUT = RADIUS( LENS_RAD, MIN_RAD)

! If the value returned is less than zero, then the fluid lens
! surrounding the larger particle does not intersect with the surface
! of the smaller particle. In this case, particle-fluild-particle
! conduction does not occur.
!!            Q_pfp = 0.0
!!            IF(RD_OUT .GT. ZERO)THEN
! Calculate the distance from the line connecting the particles' centers
! to the point of contact between the two particles. This value is
! zero if the particles are not touching and is the radius of the
! shared contact area otherwise.
!!               RD_IN = ZERO
!!               IF(CENTER_DIST_CORR < (MAX_RAD + MIN_RAD) ) &
!!               RD_IN = RADIUS(MAX_RAD, MIN_RAD)
! Calculate the rate of heat transfer between the particles through the
! fluid using adaptive Simpson's rule to manage the integral.
!!               Q_pfp = K_g(iIJK) * DeltaTp * &
!!               ADPT_SIMPSON(RD_IN,RD_OUT)

!!            ENDIF               ! RD_OUT > ZERO
!!         ENDIF ! OLD WAY

         ! Return total heat flux
!!         DES_CONDUCTION = Q_pp + Q_pfp
!! ------------------------------------------------------------------//
!  END OLD WAY
!! ------------------------------------------------------------------//

      RETURN

      CONTAINS

!......................................................................!
! Subroutine Procedure: FUNCTION RADIUS                                !
!                                                                      !
! Propose: Calculate the center line radius between two particles.     !
! Depending on what is given as input, this function calculates:       !
!   1) radius of contact area between two particles                    !
!   2) radius delineating particle-fluid-particle region.              !
!......................................................................!
      DOUBLE PRECISION FUNCTION RADIUS(R1, R2)

      IMPLICIT NONE

! Radius values
      DOUBLE PRECISION, INTENT(IN) :: R1, R2
! Distance between particle centers
      DOUBLE PRECISION BETA
! Generic value
      DOUBLE PRECISION VALUE

! Calculate
      VALUE = (R1**2 - R2**2 + CENTER_DIST_CORR**2)/(2.0d0 * R1 * CENTER_DIST_CORR)
! Check to ensure that VALUE is less than or equal to one. If VALUE is
! greater than one, the triangle inequality has been violated. Therefore
! there is no intersection between the fluid lens surrounding the larger
! particle and the surface of the smaller particle.
! Thus, there is no particle-fluid-particle heat transfer.
      IF( VALUE .GT. 1.0d0) THEN
         RADIUS = -1.0d0
      ELSE
! Calculate beta (Law of cosines)
         BETA = ACOS( VALUE )
! Calculate the radius
         RADIUS = R1 * SIN(BETA)
      ENDIF

      RETURN
      END FUNCTION RADIUS

!......................................................................!
! Subroutine Procedure: FUNCTION K_eff                                 !
!                                                                      !
! Propose: Calculate the effective thermal conductivity of two         !
! particles with conductivities K1 and K2                              !
!......................................................................!
      DOUBLE PRECISION FUNCTION K_eff(K1, K2)

      IMPLICIT NONE

! Thermal conductivities
      DOUBLE PRECISION, INTENT(IN) :: K1, K2

      K_eff = 2.0d0 * (K1*K2)/(K1 + K2)

      RETURN
      END FUNCTION K_eff

!``````````````````````````````````````````````````````````````````````!
! Function: F                                                          !
!                                                                      !
! Purpose: This function defines the region of particle-fluid-particle !
!          heat transfer. Note that this function develops a           !
!          singularity at the boundary of the contact region when two  !
!          particles are touching.  This singularity is resolved by    !
!          making the assumption that the surfaces of two particles    !
!          never directly touch (c.f. Rong and Horio, 1999). By        !
!          default, it is assumed that the particles are separated by  !
!          a minimum length of 4.0x10^(-10) meters. This value may be  !
!          modified by specifying a different value in the mfix.dat    !
!          file under the variable name DES_MIN_COND_DIST.             !
!                                                                      !
!          Note, the closer this value is to zero, the harder it will  !
!          be for the quadrature routine to converge.                  !
!......................................................................!
      DOUBLE PRECISION FUNCTION F(R)
        USE utilities
        USE des_thermo
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: R

      F = (2.0d0*Pi*R)/MAX((CENTER_DIST_CORR - SQRT(MAX_RAD**2-R**2) - &
         SQRT(MIN_RAD**2-R**2)), DES_MIN_COND_DIST)

      IF( mfix_isnan(F))PRINT*,'F IS NAN'

      END FUNCTION F

!``````````````````````````````````````````````````````````````````````!
! Function: ADPT_SIMPSON                                               !
!                                                                      !
! Purpose: This function applies adaptive Simpson's Rule to solve the  !
!          definite integral of the function F.                        !
!                                                                      !
!         *NOTE: This route will return the result of the method even  !
!          if convergence has not been achieved. If this occurs, a     !
!          message will be written to the log file to notify the user. !
!......................................................................!
      DOUBLE PRECISION FUNCTION ADPT_SIMPSON(A, B)

      IMPLICIT NONE

! Bounds of integration
      DOUBLE PRECISION, INTENT(IN) :: A, B
! Maximum recursive depth
      INTEGER, PARAMETER :: Kmax = 25
! Current depth of recursion
      INTEGER :: K
! Array storage for recursion
      DOUBLE PRECISION :: V(Kmax,6)
! Step size
      DOUBLE PRECISION :: H
! Simpson's Rule evaluated on the left and right intervals
      DOUBLE PRECISION :: lS, rS
! Value of the function F evaluate at the midpoints of the left and
! right intervals
      DOUBLE PRECISION :: rF, lF
! Local error bound
      DOUBLE PRECISION :: Err_BND
! Convergence bound
      DOUBLE PRECISION :: EPS = 10.0**(-2)

! Error indicating that an error has been specified to the user.
      LOGICAL, SAVE :: ADPT_SIMPSON_ERR = .FALSE.

! Initialize variables
      V(:,:) = 0.0d0   !
      ADPT_SIMPSON = 0.0d0   ! Integral value
      H = (B-A)/2.0d0 ! Dynamic interval length
! Calculate Simpson's Rule over the interval [A,B]
      lS = (F(A) + 4.0d0*F((A+B)/2.0d0) + F(B))*(H/3.0d0)
! Initialize the storage vector for the first pass
      V(1,:) = (/ A, H, F(A), F((A+B)/2.0d0), F(B), lS/)
      K = 1
! Enter recursion calculations
      DO
! Establish the new interval length
         H = V(K,2)/2.0d0
! Evaluate the function on the left interval
         lF = F(V(K,1) + H)
! Calculate Simpson's Rule over the left interval
         lS = (V(K,3) + 4.0d0*lF + V(K,4))*(H/3.0d0)
! Evaluate the function on the right interval
         rF = F(V(K,1) + 3.0d0*H)
! Calculate Simpson's Rule over the right interval
         rS = (V(K,4) + 4.0d0*rF + V(K,5))*(H/3.0d0)
! Evaluate the error bound
         Err_BND = (30.0d0*EPS*H)/(B-A)
! Check to see if conversion on the interval has been met
         IF( ABS(lS + rS - V(K,6)) .LT. Err_BND)THEN
! Update the integral value
            ADPT_SIMPSON = ADPT_SIMPSON + lS + rS + &
               (1.0d0/15.0d0)*(lS + rS - V(K,6))
! Decrement the recursion counter
            K = K - 1
! If K=0, then integration has been successfully completed over the
! entire interval [A,B].
            IF( K == 0) RETURN
         ELSEIF( (K .GE. Kmax) .OR. &
            (H == (B-A)*(1.0d0/2.0d0)**(Kmax+3))) THEN
! Flag that the method did not converge.
            IF(.NOT.ADPT_SIMPSON_ERR)THEN
               WRITE(*,1000)
               WRITE(UNIT_LOG,1000)
               ADPT_SIMPSON_ERR = .TRUE.
            ENDIF
! Update the integral value
            ADPT_SIMPSON = ADPT_SIMPSON + lS + rS + &
               (1.0d0/15.0d0)*(lS + rS - V(K,6))
! Decrement the recursion counter
            K = K - 1
         ELSE
! Refine the subintervals through recursive splitting of the intervals
            V(K+1,:) = (/V(K,1) + 2.0d0*H, H, V(K,4), rF, V(K,5), rS/)
            V(K,:) = (/V(K,1), H, V(K,3), lF, V(K,4), lS/)
! Increment level counter
            K = K+1
         ENDIF
      ENDDO

 1000 FORMAT(/1X,70('*'),/' From: DES_COND_EQ',/, ' Message: ',        &
         'Integration of the particle-fluid-particle equation did ',   &
         'not',/' converge! No definite bound can be placed on the ',  &
         'error.',/' Future convergence messages will be suppressed!', &
         /1X,70('*'))

      END FUNCTION ADPT_SIMPSON

     DOUBLE PRECISION FUNCTION EVAL_H_PFP(RLENS_dim,S,OLAP_dim,RP)
      USE CONSTANT
      USE PARAM1
     
      IMPLICIT NONE
      ! Note: Function inputs dimensional quantities
      DOUBLE PRECISION, intent(in) :: RLENS_dim, S, OLAP_dim, RP
      ! BELOW VARIABLES ARE NONDIMENSIONALIZED BY RP
      DOUBLE PRECISION :: RLENS, OLAP, KN
      DOUBLE PRECISION :: TERM1,TERM2,TERM3
      DOUBLE PRECISION :: Rout,Rkn
      DOUBLE PRECISION, PARAMETER :: TWO = 2.0D0
      
      RLENS = RLENS_dim/RP
      KN = S/RP
      OLAP = OLAP_dim/RP

      IF(-OLAP.ge.(RLENS-ONE))THEN
         Rout = ZERO
      ELSEIF(OLAP.le.TWO)THEN
         Rout = sqrt(RLENS**2-(ONE-OLAP)**2)
      ELSE
         WRITE(*,*)'ERROR: Extremeley excessive overlap (Olap > Diam)'
         WRITE(*,*)'OLAP = ',OLAP_dim,OLAP, RP
         WRITE(*,*)'RLENS_dim',RLENS_dim, RLENS
         write(*,*)'S, Kn', S,KN
         STOP
      ENDIF
      Rout=MIN(Rout,ONE)
      
      IF(OLAP.ge.ZERO)THEN
     !     Particle in contact (below code verified)
         TERM1 = PI*((ONE-OLAP)**2-(ONE-OLAP-KN)**2)/KN
         TERM2 = TWO*PI*(sqrt(ONE-Rout**2)-(ONE-OLAP-KN))
         TERM3 = TWO*PI*(ONE-OLAP)*log((ONE-OLAP-sqrt(ONE-Rout**2))/Kn)
         EVAL_H_PFP = TERM1+TERM2+TERM3
      ELSE
         IF(-OLAP.ge.KN)THEN
            Rkn = ZERO
         ELSE
            Rkn=sqrt(ONE-(ONE-OLAP-Kn)**2)
         ENDIF
    
         TERM1 = (Rkn**2/(TWO*KN))+sqrt(ONE-Rout**2)-sqrt(ONE-Rkn**2)
         TERM2 = (ONE-OLAP)*log((ONE-OLAP-sqrt(ONE-Rout**2))/(ONE-OLAP-sqrt(ONE-Rkn**2)))
         EVAL_H_PFP = TWO*PI*(TERM1+TERM2)
         
      ENDIF
      RETURN
      END FUNCTION EVAL_H_PFP

      END FUNCTION DES_CONDUCTION

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_CONDUCTION                                         !
!                                                                      !
!  Purpose: Compute conductive heat transfer with wall                 !
!                                                                      !
!  Author: A.Morris                                  Date: 07-Jan-16   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      FUNCTION DES_CONDUCTION_WALL(OLAP, K_sol, K_wall, K_gas, TWall, &
      TPart, Rpart, RLens, M)
      USE param1, only: zero
      USE DES_THERMO, only: des_min_cond_dist

      DOUBLE PRECISION :: DES_CONDUCTION_WALL 
      DOUBLE PRECISION, intent(in) :: OLAP, K_sol, K_wall, K_gas, TWall&
      &                               ,TPart,RPart, RLens
      DOUBLE PRECISION :: Rin, Rout
      DOUBLE PRECISION :: OLAP_ACTUAL, TAUC_CORRECTION
      DOUBLE PRECISION :: KEff
      DOUBLE PRECISION :: Q_pw, Q_pfw
      INTEGER :: M ! Solids phase index

      ! Initialize variables
      Q_pw = ZERO
      Q_pfw= ZERO
      OLAP_ACTUAL = OLAP

      ! Check to see if particle is in contact (P-W) conduction
      IF(OLAP.gt.ZERO)THEN
      ! Compute effective solids conductivity
         KEff = 2.0D0*K_sol*K_wall/(K_sol+K_wall)
      ! Correct for overlap using "area" correction terms
         IF(DO_AREA_CORRECTION)THEN
            OLAP_ACTUAL = CORRECT_OLAP(OLAP,M,-1) ! (-1 indicates wall)
         ENDIF
      ! Compute inner Rad
         Rin = sqrt(2.0D0*OLAP_ACTUAL*Rpart-OLAP_ACTUAL*OLAP_ACTUAL)
      ! Compute time correction term (commented out for now)
        ! CALL CALC_TIME_CORRECTION(TAUC_CORRECTION,M,-1)
      ! Compute heat transfer (particle-wall)
         Q_pw = 2.0D0*Rin*Keff*(TWall-TPart)
      ENDIF
      ! Compute heat transfer (particle-fluid-wall)
      Q_pfw=EVAL_H_PFW(RLens,DES_MIN_COND_DIST, OLAP_ACTUAL, RPart) * &
      &     K_gas*RPart*(TWall-TPart)

      DES_CONDUCTION_WALL=Q_pw+Q_pfw
     
      RETURN 
      END FUNCTION DES_CONDUCTION_WALL


!     THIS FUNCTION COMPUTES THE OVERLAP THAT WOULD OCCUR (for a given force)
!     IF ACTUAL MATERIAL PROPERTIES WERE USED. 
      FUNCTION CORRECT_OLAP(OLAP,M,L)
      use discretelement
      IMPLICIT NONE
      DOUBLE PRECISION :: CORRECT_OLAP
      DOUBLE PRECISION, INTENT (IN) :: OLAP
      INTEGER, INTENT (IN) :: M, L
      DOUBLE PRECISION :: KN_ACTUAL, KN_SIM
      ! L=-1 corresponds to wall
      IF(L.eq.-1)THEN
         ! WALL CONTACT
         IF(DES_COLL_MODEL_ENUM .EQ. HERTZIAN)THEN
            CORRECT_OLAP = (HERT_KWN(M)/HERT_KWN_ACTUAL(M))**(2.0D0/3.0D0)*OLAP
         ELSE
            CORRECT_OLAP = (KN_W*OLAP/HERT_KWN_ACTUAL(M))**(2.0D0/3.0D0)
         ENDIF
      ELSE
         ! PARTICLE-PARTICLE
         IF(DES_COLL_MODEL_ENUM .EQ. HERTZIAN)THEN
            CORRECT_OLAP = (HERT_KN(M,L)/HERT_KN_ACTUAL(M,L))**(2.0D0/3.0D0)*OLAP
         ELSE
            CORRECT_OLAP = (KN*OLAP/HERT_KN_ACTUAL(M,L))**(2.0D0/3.0D0)
         ENDIF
      ENDIF
      RETURN
      END FUNCTION CORRECT_OLAP


      DOUBLE PRECISION FUNCTION EVAL_H_PFW(RLENS_dim,S,OLAP_dim,RP)
      USE CONSTANT
      USE PARAM1
     
      IMPLICIT NONE
      ! Note: Function inputs dimensional quantities
      DOUBLE PRECISION, intent(in) :: RLENS_dim, S, OLAP_dim, RP
      ! BELOW VARIABLES ARE NONDIMENSIONALIZED BY RP
      DOUBLE PRECISION :: RLENS, OLAP, KN
      DOUBLE PRECISION :: TERM1,TERM2,TERM3
      DOUBLE PRECISION :: Rout,Rkn
      DOUBLE PRECISION, PARAMETER :: TWO = 2.0D0
      
      RLENS = RLENS_dim/RP
      KN = S/RP
      OLAP = OLAP_dim/RP

      IF(-OLAP.ge.(RLENS-ONE))THEN
         Rout = ZERO
      ELSEIF(OLAP.le.TWO)THEN
         Rout = sqrt(RLENS**2-(ONE-OLAP)**2)
      ELSE
         WRITE(*,*)'ERROR: Extremeley excessive overlap (Olap > Diam)'
         WRITE(*,*)'OLAP = ',OLAP_dim,OLAP, RP
         WRITE(*,*)'RLENS_dim',RLENS_dim, RLENS
         write(*,*)'S, Kn', S,KN
         STOP
      ENDIF
      Rout=MIN(Rout,ONE)
      
      IF(OLAP.ge.ZERO)THEN
     !     Particle in contact (below code verified)
         TERM1 = PI*((ONE-OLAP)**2-(ONE-OLAP-KN)**2)/KN
         TERM2 = TWO*PI*(sqrt(ONE-Rout**2)-(ONE-OLAP-KN))
         TERM3 = TWO*PI*(ONE-OLAP)*log((ONE-OLAP-sqrt(ONE-Rout**2))/Kn)
         EVAL_H_PFW = TERM1+TERM2+TERM3
      ELSE
         IF(-OLAP.ge.KN)THEN
            Rkn = ZERO
         ELSE
            Rkn=sqrt(ONE-(ONE-OLAP-Kn)**2)
         ENDIF
    
         TERM1 = (Rkn**2/(TWO*KN))+sqrt(ONE-Rout**2)-sqrt(ONE-Rkn**2)
         TERM2 = (ONE-OLAP)*log((ONE-OLAP-sqrt(ONE-Rout**2))/(ONE-OLAP-sqrt(ONE-Rkn**2)))
         EVAL_H_PFW = TWO*PI*(TERM1+TERM2)
         
      ENDIF
      RETURN
      END FUNCTION EVAL_H_PFW


      SUBROUTINE CALC_TIME_CORRECTION(time_corr, phaseI, phaseJ)
      use discretelement, only: tau_c_base_actual, tauw_c_base_actual
      use discretelement, only: tau_c_base_sim, tauw_c_base_sim
      use discretelement, only: des_coll_model_enum, HERTZIAN
      IMPLICIT NONE
      INTEGER, intent (in) :: phaseI, phaseJ
      DOUBLE PRECISION :: time_corr, tau_actual, tau_sim
      DOUBLE PRECISION :: vimp ! impact veloctiy
      INTEGER :: M, N
      vimp = 1.0D0
      time_corr = 1.0D0
      if (phaseJ == -1) then
         ! Wall contact
         M = phaseI
         IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
            time_corr = tauw_c_base_actual(M) / tauw_c_base_sim(M)
         ELSE
            vimp = 1.0D0
            if (vimp .le. 0.0D0)then
               time_corr = 1.0D0
            else
               time_corr = tauw_c_base_actual(M)*vimp**(-0.2D0) / tauw_c_base_sim(M)
            endif
         ENDIF

      ELSE
         ! particle-particle contact
         if (phaseI .le. phaseJ)then
            M = phaseI
            N = phaseJ
         else
            M = phaseJ
            N = phaseI
         endif
         
         IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
            time_corr = tau_c_base_actual(M,N) / tau_c_base_sim(M,N)
         ELSE
            vimp = 1.0D0
            if (vimp .le. 0.0D0)then
               time_corr = 1.0D0
            else
               time_corr = tau_c_base_actual(M,N)*vimp**(-0.2D0) / tau_c_base_sim(M,N)
            endif
         ENDIF
      ENDIF
      time_corr = time_corr **(2.0D0/3.0D0)
      ! TEMPORARY 
      time_corr = 1.0D0
      RETURN
      END SUBROUTINE CALC_TIME_CORRECTION


      
  END MODULE DES_THERMO_COND
