!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine name: ghd_model
!
!  Purpose: find transport coefficients of GHD polydisperse KT
!           for known inputs.
!
!  Literature / References
!     C. Hrenya handnotes and Garzo, Hrenya, Dufty papers (PRE, 2007)
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE  GHD_MODEL (S, SIGMAI,IJK, alpha, MI, phii, T, Zeta0, &
                             zetau, Ti, P, Kappa, Eta, DT, DF, Lambda, &
                             Lij, Dij, Dq)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE drag
      Implicit NONE
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      double precision :: pi
      parameter (pi=3.14159265458979323846d0)
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      integer :: s                           !number of species
      double precision sigmai(s)             !diameter of each species
      integer :: ijk                         !fluid cell index
      double precision :: alpha(s,s)         !restitution coeff of each pair
      double precision :: mi(s)              !mass of each species
      double precision :: phii(s)            !solids volume fraction of each species
      double precision :: T                  !granular (mixture) temperature
      double precision :: zeta0              !zeta0 is cooling rate
      double precision :: zetau              !cooling rate transport coefficient (1st order)
      double precision :: Ti(s)              !species temperature
      double precision :: p                  !pressure
      double precision :: kappa              !bulk viscosity
      double precision :: eta                !shear viscosity
      double precision :: DT(s)              !thermal diffusivity
      double precision :: DF(s,s)            !mass mobility coefficient
      double precision :: lambda             !thermal conductivity
      double precision :: Lij(s,s)           !thermal mobility
      double precision :: Dij(s,s)           !ordinary diffusion
      double precision :: Dq(s,s)            !Dufour coefficient
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer i, j, k
      double precision ni(s)                 !number density of each species
      double precision phi                   !total solids fraction
      double precision n                     !total number density
      double precision rhoi(s)               !mi(i)*ni(i) of each species
      double precision rho                   !m*n
      double precision m                     !average mass
      double precision v0                    !commonly used quantities: v0, mu, sigma
      double precision mu(s,s)
      double precision sigma(s,s)
      double precision chi_ij
      double precision chi(s,s)              !radial distribution function at contact for each pair
      double precision group1(s,s), group2(s,s)
      double precision theta(s)
      double precision beta(s,s)             !defined on p11 CMH notes
      double precision Beta_tot, niTi, scale_fac
      double precision nu(s,s)               !commonly used quantity: p7 cmh notes
      double precision I_ilj(s,s,s)          !commonly used quantity: p6.1 CMH notes
      double precision dTl_dnj(s,s)          !partial derivative of Tl wrt nj: p 6 CMH notes
      double precision dzeta0_dnj(s)         !partial derivative of zeta 0 wrt nj: p 6 CMH notes
      double precision dchi0il_dnj(s,s,s)    !partial derivative of chi_il wrt nj: p 6 CMH notes
      double precision :: omega(s,s)         !commonly used quantity: p16 CMH notes
      double precision :: gammaij(s,s)       !commonly used quantity: p13 CMH notes
      double precision :: lambdai(s)         !commonly used quantity: p13 CMH notes
!-----------------------------------------------


! Initializing
      phi = 0.d0
      n   = 0.d0
      rho = 0.d0
      m   = 0.d0
      niTi = 0.d0
      Beta_tot = 0.d0
      do i=1,s
         phi = phi + phii(i)
         ni(i) = 6.d0*phii(i) / (pi*sigmai(i)**3)
         n = n+ ni(i)
         rhoi(i) = mi(i)*ni(i)
         rho = rho + rhoi(i)
         m = m + mi(i)
         Beta_tot = Beta_tot + F_GS(IJK,s)
      enddo
      do i=1,s
        if(n .eq. 0) then  ! do not do any calculation if total solids concentration is zero.
          Ti(:)      = T
          zeta0      = 0.d0
          zetau      = 0.d0
          p          = 0.d0
          kappa      = 0.d0
          eta        = 0.d0
          DT(:)      = 0.d0
          DF(:,:)    = 0.d0
          lambda     = 0.d0
          Lij(:,:)   = 0.d0
          Dij(:,:)   = 0.d0
          RETURN
        endif
      enddo
      m = m/dble(s)
      v0 = dsqrt(2.d0*T/m)

      do i=1,s
         do j=1,s
            mu(i,j) = mi(i)/(mi(i)+mi(j))
            sigma(i,j) = 0.5d0*(sigmai(i)+sigmai(j))
            call chi_ij_GHD(s,i,j,sigmai,phi,ni,chi_ij)
            chi(i,j) = chi_ij
         enddo
      enddo
!---------------------------------------------------------------------
! ZEROTH-ORDER COOLING RATE (zeta0)
!---------------------------------------------------------------------
! This subroutine solves the nonlinear algebraic equations for theta
      call cooling_rate(s,mi,ni,n,m,Ti,T,chi,sigmai,alpha,rhoi,theta)

! Ti and zeta0 are calculated from theta
      do i=1,s
         Ti(i) = mi(i)*T/(m*theta(i))
         niTi = niTi + ni(i)*Ti(i)
      enddo

      do j=1,s
           group1(1,j)  = chi(1,j)*ni(j)*mu(j,1)* &
                            sigma(1,j)**2*(1d0+alpha(1,j))
           group2(1,j)  = mu(j,1)/2.d0*(1d0+alpha(1,j))
      enddo
      zeta0 = 0d0;
      do k=1,s
         zeta0 = zeta0 + group1(1,k)*dsqrt((theta(1)+theta(k)) &
                       /(theta(1)*theta(k))) &
                * ((1.d0-group2(1,k)*(theta(1)+theta(k))/theta(k)))
      enddo
      zeta0 = 8.d0/3.d0*dsqrt(2.d0*pi*T/m)*zeta0

! this quantity is needed for determination of several other transport coefficients
      do i=1,s
         do j=1,s
            beta(i,j) = mu(i,j)*theta(j)-mu(j,i)*theta(i);  !p 11 CMH notes
         enddo
      enddo

!---------------------------------------------------------------------
! COOLING RATE TRANSPORT COEFFICIENT (zetaU)
!---------------------------------------------------------------------
      call cooling_rate_tc(s,mi,sigmai,alpha,ni,n,v0,mu,sigma,chi,T, &
                          zeta0,theta,Ti,zetau)
!---------------------------------------------------------------------
! PRESSURE (p)
!---------------------------------------------------------------------
      call pressure (s,alpha,ni,n,mu,sigma,chi,T,Ti,p)

!---------------------------------------------------------------------
! BULK VISCOSITY (kappa)
!---------------------------------------------------------------------
      call bulk_viscosity(s,mi,alpha,ni,v0,mu,sigma,chi,theta,kappa)

!---------------------------------------------------------------------
! SHEAR VISCOSITY (eta)
!---------------------------------------------------------------------
      call shear_viscosity(s,mi,sigmai,alpha,ni,v0,mu,sigma,chi, &
                           beta,zeta0,theta,Ti,kappa,eta)

!---------------------------------------------------------------------
! THERMAL DIFFUSION COEFFICIENT (DT) & nu
!---------------------------------------------------------------------
      call thermal_diffusivity(s,alpha,ni,mi,rho,v0,mu,sigma,chi, &
                               zeta0,theta,Ti,p,DT,nu)

!---------------------------------------------------------------------
! MASS MOBILITY COEFFICIENT (DF)
!---------------------------------------------------------------------
      call mass_mobility(s,mi,ni,rho,zeta0,theta,nu,DF)

!---------------------------------------------------------------------
! ORDINARY DIFFUSION (Dij)
!---------------------------------------------------------------------
      call ordinary_diff(s,mi,sigmai,alpha,phii,T,phi,ni,n,rhoi,rho, &
             m,mu,sigma,chi,zeta0,theta,Ti,DT,nu,Dij,I_ilj,dTl_dnj, &
             dzeta0_dnj,dchi0il_dnj)

!---------------------------------------------------------------------
! CONDUCTIVITY (lambda)
!---------------------------------------------------------------------
      call thermal_conductivity(s,mi,T,sigmai,alpha,ni,rho,v0,mu,&
                         sigma,chi,beta,zeta0,theta,Ti,DT,&
                         lambda,omega,gammaij,lambdai)

!---------------------------------------------------------------------
! THERMAL MOBILITY COEFFICIENT (Lij)
!---------------------------------------------------------------------
      call thermal_mobility(s,mi,alpha,ni,mu,sigma,chi,zeta0,&
                            theta,Ti,DF,gammaij,omega,Lij)

!---------------------------------------------------------------------
! DUFOUR COEFFICIENT (Dq)
!---------------------------------------------------------------------
      call dufour_coeff(s,mi,alpha,T,ni,rho,v0, &
             mu,sigma,chi,beta,zeta0,theta,Ti,Dij,lambdai,gammaij, &
             omega,I_ilj,dTl_dnj,dzeta0_dnj,dchi0il_dnj,Dq)


! SCALE ALL TRANSPORT COEFFICIENTS FOR THE LOW NUMBER DENSITY LIMIT
      scale_fac = 1d0+2d0*eta*Beta_tot/(rho*niTi)


      lambda = lambda/scale_fac
      eta = eta/scale_fac
      kappa = kappa/scale_fac
      do i=1,s
          DT(i) = DT(i)/scale_fac
        do j=1,s
          Lij(i,j) = Lij(i,j)/scale_fac
          Dq(i,j) = Dq(i,j)/scale_fac
          Dij(i,j) = Dij(i,j)/scale_fac

! do not scale mass mobility as fluxes need be computed correctly in dilute limit.
!          DF(i,j) = DF(i,j)/scale_fac
        enddo
      enddo

      RETURN
      END SUBROUTINE GHD_MODEL

