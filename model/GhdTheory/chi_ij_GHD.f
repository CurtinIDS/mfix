      subroutine chi_ij_GHD(s,i,j,sigmai,phi,ni,chi_ij)
      Implicit NONE

      integer s
      double precision pi
      parameter (pi=3.14159265458979323846d0)

      integer i,j
      double precision sigmai(s)
      double precision phi
      double precision ni(s)

      integer kk
      double precision beta
      double precision phi_g
      double precision sigma_ij
      double precision group1
      double precision chi_ij

      beta = 0.d0
      do kk=1,s
         beta = beta + pi/6.d0*ni(kk)*sigmai(kk)**2
      enddo
      phi_g = 1d0-phi
      sigma_ij = 0.5d0*(sigmai(i)+sigmai(j))
      group1 = (sigmai(i)*sigmai(j))/sigma_ij

      chi_ij = 1.d0/phi_g + 1.5d0*beta/phi_g**2*group1 &
               + 0.5d0*beta**2/phi_g**3*group1**2

      return
      end subroutine chi_ij_GHD
