!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  subroutine name: pressure(s,alpha,ni,n,mu,sigma,chi,T,Ti,p)
!
!  author:  C. Hrenya, Jan 2009
!
!  Purpose: find pressure according to GHD polydisperse KT
!
!  Literature/References:
!     C. Hrenya handwritten notes & Garzo, Hrenya, Dufty papers (PRE, 2007)
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine pressure(s,alpha,ni,n,mu,sigma,chi,T,Ti,p)
      Implicit NONE

      integer s

      double precision alpha(s,s),ni(s),n,mu(s,s),sigma(s,s), &
                       chi(s,s),T,Ti(s),p

      integer i,j
      double precision pkin, pcol

      double precision pi
      parameter (pi=3.14159265458979323846d0)

      pkin = n*T                   !kinetic pressure
      pcol = 0.d0
      do i=1,s
         do j=1,s
           pcol = pcol + mu(j,i)*(1.d0+alpha(i,j))*sigma(i,j)**3* &
                            chi(i,j)*ni(i)*ni(j)*Ti(i)
         enddo
      enddo
      pcol = 2.d0*pi/3.d0*pcol     !collisional pressure
      p = pkin + pcol              !total pressure (p 9 CMH notes)

      return
      end subroutine pressure

