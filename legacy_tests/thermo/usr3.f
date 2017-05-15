!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR3                                                   !
!  Author: J. Musser                                  Date: 31-MAR-14  !
!                                                                      !
!  Purpose: Write out the heat of reactions calculated in RATES0.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR3

      use energy
      use rxns
      use usr

      IMPLICIT NONE

! Error index
      INTEGER :: IER
! fluid cell index and loop counter
      INTEGER :: IJK, LC
! total heat of reaction.
      DOUBLE PRECISION :: tHORg, tHORs(1:3)
! File unit
      INTEGER, parameter :: lUnit = 668


! Initialize the 'sum' variables for the gas and solids phase heates
! of reaction.
      tHORg = 0.0d0
      tHORs = 0.0d0

! Calculate user defined reaction rates.
      CALL RRATES0(IER)

! Open the POST file again to append the heat of reaction data.
      OPEN(lUnit,FILE='POST_Thermo.dat',access='APPEND')

! Write the table header.
      WRITE(lUnit,1100) 'Reaction','HORg', 'HORs1', 'HORs2', 'HORs3'

! This case has ten fulid cells and ten reactions. Each fluid cell has
! a reaction rate of "ONE" for a single reaction and "ZERO" for all
! other reactions. Note that the first 49 IJK values are for ghost
! cells and therefore 49 is added to the reaction loop counter to
! map between the rates array and the fluid cells.
      DO LC=1, NO_OF_RXNS
         IJK = LC + 49

! Write out the calculated heat of reaction for reaction LC.
         WRITE(lUnit,1200) Reaction(LC)%Name(1:18), &
            HOR_g(IJK), HOR_s(IJK,1:3)
         tHORg = tHORg + HOR_g(IJK)
         tHORs = tHORs + HOR_s(IJK,1:3)
      ENDDO

! Close out the file.
      WRITE(lUnit,1300)'Total', tHORg, tHORs(1:3)
      CLOSE(lUnit)

      RETURN

! Formatted output statements.
 1100 FORMAT(3X,A,17X,A,3(11X,A))
 1200 FORMAT(3X,A,4(3X,G12.3))
 1300 FORMAT(90('-'),/3X,A,16X,G12.3,3(3X,G12.3))

      END SUBROUTINE USR3
