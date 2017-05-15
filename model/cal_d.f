
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE CAL_D calculates shear velocity
!

      SUBROUTINE CAL_D(V_sh)

      USE param
      USE param1
      USE parallel
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_s
      USE bc
      USE vshear
      USE compar
      USE fun_avg
      USE functions

      DOUBLE PRECISION V_sh,dis
!     DOUBLE PRECISION xdist(IMAX2,JMAX2)
!     DOUBLE PRECISION xdist3(IMAX2,JMAX2,KMAX2),cnter3(IMAX2,JMAX2,KMAX2)

!//SP Note the rational behind using global direction in I direction. This is to ensure
!     correctness in the way distances are calculated in the serial version and also
!     to give the capability to perform calculations over additional ghost layers in
!     to avoid if checks and communication
      DOUBLE PRECISION xdist(IMIN3:IMAX3,JSTART3:JEND3)
      DOUBLE PRECISION xdist3(IMIN3:IMAX3,JSTART3:JEND3,KSTART3:KEND3)

      INTEGER IJK,I1,J1,K1,I,J,K

      IF (NO_K) THEN

!calculate distances
      DO  J1= JSTART3, JEND3
        xdist(1,J1)=1d0/(ODX(1))
        if(imin3.ne.imin2) xdist(IMIN3,J1)=-1d0/(ODX(IMIN3))
        DO  I1 = 2, IMAX3
        xdist(I1,J1)=1d0/(ODX(I1))+xdist((I1-1),J1)
        END DO
      END DO

      DO  IJK= ijkstart3, ijkend3
          I = I_OF(IJK)
          J = J_OF(IJK)

        dis=xdist(I,J)

!shear velocity alligned u momentum cells

        VSHE(IJK)=V_sh-2d0*(1d0-(dis-(1d0/ODX(I))/xlength))*V_sh

!shear velocity alligned with scalar, v, w momentum cells

        VSH(IJK)=V_sh-2d0*(1d0-(dis-(1d0/ODX(I))/xlength))*V_sh&
        -2d0*(1d0/(2d0*ODX(I)*xlength))*V_sh

      END DO


      ELSE

!calculate distances
      DO K1=KSTART3,KEND3
        DO J1= JSTART3, JEND3
          IF (DEAD_CELL_AT(1,J1,K1)) CYCLE  ! skip dead cells
          xdist3(1,J1,K1)=1d0/(ODX(1))
          if(imin3.ne.imin2) xdist(IMIN3,J1)=-1d0/(ODX(IMIN3))
          DO  I1 = 2, IMAX3
           IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
           xdist3(I1,J1,K1)=1d0/(ODX(I1))+xdist3((I1-1),J1,K1)
          END DO
        END DO
      END DO


      DO  IJK= ijkstart3, ijkend3
          I = I_OF(IJK)
          J = J_OF(IJK)
          K = K_OF(IJK)

          dis=xdist3(I,J,K)

!shear velocity alligned u momentum cells

        VSHE(IJK)=V_sh-2d0*(1d0-(dis-(1d0/ODX(I))/xlength))*V_sh

!shear velocity alligned with scalar, v, w momentum cells

        VSH(IJK)=V_sh-2d0*(1d0-(dis-(1d0/ODX(I))/xlength))*V_sh&
        -2d0*(1d0/(2d0*ODX(I)*xlength))*V_sh

      END DO

      END IF
      RETURN
      END SUBROUTINE CAL_D

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 350 1206 change do loop limits: 1,kmax2->kstart3,kend3

