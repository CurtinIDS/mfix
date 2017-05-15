!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_DQMOM                                              !
!  Author: rong fan                                   Date: dd-mmm-yy  !
!                                                                      !
!  Purpose: This routine is called from the time loop and is           !
!           user-definable tosolve teh population equation             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_DQMOM

      USE run, only: CALL_DQMOM, TIME, DT
      USE scalars, only: NSCALAR
      USE scalars, only: IJK_INDEX, YSTART
      USE fldvar, only: ROP_s, RO_s, SCALAR
      use param, only: dim_scalar2

      USE geometry
      USE indices
      USE compar
      USE usr
      USE functions

      IMPLICIT NONE


      INTEGER :: IJK
! phase index
      double precision t1,t2
! beginning time and stop time
      double precision eps
! error
      double precision h1,hmin
      integer nok, nbad,I,K
      double precision max1, min1

      IF(.NOT.CALL_DQMOM) RETURN

!  Insert user-defined code here
      IF(time<= 1E-15) THEN
         t1 =time
         t2= time
      ELSE
         t1= time-dt
         t2= time
      ENDIF


      eps=1.0E-3
      h1=1.0E-4
      hmin=0
      nok=0
      nbad=0
!      n=2*Nscalar



      DO IJK = ijkstart3, ijkend3

         IF(.NOT.FLUID_AT(IJK)) CYCLE

         DO I=1,Nscalar
            ystart(I)=ROP_s(IJK,I)/RO_S(IJK,I)
         ENDDO

         DO I=Nscalar+1,2*Nscalar
            ystart(I)=Scalar(IJK,I-Nscalar)
         ENDDO

         IJK_INDEX=IJK
         max1=ystart(1)

         DO K=2,Nscalar
            max1=MAX(max1,ystart(K))
         ENDDO

         min1=ystart(1)

         DO K=2,Nscalar
            min1=MIN(min1,ystart(K))
         ENDDO


         IF(max1>1.0e-3) THEN
            call odeint(ystart,DIM_Scalar2,t1,t2,eps,h1, hmin,nok,nbad)
         ENDIF

         DO I=1,Nscalar
            ROP_s(IJK,I) = ystart(I)*RO_S(IJK,I)
            Scalar(IJK,I) = ystart(I+Nscalar)
         ENDDO

      ENDDO

      RETURN
      END SUBROUTINE USR_DQMOM
