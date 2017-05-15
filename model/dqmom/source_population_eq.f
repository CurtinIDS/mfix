
      SUBROUTINE Source_population_eq(x,y,dydx)

      use param1, only: zero, small_number
      USE constant
      USE fldvar
      USE physprop
      USE rdf
      USE scalars

      IMPLICIT NONE

      double precision x,y(*),dydx(*)
      double precision K_v
      double precision m1,m2, dav,theta,c11
      double precision epstotal
      integer L,n,np,i,j,k,IJK

      k_v=PI/6
      epstotal=0.0
      theta=0.0


      DO I=1,Nscalar
      A(I)= 0.0
      omega(I)=y(I)/(k_v*(y(Nscalar+I)**3))
      ENDDO

      n=2*Nscalar
      np=2*Nscalar

      IJK=IJK_INDEX

! initinalize variables
      DO i=1,n
       DO j=1,n
         matrix_a(i,j)=ZERO
         matrix_b(i,j)=ZERO
         matrix_c(i,j)=ZERO
         inv_a(i,j)=ZERO
       ENDDO
      ENDDO

      DO i=1,n
       dydx(I)=ZERO
      ENDDO

      DO i=0,2*Nscalar-1
       S_bar(I)=ZERO
      ENDDO

! calculate BETA_A(I,J) and A(I)

      DO I=1, Nscalar
        epstotal=ROP_s(IJK,I)/RO_S(IJK,I)+epstotal
      ENDDO

       If(epstotal>=small_number) THEN
           Do J=1,Nscalar
            theta= theta+(ROP_s(IJk,J)/RO_S(IJK,J))*Theta_m(IJK,J)
           ENDDO
            theta= theta/epstotal
        else
           theta=0.0
        endif


       DO I=1,Nscalar
       DO J=1,Nscalar

        m1=PI*(Scalar(IJK,I)**3)*RO_S(IJK,I)/6.0
        m2=PI*(Scalar(IJK,J)**3)*RO_S(IJK,J)/6.0
        dav=(Scalar(IJk,I)+Scalar(IJk,J))/2.0
        c11=((theta*(m1+m2)**2)/(4*PI*m1*m2))**(0.5)*(4.0/dav)
        BETA_A(I,J)=aggregation_eff*PI*(dav**3)*G_0(IJK,I,J)*c11
        A(I)=A(I)+PI*(dav**3)*G_0(IJK,I,J)*c11*omega(J)*breakage_eff

      ENDDO
      ENDDO


! calculate S_bar
!
       DO L=0,2*Nscalar-1
        DO I=1,Nscalar
         DO J=1,Nscalar
           If(y(I)<=1.0E-6) THEN
              IF(y(J)<=1.0E-6) THEN
                S_bar(L)=S_bar(L)+0.0
              ELSE
                S_bar(L)=S_bar(L)+&
                BETA_a(I,J)*(1.0/2.0)*omega(I)*omega(J)*&
                ((y(Nscalar+J)**3)**(L/3.0))
              ENDIF
           ELSE
              IF(y(J)<=1.0E-6)THEN
                S_bar(L)=S_bar(L)+&
                BETA_a(I,J)*(1.0/2.0)*omega(I)*omega(J)*&
                ((y(Nscalar+I)**3)**(L/3.0))-&
                BETA_a(I,J)*(y(Nscalar+I)**L)*omega(I)*omega(J)
              ELSE
                S_bar(L)=S_bar(L)+&
                BETA_a(I,J)*(1.0/2.0)*omega(I)*omega(J)*&
                (((y(Nscalar+I)**3)+(y(Nscalar+J)**3))**(L/3.0))-&
                BETA_a(I,J)*(y(Nscalar+I)**L)*omega(I)*omega(J)
              ENDIF
            ENDIF
          ENDDO
        ENDDO

          DO I=1,Nscalar
              IF(y(I)<1.0E-6) THEN
                S_bar(L)=S_bar(L)+0.0
              ELSE
                S_bar(L)=S_bar(L)+&
                A(I)*(2**((3.0-L)/3.0))*(y(Nscalar+I)**L)*omega(I)-&
                A(I)*(y(Nscalar+I)**L)*omega(I)
              ENDIF
          ENDDO
       ENDDO

! calcalute inv(A)


       IF(NSCALAR==2) THEN
         IF(abs(y(3)-y(4))<=1.0E-4) THEN
          y(3)=y(3)*0.999
          y(4)=y(4)*1.001
         ENDIF
       ENDIF

       IF(NSCALAR==3) THEN
         IF(abs(y(4)-y(5))<=1.0E-4) THEN
          y(4)=y(4)*0.999
          y(5)=y(4)*1.001
         ENDIF
         IF(abs(y(5)-y(6))<=1.0E-4) THEN
          y(5)=y(5)*0.999
          y(6)=y(5)*1.001
         ENDIF
         IF(abs(y(4)-y(6))<=1.0E-4) THEN
          y(4)=y(4)*0.999
          y(6)=y(4)*1.001
         ENDIF
       ENDIF

       DO j=1,Nscalar
            matrix_a(1,j)=1.0
            matrix_a(2,j)=0.0
       ENDDO

       DO j=Nscalar+1,n
            matrix_a(1,j)=0.0
            matrix_a(2,j)=1.0
       ENDDO

       DO i=3,n
         DO j=1,Nscalar
          matrix_a(i,j)=(1.0-i+1.0)*(y(j+Nscalar)**((i-1)*1.0))
         ENDDO

         DO j=Nscalar+1,n
          matrix_a(i,j)=(i-1.0)*(y(j)**((i-2)*1.0))
         ENDDO
       ENDDO




       DO i=1,n
         inv_a(i,i)=1.0
       ENDDO

       call gaussj(matrix_a,n,np,inv_a,n,n)

       DO i=1,n
         IF(i<=Nscalar) THEN
           IF(y(i)>1.0E-3) matrix_c(i,i)=-2*(y(i+Nscalar)**3)
         ELSE
           IF(y(i-Nscalar)>1.0E-3) matrix_c(i,i)= 4*(y(i)**3)
         ENDIF
       ENDDO

       DO i=1,Nscalar
         IF(y(i)>1.0E-3) THEN
           matrix_c(i,i+Nscalar)=3*(y(i+Nscalar)**2)
           matrix_c(i+Nscalar,i)=-3*(y(i+Nscalar)**4)
         ENDIF
       ENDDO


       DO i = 1, n
          DO j = 1, n
            DO k = 1, n
            matrix_b(i,j) = matrix_b(i,j) + &
                      matrix_c(i,k)*inv_a(k,j)
            ENDDO
           ENDDO
       ENDDO

       DO i=1,Nscalar
         DO j=1,n
          dydx(i)=dydx(i)+k_v*matrix_b(i,j)*s_bar(j-1)
         ENDDO
       ENDDO

       DO i=Nscalar+1,n
           IF(y(i-Nscalar)>1.0E-6) THEN
           dydx(i)=-dydx(i-Nscalar)*y(i)/y(i-Nscalar)
          ELSE
           dydx(i)=dydx(i)+0.0
          ENDIF
         DO j=1,n
           IF(y(i-Nscalar)>1.0E-6) THEN
          dydx(i)=dydx(i)+(k_v/y(i-Nscalar))*(matrix_b(i,j)*s_bar(j-1))
           ELSE
            dydx(i)=dydx(i)+0.0
           ENDIF
         ENDDO
       ENDDO


      return
      END




