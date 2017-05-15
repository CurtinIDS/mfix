!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_PS_GRAD_PIC                                        !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose:                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INTERPOLATE_PIC

      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_GARG

      IMPLICIT NONE

      SELECT CASE(DES_INTERP_SCHEME_ENUM)
      CASE(DES_INTERP_GARG) ; CALL INTERPOLATE_PIC0
      CASE DEFAULT; CALL INTERPOLATE_PIC1
      END SELECT

      RETURN
      END SUBROUTINE INTERPOLATE_PIC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: INTERPOLATE_PIC_NONE                                    !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose:                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INTERPOLATE_PIC1

      use discretelement, only: MAX_PIP, PIJK
      USE geometry, only: DO_K
      use fldvar, only: EP_G, u_s, v_s, w_s
      use functions, only: FLUID_AT
      use functions, only: IS_NONEXISTENT
      use mfix_pic, only: AVGVEL => AVGSOLVEL_P
      use mfix_pic, only: PS_GRAD, PS_FORCE_PIC
      use mfix_pic, only: EPG_P
! Flag to use interpolation
      use particle_filter, only: DES_INTERP_ON
! Interpolation cells and weights
      use particle_filter, only: FILTER_CELL, FILTER_WEIGHT
      use physprop, only: mmax

      implicit none

! Local Variables
!---------------------------------------------------------------------//
! Interpolation weight
      DOUBLE PRECISION :: WEIGHT
! temporary variables used to calculate pressure at scalar cell edge
      INTEGER :: IJK, NP, LC
! Loop bound for filter
      INTEGER :: LP_BND
!......................................................................!

! Loop bounds for interpolation.
      LP_BND = merge(27,9,DO_K)


      DO NP=1,MAX_PIP

! Skip parcels that don't exist or are in a non-fluid cell.
         IF(IS_NONEXISTENT(NP) .OR. .NOT.FLUID_AT(PIJK(NP,4))) THEN
            AVGVEL(:,NP) = 0.0d0
            PS_GRAD(:,NP) = 0.0d0
            EPG_P(NP) = 1.0d0

! Interpolated solids pressure and average solids velocity.
         ELSEIF(DES_INTERP_ON) THEN
            AVGVEL(:,NP) = 0.0d0
            PS_GRAD(:,NP) = 0.0d0
            EPG_P(NP) = 0.0d0
            DO LC=1,LP_BND
               IJK = FILTER_CELL(LC,NP)
               WEIGHT = FILTER_WEIGHT(LC,NP)
! Gas phase velocity.
               AVGVEL(1,NP) = AVGVEL(1,NP) + U_s(IJK,MMAX+1)*WEIGHT
               AVGVEL(2,NP) = AVGVEL(2,NP) + V_s(IJK,MMAX+1)*WEIGHT
               AVGVEL(3,NP) = AVGVEL(3,NP) + W_s(IJK,MMAX+1)*WEIGHT
! Particle normal stress.
               PS_GRAD(:,NP) = PS_GRAD(:,NP)+PS_FORCE_PIC(:,IJK)*WEIGHT
! TO BE REMOVED
               EPG_P(NP) = EPG_P(NP) + EP_G(IJK)*WEIGHT
            ENDDO

! Centroid method. 
        ELSE
            IJK = PIJK(NP,4)
! Gas phase velocity.
            AVGVEL(1,NP) = U_s(IJK,MMAX+1)
            AVGVEL(2,NP) = V_s(IJK,MMAX+1)
            AVGVEL(3,NP) = W_s(IJK,MMAX+1)
! Particle normal stress.
            PS_GRAD(:,NP) = PS_FORCE_PIC(:,IJK)
! TO BE REMOVED
            EPG_P(NP) = EP_G(IJK)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE INTERPOLATE_PIC1


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_PS_GRAD_PIC                                        !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose:                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INTERPOLATE_PIC0

      USE bc
      USE compar
      USE constant
      USE cutcell
      USE derived_types, only: pic
      USE discretelement
      USE fldvar, only: ep_g, u_g, v_g, w_g
      USE fun_avg
      USE functions
      USE geometry
      USE indices
      USE interpolation
      USE mfix_pic
      USE parallel
      USE param
      USE param1
      USE physprop
      USE run
      USE sendrecv

      implicit none

      ! general i, j, k indices
      INTEGER I, J, K, IJK, IPJK, IJPK, IJKP, IDIM, M

      ! temporary variables used to calculate pressure at scalar cell edge
      DOUBLE PRECISION avg_factor, VOL_TOT_VEC(3), VOL_TOT_SCAL

      integer :: korder, iw,ie,js,jn,kb,ktp, onew, pcell(3), cur_ijk, NP, nindx

      integer :: ii,jj,kk, ipjpk, ijpkp, ipjkp, ipjpkp

      double precision :: vol_ijk, vol_ipjk, vol_ijpk, vol_ipjpk
      double precision :: vol_ijkp, vol_ipjkp, vol_ijpkp, vol_ipjpkp


      CALL SET_INTERPOLATION_SCHEME(2)

      KORDER = merge ( 1, 2, no_k) !1+(DIMN-2)

! avg_factor=0.25 (in 3D) or =0.5 (in 2D)
      !avg_factor = 0.25d0*(dimn-2) + 0.5d0*(3-dimn)
      AVG_FACTOR = merge(0.5d0, 0.25D0, NO_K)

      do ijk = ijkstart3,ijkend3

         if(.not.fluid_at(ijk) .or. pinc(ijk).eq.0) cycle
         i = i_of(ijk)
         j = j_of(ijk)
         k = k_of(ijk)

         pcell(1) = i-1
         pcell(2) = j-1
         pcell(3) = merge(1, k-1, no_k) ! =k-1 (in 3d) or =1 (in 2d)

         CALL SET_INTERPOLATION_STENCIL(PCELL,IW,IE,JS,JN,KB,          &
            KTP,INTERP_SCHEME,MERGE(2,3,NO_K),ORDERNEW = ONEW)

! Compute velocity at grid nodes and set the geometric stencil
         DO K = 1, merge(1, ONEW, NO_K)!  (3-DIMN)*1+(DIMN-2)*ONEW
         DO J = 1,ONEW
         DO I = 1,ONEW

            II = IW + I-1
            JJ = JS + J-1
            KK = KB + K-1

            CUR_IJK = FUNIJK_MAP_C(II,JJ,KK)

            IPJK    = FUNIJK_MAP_C(II+1,JJ,KK)
            IJPK    = FUNIJK_MAP_C(II,JJ+1,KK)
            IPJPK   = FUNIJK_MAP_C(II+1,JJ+1,KK)

            VOL_IJK = ZERO
            VOL_IPJK = ZERO
            VOL_IJPK = ZERO
            VOL_IPJPK = ZERO

            VOL_IJKP = ZERO
            VOL_IPJKP = ZERO
            VOL_IJPKP = ZERO
            VOL_IPJPKP = ZERO

            IF(FLUID_AT(CUR_IJK)) VOL_IJK   = VOL(CUR_IJK)
            IF(FLUID_AT(IPJK))    VOL_IPJK  = VOL(IPJK)
            IF(FLUID_AT(IJPK))    VOL_IJPK  = VOL(IJPK)
            IF(FLUID_AT(IPJPK))   VOL_IPJPK = VOL(IPJPK)

            IF(DO_K) THEN
               IJKP   = FUNIJK_MAP_C(II,JJ,KK+1)
               IJPKP  = FUNIJK_MAP_C(II,JJ+1,KK+1)
               IPJKP  = FUNIJK_MAP_C(II+1,JJ,KK+1)
               IPJPKP = FUNIJK_MAP_C(II+1,JJ+1,KK+1)

               IF(FLUID_AT(IJKP))   VOL_IJKP   = VOL(IJKP)
               IF(FLUID_AT(IPJKP))  VOL_IPJKP  = VOL(IPJKP)
               IF(FLUID_AT(IJPKP))  VOL_IJPKP  = VOL(IJPKP)
               IF(FLUID_AT(IPJPKP)) VOL_IPJPKP = VOL(IPJPKP)
            ENDIF

            GSTENCIL(I,J,K,1) = XE(II)
            GSTENCIL(I,J,K,2) = YN(JJ)
            GSTENCIL(I,J,K,3) = merge(DZ(1), ZT(KK), NO_K)

            VOL_TOT_SCAL = ZERO

            VOL_TOT_SCAL = VOL_IJK + VOL_IPJK + VOL_IJPK + VOL_IPJPK + &
               VOL_IJKP + VOL_IPJKP + VOL_IJPKP + VOL_IPJPKP

            VOL_TOT_VEC = ZERO

            VOL_TOT_VEC(1) = VOL(CUR_IJK) + VOL(IJPK)
            VOL_TOT_VEC(2) = VOL(CUR_IJK) + VOL(IPJK)

            DO M = MMAX+1,DES_MMAX+MMAX
               VEL_SOL_STENCIL(I,J,K,1,M) = &
                  PIC_U_S(CUR_IJK,M)*VOL(CUR_IJK) + &
                  PIC_U_S(IJPK,M)*VOL(IJPK)

               VEL_SOL_STENCIL(I,J,K,2,M) = &
                  PIC_V_S(CUR_IJK,M)*VOL(CUR_IJK) + &
                  PIC_V_S(IPJK,M)*VOL(IPJK)
            ENDDO

            SSTENCIL(I,J,K) = EP_G(CUR_IJK)*VOL_IJK + &
               EP_G(IPJK)*VOL_IPJK + EP_G(IJPK)*VOL_IJPK + &
               EP_G(IPJPK)*VOL_IPJPK

            PSGRADSTENCIL(I,J,K,1) = &
               PS_FORCE_PIC(1,CUR_IJK)*VOL(CUR_IJK) + &
               PS_FORCE_PIC(1,IJPK)*VOL(IJPK)

            PSGRADSTENCIL(I,J,K,2) = &
               PS_FORCE_PIC(2,CUR_IJK)*VOL(CUR_IJK) + &
               PS_FORCE_PIC(2,IPJK)*VOL(IPJK)

            VSTENCIL(I,J,K,1) = &
               U_G(CUR_IJK)*VOL(CUR_IJK) + U_G(IJPK)*VOL(IJPK)

            VSTENCIL(I,J,K,2) = &
               V_G(CUR_IJK)*VOL(CUR_IJK) + V_G(IPJK)*VOL(IPJK)

            IF(DO_K) THEN
               VOL_TOT_VEC(1) = VOL_TOT_VEC(1) + VOL(IJKP) + VOL(IJPKP)
               VOL_TOT_VEC(2) = VOL_TOT_VEC(2) + VOL(IJKP) + VOL(IPJKP)
               VOL_TOT_VEC(3) = &
                  VOL(CUR_IJK) + VOL(IPJK) + VOL(IJPK) + VOL(IPJPK)

               SSTENCIL(I,J,K) = SSTENCIL(I,J,K) + &
                  EP_G(IJKP)*VOL_IJKP + EP_G(IPJKP)*VOL_IPJKP + &
                  EP_G(IJPKP)*VOL_IJPKP + EP_G(IPJPKP)*VOL_IPJPKP

               PSGRADSTENCIL(I,J,K,1) = PSGRADSTENCIL(I,J,K,1) + &
                  PS_FORCE_PIC(1,IJKP)*VOL(IJKP) + &
                  PS_FORCE_PIC(1,IJPKP)*VOL(IJPKP)

               PSGRADSTENCIL(I,J,K,2) = PSGRADSTENCIL(I,J,K,2) + &
                  PS_FORCE_PIC(2,IJKP)*VOL(IJKP) + &
                  PS_FORCE_PIC(2,IPJKP)*VOL(IPJKP)

               PSGRADSTENCIL(I,J,K,3) = &
                  PS_FORCE_PIC(3,CUR_IJK)*VOL(CUR_IJK)+ &
                  PS_FORCE_PIC(3,IJPK)*VOL(IJPK)+       &
                  PS_FORCE_PIC(3,IPJK)*VOL(IPJK)+       &
                  PS_FORCE_PIC(3,IPJPK)*VOL(IPJPK)

               VSTENCIL(I,J,K,1) = VSTENCIL(I,J,K,1) + &
                  U_G(IJKP)*VOL(IJKP) + U_G(IJPKP)*VOL(IJPKP)

               VSTENCIL(I,J,K,2) = VSTENCIL(I,J,K,2) + &
                  V_G(IJKP)*VOL(IJKP) + V_G(IPJKP)*VOL(IPJKP)

               VSTENCIL(I,J,K,3) = W_G(CUR_IJK)*VOL(CUR_IJK)+&
                  W_G(IJPK)*VOL(IJPK) + W_G(IPJK)*VOL(IPJK) + &
                  W_G(IPJPK)*VOL(IPJPK)

               DO M = MMAX+1, DES_MMAX+MMAX
                  VEL_SOL_STENCIL(I,J,K,1,M) = &
                     VEL_SOL_STENCIL(I,J,K,1,M) + &
                     PIC_U_S(IJKP,M)*VOL(IJKP) + & 
                     PIC_U_S(IJPKP,M)*VOL(IJPKP)

                  VEL_SOL_STENCIL(I,J,K,2, M) = &
                     VEL_SOL_STENCIL(I,J,K,2,M) + &
                     PIC_V_S(IJKP,M)*VOL(IJKP) +&
                     PIC_V_S(IPJKP,M)*VOL(IPJKP)

                  VEL_SOL_STENCIL(I,J,K,3, M) = &
                     PIC_W_S(CUR_IJK,M)*VOL(CUR_IJK) + &
                     PIC_W_S(IJPK,M)*VOL(IJPK) + &
                     PIC_W_S(IPJK,M)*VOL(IPJK) + &
                     PIC_W_S(IPJPK,M)*VOL(IPJPK)
               ENDDO
            ELSE
               PSGRADSTENCIL(I,J,K,3) = 0.D0
               VEL_SOL_STENCIL(I,J,K,3, MMAX+1:DES_MMAX+MMAX) = 0.D0
               VSTENCIL(I,J,K,3) = 0.D0
            ENDIF




            DO IDIM = 1, merge(2,3,NO_K)
               IF(VOL_TOT_VEC(IDIM).GT.ZERO)  THEN
                  psgradstencil(i,j,k,idim) = &
                     psgradstencil(i,j,k,idim)/VOL_TOT_VEC(idim)

                  VEL_SOL_STENCIL(i,j,k,idim, MMAX+1:DES_MMAX+MMAX) = &
                     VEL_SOL_STENCIL(i,j,k,idim, MMAX+1:DES_MMAX+MMAX)/&
                     VOL_TOT_VEC(idim)

                  vstencil(i,j,k,idim) = &
                     vstencil(i,j,k,idim)/VOL_TOT_VEC(idim)
               ENDIF
            ENDDO


            IF(VOL_TOT_SCAL.GT.ZERO) &
               SSTENCIL(I,J,K) = SSTENCIL(I,J,K)/VOL_TOT_SCAL

         enddo
         enddo
         enddo



! Loop through particles in the cell
         do nindx = 1,pinc(ijk)
            np = pic(ijk)%p(nindx)
            m = pijk(np,5)

! Interpolate particle stress to particle position
            if(NO_K) then !2-D
               call interpolator(gstencil(1:onew,1:onew,1,1:dimn), &
               psgradstencil(1:onew,1:onew,1,1:dimn), &
               des_pos_new(np,1:dimn),PS_GRAD(1:dimn,np),  &
               onew,interp_scheme,weightp)
            else !3-D, diff in psgradstencil size
               call interpolator(gstencil(1:onew,1:onew,1:onew,1:dimn), &
               psgradstencil(1:onew,1:onew,1:onew,1:dimn), &
               des_pos_new(np,1:dimn),PS_GRAD(1:dimn,np),  &
               onew,interp_scheme,weightp)
            endif

            do idim = 1,  merge(2,3,NO_K)
               AVGSOLVEL_P(IDIM,NP) = ARRAY_DOT_PRODUCT(               &
                  VEL_SOL_STENCIL(:,:,:,IDIM,M),WEIGHTP(:,:,:))
            ENDDO

            EPG_P(NP) = ARRAY_DOT_PRODUCT(SSTENCIL(:,:,:),WEIGHTP(:,:,:))

         ENDDO
      ENDDO


      END SUBROUTINE INTERPOLATE_PIC0
