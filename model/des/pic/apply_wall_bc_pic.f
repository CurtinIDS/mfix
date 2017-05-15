!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: APPLY_WALL_BC_PIC                                       !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose: Detect collisions between PIC particles and STLs.          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE APPLY_WALL_BC_PIC

! Global Variables:
!---------------------------------------------------------------------//
! Particle postions: Current/Previous
      use discretelement, only: DES_POS_NEW
! Paricle velocities
      use discretelement, only: DES_VEL_NEW
! Particle radius
      use discretelement, only: DES_RADIUS
! Max number of particles on this process
      use discretelement, only: MAX_PIP
! Map from particle to DES grid cell.
      use discretelement, only: DG_PIJK, DTSOLID
! Flag indicating that the index cell contains no STLs
      use stl, only: FACETS_AT_DG
! STL Vertices
      use stl, only: VERTEX
! STL Facet normals
      use stl, only: NORM_FACE
! Minimum velocity to offset gravitational forces
      use pic_bc, only: minVEL, minVEL_MAG, OoMinVEL_MAG
! Solids time step size
      use discretelement, only: DTSOLID
! Gravitational force vector
      use discretelement, only: GRAV

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ZERO, SMALL_NUMBER, ONE, UNDEFINED
      use functions, only: IS_NONEXISTENT

! Module Procedures:
!---------------------------------------------------------------------//
      use stl_functions_des

      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop counters.
      INTEGER NF, NP, LC1
! STL normal vector
      DOUBLE PRECISION :: NORM_PLANE(3)
! line is parameterized as p = p_ref + t * dir_line, t is line_param
      DOUBLE PRECISION :: LINE_T

! distance parcel travelled out of domain along the facet normal
      DOUBLE PRECISION :: DIST(3), tPOS(3)

      DOUBLE PRECISION :: DT_RBND, RADSQ
!......................................................................!

      double precision :: vv(3)

! Minimum velocity needed to offset gravity.
      minVEL = -DTSOLID*GRAV(:)
      minVEL_MAG = dot_product(minVEL,minVEL)
      OoMinVEL_MAG = 1.0d0/minVEL_MAG

      DT_RBND = 1.01d0*DTSOLID

      DO NP = 1, MAX_PIP

! Skip non-existent particles
         IF(IS_NONEXISTENT(NP)) CYCLE

! If no neighboring facet in the surrounding 27 cells, then exit
         IF(FACETS_AT_DG(DG_PIJK(NP))%COUNT < 1) CYCLE

         RADSQ = DES_RADIUS(NP)*DES_RADIUS(NP)

! Loop through the STLs associated with the DES grid cell.
         LC1_LP: DO LC1=1, FACETS_AT_DG(DG_PIJK(NP))%COUNT

! Get the index of the neighbor facet
            NF = FACETS_AT_DG(DG_PIJK(NP))%ID(LC1)

            NORM_PLANE = NORM_FACE(:,NF)

!            IF(dot_product(NORM_PLANE, DES_VEL_NEW(NP,:)) > 0.d0) &
!               CYCLE LC1_LP

            vv = VERTEX(1,:,NF)
            CALL INTERSECTLNPLANE(DES_POS_NEW(NP,:), DES_VEL_NEW(NP,:),&
                vv, NORM_FACE(:,NF), LINE_T)

            IF(abs(LINE_T) <= DT_RBND) THEN

! Project the parcel onto the facet.
               DIST = LINE_T*DES_VEL_NEW(NP,:)
               tPOS = DES_POS_NEW(NP,:) + DIST

! Avoid collisions with STLs that are not next to the parcel
               IF(HIT_FACET(VERTEX(:,:,NF), tPOS)) THEN

! Correct the position of particles found too close to the STL.
                  IF(dot_product(DIST,DIST) <= RADSQ .OR. &
                     LINE_T <= ZERO) THEN
                     DES_POS_NEW(NP,:) = DES_POS_NEW(NP,:) + DIST(:) + &
                        DES_RADIUS(NP)*NORM_PLANE
                  ENDIF

! Reflect the parcel.
                  CALL PIC_REFLECT_PART(NP, NORM_PLANE(:))
               ENDIF
            ENDIF

         ENDDO LC1_LP
      END DO

      RETURN

      contains


!......................................................................!
      LOGICAL FUNCTION HIT_FACET(VERTS, POINT)

      DOUBLE PRECISION, INTENT(IN) :: VERTS(3,3)
      DOUBLE PRECISION, INTENT(IN) :: POINT(3)

      DOUBLE PRECISION :: V0(3), V1(3), V2(3)
      DOUBLE PRECISION :: d00, d01, d02, d11, d12

      DOUBLE PRECISION :: OoDEN, u, v

      V0 = VERTS(3,:) - VERTS(1,:)
      V1 = VERTS(2,:) - VERTS(1,:)
      V2 = POINT - VERTS(1,:)

      d00 = dot_product(V0,V0)
      d01 = dot_product(V0,V1)
      d02 = dot_product(V0,V2)

      d11 = dot_product(V1,V1)
      d12 = dot_product(V1,V2)

      OoDEN = 1.0d0/(d00*d11 - d01*d01)
      u = (d11*d02 - d01*d12)*OoDEN
      v = (d00*d12 - d01*d02)*OoDEN

      HIT_FACET = ((u>=0.0d0) .AND. (v>=0.0d0) .AND. (u+v < 1.0d0))

      RETURN
      END FUNCTION HIT_FACET

      END SUBROUTINE APPLY_WALL_BC_PIC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PIC_REFLECT_PARTICLE                                    C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PIC_REFLECT_PART(LL, WALL_NORM)

      USE discretelement, only : dimn, DES_VEL_NEW
      USE mfix_pic, only : MPPIC_COEFF_EN_WALL, MPPIC_COEFF_ET_WALL

! Minimum velocity to offset gravitational forces
      use pic_bc, only: minVEL, minVEL_MAG, OoMinVEL_MAG


      use param1, only: SMALL_NUMBER

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, INTENT(IN) :: LL
      DOUBLE PRECISION, INTENT(IN) ::  WALL_NORM(DIMN)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      !magnitude of pre-collisional normal and tangential velocity components
      DOUBLE PRECISION :: VEL_NORMMAG_APP

      !pre collisional normal and tangential velocity components in vector format
      !APP ==> approach
      DOUBLE PRECISION :: VEL_NORM_APP(DIMN), VEL_TANG_APP(DIMN)


      !post collisional normal and tangential velocity components in vector format
      !SEP ==> separation
      DOUBLE PRECISION :: VEL_NORM_SEP(DIMN), VEL_TANG_SEP(DIMN)

      DOUBLE PRECISION :: COEFF_REST_EN, COEFF_REST_ET

! Minimum particle velocity with respect to gravity. This ensures that
! walls "push back" to offset gravitational forces.
      DOUBLE PRECISION :: projGRAV(3)

      INTEGER :: LC

!-----------------------------------------------





      VEL_NORMMAG_APP = DOT_PRODUCT(WALL_NORM(:), DES_VEL_NEW(LL,:))

! Currently assuming that wall is at rest. Needs improvement for moving wall

      VEL_NORM_APP(:) = VEL_NORMMAG_APP*WALL_NORM(:)
      VEL_TANG_APP(:) = DES_VEL_NEW(LL,:) - VEL_NORM_APP(:)


!post collisional velocities

      COEFF_REST_EN = MPPIC_COEFF_EN_WALL
      COEFF_REST_ET = MPPIC_COEFF_ET_WALL

      VEL_NORM_SEP(:) = -COEFF_REST_EN*VEL_NORM_APP(:)
      VEL_TANG_SEP(:) =  COEFF_REST_ET*VEL_TANG_APP(:)

      DES_VEL_NEW(LL,:) = VEL_NORM_SEP(:) + VEL_TANG_SEP(:)

      IF(dot_product(WALL_NORM, minVEL) > 0.0d0)THEN

! Projection of rebound velocity onto mininum velocity.
         projGRAV = dot_product(minVEL, DES_VEL_NEW(LL,:))*OoMinVEL_MAG
         projGRAV = minVEL*projGRAV

         IF(dot_product(projGRAV, projGRAV) < minVEL_MAG) THEN
            DO LC=1,3
               IF(minVEL(LC) > SMALL_NUMBER) THEN
                  DES_VEL_NEW(LL,LC)=max(DES_VEL_NEW(LL,LC),minVEL(LC))
               ELSEIF(minVEL(LC) < -SMALL_NUMBER) THEN
                  DES_VEL_NEW(LL,LC)=min(DES_VEL_NEW(LL,LC),minVEL(LC))
               ENDIF
            ENDDO

         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE PIC_REFLECT_PART


!
      SUBROUTINE OUT_THIS_STL(LC)

      Use usr

! STL Vertices
      use stl, only: VERTEX
! STL Facet normals
      use stl, only: NORM_FACE


      IMPLICIT NONE
!-----------------------------------------------
      integer, INTENT(IN) :: lc

      logical :: lExists
      character(len=8) :: IDX


      write(idx,"(I8.8)") LC
      inquire(file='geo_'//idx//'.stl', EXIST=lExists)

      IF(lExists) RETURN

      open(unit=555, file='geo_'//idx//'.stl', status='UNKNOWN')
      write(555,*) 'solid vcg'
      write(555,*) '   facet normal ', NORM_FACE(:,LC)
      write(555,*) '      outer loop'
      write(555,*) '         vertex ', VERTEX(1,1:3,LC)
      write(555,*) '         vertex ', VERTEX(2,1:3,LC)
      write(555,*) '         vertex ', VERTEX(3,1:3,LC)
      write(555,*) '      endloop'
      write(555,*) '   endfacet'
      close(555)

      RETURN
      END SUBROUTINE OUT_THIS_STL




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: CHECK_IF_PARCEL_OVERLAPS_STL                            C
!  Authors: Rahul Garg                               Date: 21-Mar-2014 C
!                                                                      C
!  Purpose: This subroutine is special written to check if a particle  C
!          overlaps any of the STL faces. The routine exits on         C
!          detecting an overlap. It is called after initial            C
!          generation of lattice configuration to remove out of domain C
!          particles                                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_IF_PARCEL_OVERLAPS_STL(POS, OVERLAP_EXISTS)

      USE discretelement, only: MAX_RADIUS
      USE geometry, only: do_k

      USE stl_functions_des
      USE desgrid
      USE stl

      Implicit none

      DOUBLE PRECISION, INTENT(IN) :: POS(3)
      LOGICAL, INTENT(OUT) :: OVERLAP_EXISTS

      INTEGER I, J, K, IJK, NF, LC

! line is parameterized as p = p_ref + t * dir_line, t is line_param
      DOUBLE PRECISION :: LINE_T
      DOUBLE PRECISION :: DIST(3), RADSQ
      double precision :: vv(3)

      K = 1
      IF(DO_K) K = min(DG_KEND2,max(DG_KSTART2,KOFPOS(POS(3))))
      J = min(DG_JEND2,max(DG_JSTART2,JOFPOS(POS(2))))
      I = min(DG_IEND2,max(DG_ISTART2,IOFPOS(POS(1))))

      IJK = DG_FUNIJK(I,J,K)
      IF(FACETS_AT_DG(IJK)%COUNT < 1) RETURN

      OVERLAP_EXISTS = .TRUE.

! Pad the max radius to keep parcels sufficiently away from walls.
      RADSQ = (1.1d0*MAX_RADIUS)**2

! Parametrize a line as p = p_0 + t normal and intersect with the
! triangular plane.

! If t>0, then point is on the non-fluid side of the plane, if the
! plane normal is assumed to point toward the fluid side.
      DO LC = 1, FACETS_AT_DG(IJK)%COUNT
         NF = FACETS_AT_DG(IJK)%ID(LC)

         vv = VERTEX(1,:,NF)
         CALL INTERSECTLNPLANE(POS, NORM_FACE(:,NF), &
            vv, NORM_FACE(:,NF), LINE_T)

! Orthogonal projection puts the point outside of the domain or less than
! one particle radius to the facet.
         DIST = LINE_T*NORM_FACE(:,NF)
         IF(LINE_T > ZERO .OR. dot_product(DIST,DIST) <= RADSQ) RETURN

      ENDDO

      OVERLAP_EXISTS = .FALSE.

      RETURN
      END SUBROUTINE CHECK_IF_PARCEL_OVERLAPS_STL


!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE write_this_facet_and_parcel(FID, position, velocity)
      USE run
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE constant
      USE cutcell
      USE funits
      USE indices
      USE physprop
      USE parallel
      USE stl
      USE stl_functions_des
      Implicit none
      !facet id and particle id
      double precision, intent(in), dimension(dimn) :: position, velocity
      Integer, intent(in) :: fid
      Integer :: stl_unit, vtp_unit , k
      CHARACTER(LEN=100) :: stl_fname, vtp_fname
      real :: temp_array(3)

      stl_unit = 1001
      vtp_unit = 1002

      WRITE(vtp_fname,'(A,"_OFFENDING_PARTICLE",".vtp")') TRIM(RUN_NAME)
      WRITE(stl_fname,'(A,"_STL_FACE",".stl")') TRIM(RUN_NAME)

      open(vtp_unit, file = vtp_fname, form='formatted',convert='big_endian')
      open(stl_unit, file = stl_fname, form='formatted',convert='big_endian')

      write(vtp_unit,"(a)") '<?xml version="1.0"?>'
      write(vtp_unit,"(a,es24.16,a)") '<!-- time =',s_time,'s -->'
      write(vtp_unit,"(a,a)") '<VTKFile type="PolyData"',&
           ' version="0.1" byte_order="LittleEndian">'
      write(vtp_unit,"(3x,a)") '<PolyData>'
      write(vtp_unit,"(6x,a,i10.10,a,a)")&
           '<Piece NumberOfPoints="',1,'" NumberOfVerts="0" ',&
           'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
      write(vtp_unit,"(9x,a)")&
           '<PointData Scalars="Diameter" Vectors="Velocity">'
      write(vtp_unit,"(12x,a)")&
           '<DataArray type="Float32" Name="Diameter" format="ascii">'
      write (vtp_unit,"(15x,es13.6)") (1.d0)
      write(vtp_unit,"(12x,a)") '</DataArray>'

      temp_array = zero
      temp_array(1:DIMN) = velocity(1:dimn)
      write(vtp_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
           'Name="Velocity" NumberOfComponents="3" format="ascii">'
      write (vtp_unit,"(15x,3(es13.6,3x))")&
           ((temp_array(k)),k=1,3)
      write(vtp_unit,"(12x,a,/9x,a)") '</DataArray>','</PointData>'
      ! skip cell data
      write(vtp_unit,"(9x,a)") '<CellData></CellData>'

      temp_array = zero
      temp_array(1:dimn) = position(1:dimn)
      write(vtp_unit,"(9x,a)") '<Points>'
      write(vtp_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
           'Name="Position" NumberOfComponents="3" format="ascii">'
      write (vtp_unit,"(15x,3(es13.6,3x))")&
           ((temp_array(k)),k=1,3)
      write(vtp_unit,"(12x,a,/9x,a)")'</DataArray>','</Points>'
      ! Write tags for data not included (vtp format style)
      write(vtp_unit,"(9x,a,/9x,a,/9x,a,/9x,a)")'<Verts></Verts>',&
           '<Lines></Lines>','<Strips></Strips>','<Polys></Polys>'
      write(vtp_unit,"(6x,a,/3x,a,/a)")&
           '</Piece>','</PolyData>','</VTKFile>'

      !Now write the facet info

      write(stl_unit,*)'solid vcg'

      write(stl_unit,*) '   facet normal ', NORM_FACE(1:3,FID)
      write(stl_unit,*) '      outer loop'
      write(stl_unit,*) '         vertex ', VERTEX(1,1:3,FID)
      write(stl_unit,*) '         vertex ', VERTEX(2,1:3,FID)
      write(stl_unit,*) '         vertex ', VERTEX(3,1:3,FID)
      write(stl_unit,*) '      endloop'
      write(stl_unit,*) '   endfacet'

      write(stl_unit,*)'endsolid vcg'

      close(vtp_unit, status = 'keep')
      close(stl_unit, status = 'keep')
      write(*,*) 'wrote a facet and a parcel. now waiting'
      read(*,*)
      end SUBROUTINE write_this_facet_and_parcel

