!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: READ_PAR_INPUT                                           !
!                                                                      !
! Purpose: Read the particle input and broadcasts the particle data to !
! respective processors.                                               !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_PAR_INPUT

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      use cdist
      use compar
      use desmpi
      use error_manager
      use functions
      use funits
      use geometry, only: NO_K
      use mpi_init_des, only: des_scatter_particle
      use mpi_utility

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: k
! index of particle
      INTEGER :: lcurpar
! local unit
      INTEGER, PARAMETER :: lunit=10
! local filename
      character(255) lfilename
! IO Status:
      INTEGER :: IOS
! Flag to indicate if file exists.
      LOGICAL :: lEXISTS
! Read dimension: 2D vs 3D data
      integer :: RDMN
!-----------------------------------------------


      CALL INIT_ERR_MSG("READ_PAR_INPUT")


      IOS = 0
      RDMN = merge(2,3,NO_K)

! Setup the file name based on distributed or serial IO.
      IF(bDIST_IO) THEN
         lFILENAME = ''
         WRITE(lFILENAME,'("particle_input_",I4.4,".dat")') myPE
      ELSE
         lFILENAME= "particle_input.dat"
      ENDIF

! Check the the file exists and open it.
      IF(bDIST_IO .OR. myPE == PE_IO) THEN
         INQUIRE(FILE=lFILENAME, EXIST=lEXISTS)
         IF(.NOT.LEXISTS) THEN
            WRITE(ERR_MSG, 1100)
            CALL FLUSH_ERR_MSG
            IOS = 1
         ELSE
            OPEN(CONVERT='BIG_ENDIAN',UNIT=lUNIT, FILE=lFILENAME, FORM="FORMATTED")
         ENDIF
      ENDIF

! Collect the error message and quit.
      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) CALL MFIX_EXIT(myPE)

 1100 FORMAT('Error 1100: FATAL - DEM particle input file not found!')

! Read the file
!----------------------------------------------------------------->>>
! In distributed IO the first line of the file will be number of
! particles in that processor
      IF (bdist_io) then
         read(lunit,*) pip
         DO lcurpar = 1,pip
            call set_normal(lcurpar)
            read (lunit,*) (des_pos_new(lcurpar,k),k=1,RDMN),&
               des_radius(lcurpar), ro_sol(lcurpar),&
               (des_vel_new(lcurpar,k),k=1,RDMN)
         ENDDO

! Serial IO (not bDIST_IO)
      ELSE
!----------------------------------------------------------------->>>

! Read into temporary variable and scatter
         IF (myPE .eq. PE_IO) THEN

! Allocate and initialize temporary variables.
            ALLOCATE (dpar_pos(particles,3)); dpar_pos=0.0
            ALLOCATE (dpar_vel(particles,3)); dpar_vel=0.0
            ALLOCATE (dpar_rad(particles));   dpar_rad=0.0
            ALLOCATE (dpar_den(particles));   dpar_den = 0.0
! Loop through the input file.
            DO lcurpar = 1, particles
               read (lunit,*,IOSTAT=IOS)                               &
               (dpar_pos(lcurpar,k),k=1,RDMN),dpar_rad(lcurpar),       &
               dpar_den(lcurpar),(dpar_vel(lcurpar,k),k=1,RDMN)

! Report read errors.
               IF(IOS > 0) THEN
                  WRITE(ERR_MSG,1200)
                  CALL FLUSH_ERR_MSG
                  EXIT
 1200 FORMAT('Error 1200: Error reported when reading particle input ',&
         'file.',/'A common error is 2D input for 3D cases.')

! Report End-of-File errors.
               ELSEIF(IOS < 0) THEN
                  WRITE(ERR_MSG,1201) &
                     trim(iVal(lcurpar)), trim(iVal(Particles))
                  CALL FLUSH_ERR_MSG
                  EXIT
 1201 FORMAT('Error 1201: Error reported when reading particle input ',&
         'file.',/'End-of-File found for particle ',A,' and ',A,1X,    &
         'entries are expected.')

               ENDIF

            ENDDO

         ENDIF

         CALL GLOBAL_ALL_SUM(IOS)
         IF(IOS /= 0) CALL MFIX_EXIT(myPE)

         CALL DES_SCATTER_PARTICLE

         IF(myPE == PE_IO) &
            deallocate (dpar_pos,dpar_vel,dpar_rad,dpar_den)

      ENDIF   ! end if/else bdist_io
!-----------------------------------------------------------------<<<

      IF(bDIST_IO .OR. myPE == PE_IO) CLOSE(lUNIT)

      CALL FINL_ERR_MSG()

      RETURN

      CALL MFIX_EXIT(myPE)

      END SUBROUTINE READ_PAR_INPUT

