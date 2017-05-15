!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_CLUSTER                                            !
!                                                                      !
!  Purpose: Common elements for MFIX-DEM cluster identification        !
!  condition.                                                          !
!                                                                      !
!  Author: J.Galvin, J.Musser                         Date:  Nov-12    !
!  Modified: S. Benyahia                              Date:  Dec-12    !
!                                                                      !
!  Comments: Info on clusters such as average void fraction, Re and    !
!  Comments: cluster size can now be printed from this file            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE DES_CLUSTER


!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE compar
      USE fldvar
      USE physprop
      USE discretelement

      use mpi_utility
      use desmpi_wrapper
      use desmpi
      use compar
      use desmpi
      use parallel
      use sendrecv


      IMPLICIT NONE
!-----------------------------------------------


! Number of clusters.
      INTEGER ClusterCount

! Number of base particles held in current search list
      INTEGER pSearchHistoryCount

!----------------------------------------------->>>
      TYPE PARTICLE_TYPE
! Global index of particle.
         INTEGER :: ID
! Successor in particle list
         TYPE (PARTICLE_TYPE), POINTER :: next_particle
      END TYPE PARTICLE_TYPE
!-----------------------------------------------<<<

! Basic particle link-list.
      TYPE(PARTICLE_TYPE), POINTER :: PSEARCH_HISTORY_LL

!----------------------------------------------->>>
      TYPE CLUSTER_TYPE
! Cluster index.
         INTEGER :: ID
! Number of particles in cluster.
         INTEGER :: ParticleCount
! Successor in cluster link-list.
         TYPE (CLUSTER_TYPE), POINTER :: next_cluster
! Particle link-list
         TYPE (PARTICLE_TYPE), POINTER :: PARTICLE_LL
      END TYPE CLUSTER_TYPE
!-----------------------------------------------<<<

! Cluster link-list.
      TYPE(CLUSTER_TYPE), POINTER :: CLUSTER_LL

!<<<-------------------  DMP Related Variables  -------------------->>>!

! Process that does global cluster calculations.
      integer, parameter :: clusterPE = 0

      TYPE pType
         INTEGER :: map
         INTEGER :: id
! Successor in particle list
         TYPE (pType), POINTER :: next
      END TYPE pType


      TYPE cType
! Number of particles in cluster.
         INTEGER :: size
! Particle link-list
         TYPE (pType), POINTER :: particle
      END TYPE cType

! Cluster array.
      TYPE(cType), dimension(:), allocatable, target :: clusters
      INTEGER clusterCount_all

! Send/Recv information:
!----------------------->>

! Number of data entrires being sent from the local process to Root
! during a call to cluster_gather.
      INTEGER send_cnt

! Number of data entries being received by root from each process
! during a call to cluster_gather.
      INTEGER, dimension(:), allocatable :: recv_cnt

! Displacement of data from all processes during a cluster_gather.
      INTEGER, dimension(:), allocatable :: recv_dsp

! Size of array needed by Root to accept all data during a gather call.
      INTEGER recv_sum


! Module interfaces:
!----------------------->>

      interface getClusterParticleData
         module procedure getClusterParticleData_1i
         module procedure getClusterParticleData_2i
         module procedure getClusterParticleData_1d
         module procedure getClusterParticleData_2d
      end interface

      interface getClusterFieldData
         module procedure getClusterFieldData_1d
         module procedure getClusterFieldData_3d
      end interface

      interface sendClusterData
         module procedure sendClusterData_1d
      end interface


      contains

!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE CREATE_CLUSTER(cluster)

      TYPE(CLUSTER_TYPE), INTENT(OUT), POINTER :: cluster

      ALLOCATE(cluster)

      NULLIFY(cluster%next_cluster)
      NULLIFY(cluster%PARTICLE_LL)

      if(ClusterCount == 0) then
! no clusters have been created/identified
         if(associated(CLUSTER_LL)) then
            print*, ' Error - cluster pointer already associated!'
            CALL MFIX_EXIT(myPE)
         else
            ClusterCount = 1
            Cluster%ParticleCount = 0
            cluster%ID = ClusterCount
! create first cluster of linked list of clusters. with cluster_ll
! always being the 'first' in the list
            CLUSTER_LL => cluster
         endif
      else
         if(.NOT.associated(CLUSTER_LL)) then
            print*, ' Error - cluster pointer is not associated!'
            CALL MFIX_EXIT(myPE)
         else
            ClusterCount = ClusterCount + 1
            Cluster%ParticleCount = 0
            cluster%ID = ClusterCount
! establish the link between the new cluster and existing cluster list
            cluster%next_cluster => CLUSTER_LL
! reassign/point cluster_ll to be 'first' in the linked list
            CLUSTER_LL => cluster
! as a result the linked list of clusters is created from 'bottom-up'
! with the new cluster always being inserted before any existing
! clusters
         endif
      ENDIF

      END SUBROUTINE CREATE_CLUSTER


!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE GetTopCluster(cluster)

      TYPE(CLUSTER_TYPE), INTENT(INOUT), POINTER :: cluster

      NULLIFY(cluster)

! first check that a list has been formed
      if(.NOT.associated(CLUSTER_LL) .AND. ClusterCount == 0) then
         write(*,"(//,3x,A,//)") 'No clusters to delete!'
      elseif(.NOT.associated(CLUSTER_LL) .AND. ClusterCount /= 0) then
         write(*,"(//,3x,A,//)") ' Error with ClusterCount and pointer - Er:1'
      elseif(associated(CLUSTER_LL) .AND. ClusterCount == 0) then
         write(*,"(//,3x,A,//)") ' Error with ClusterCount and pointer - Er:2'
      else
! now identify top cluster
         CALL getNextCluster(cluster)
      endif
      END SUBROUTINE GetTopCluster


!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE DeleteTopCluster(cluster)

      TYPE(CLUSTER_TYPE), INTENT(INOUT), POINTER :: cluster

! dbg
!      write(*,"(/A,I7)")' Deleting top cluster: ',cluster%ID
      CALL DELETE_PARTICLES_IN_CLUSTER(cluster)

      if(associated(cluster%next_cluster)) then
         CLUSTER_LL => cluster%next_cluster
      else
         NULLIFY(CLUSTER_LL)
      endif

      ClusterCount = ClusterCount - 1
      NULLIFY(cluster%next_cluster)
      NULLIFY(cluster%PARTICLE_LL)
      DEALLOCATE(cluster)
      NULLIFY(cluster)

      if(ClusterCount <0) then
         write(*,"()")' ClusterCount < 0'
         CALL MFIX_EXIT(myPE)
      endif

      END SUBROUTINE DeleteTopCluster


!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE DELETE_CLUSTERS()

      TYPE(CLUSTER_TYPE), POINTER :: cluster
      INTEGER cL, cDeleted

      NULLIFY(cluster)
      cDeleted = 0

      if(.NOT.associated(CLUSTER_LL) .AND. ClusterCount == 0) then
!         write(*,"(//,3x,A,//)") 'No clusters to delete!'
      elseif(.NOT.associated(CLUSTER_LL) .AND. ClusterCount /= 0) then
         write(*,"(//,3x,A,//)") ' Error with ClusterCount and pointer - Er:1'
      elseif(associated(CLUSTER_LL) .AND. ClusterCount == 0) then
         write(*,"(//,3x,A,//)") ' Error with ClusterCount and pointer - Er:2'
      else
         do cL =1, ClusterCount
            CALL getNextCluster(cluster)
! dbg
!            write(*,"(/A,I7)")' Deleting cluster: ',cluster%ID
            CALL DELETE_PARTICLES_IN_CLUSTER(cluster)

            if(associated(cluster%next_cluster)) then
               CLUSTER_LL => cluster%next_cluster
            else
               NULLIFY(CLUSTER_LL)
            endif
            cDeleted = cDeleted + 1
            NULLIFY(cluster%next_cluster)
            NULLIFY(cluster%PARTICLE_LL)
            DEALLOCATE(cluster)
            NULLIFY(cluster)
         enddo

         if(cDeleted == ClusterCount) then
! dbg
!            write(*,"(4X,A,I7)")'Number of clusters deleted: ',cDeleted
            ClusterCount = 0
         else
            write(*,"()")' cDeleted /= ClusterCount'
            CALL MFIX_EXIT(myPE)
         endif

      endif

      END SUBROUTINE DELETE_CLUSTERS



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE DELETE_PARTICLES_IN_CLUSTER(cluster)

      TYPE(CLUSTER_TYPE), INTENT(INOUT), POINTER :: cluster

      TYPE(PARTICLE_TYPE), POINTER :: particle
      INTEGER pL, pDeleted

      NULLIFY(particle)
      pDeleted = 0

      if(.NOT.associated(cluster%PARTICLE_LL) .AND. cluster%ParticleCount == 0) then
         write(*,"(//,3x,A,//)") 'No particles to delete!'
      elseif(.NOT.associated(cluster%PARTICLE_LL) .AND. cluster%ParticleCount /= 0) then
         write(*,"(//,3x,A,//)") ' Error with ParticleCount and pointer - Er:1'
      elseif(associated(cluster%PARTICLE_LL) .AND. cluster%ParticleCount == 0) then
         write(*,"(//,3x,A,//)") ' Error with ParticleCount and pointer - Er:2'
      else
         do pL =1, cluster%ParticleCount
            CALL GetNextParticle(cluster, particle)
! dbg
!            write(*,"(6X,A,I7)")' Deleting particle: ',particle%ID

            if(associated(particle%next_particle)) then
               cluster%PARTICLE_LL => particle%next_particle
            else
               NULLIFY(cluster%PARTICLE_LL)
            endif
            pDeleted = pDeleted + 1
            NULLIFY(particle%next_particle)
            DEALLOCATE(particle)
            NULLIFY(particle)
         enddo

         if(pDeleted == cluster%ParticleCount) then
! dbg
!            if(pDeleted >1 ) &
!               write(*,"(4X,A,I7)")'Number of deleted particles: ',&
!               pDeleted
            cluster%ParticleCount = 0
         else
            write(*,"()")' pDeleted /= cluster%ParticleCount'
            CALL MFIX_EXIT(myPE)
         endif

      endif


      END SUBROUTINE DELETE_PARTICLES_IN_CLUSTER



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE ADD_PARTICLE_TO_CLUSTER(cluster, pID)

      TYPE(CLUSTER_TYPE), INTENT(INOUT), POINTER :: cluster
      INTEGER, INTENT(IN) :: pID

      TYPE (PARTICLE_TYPE), POINTER :: particle

! Create a particle instance.
      ALLOCATE(particle)
      NULLIFY(particle%next_particle)
! Set the particle id.
      particle%ID = pID
! Set the cluster pointer to particle.
      IF(cluster%ParticleCount == 0) THEN
         cluster%ParticleCount = 1
         cluster%PARTICLE_LL => particle
      ELSE
         cluster%ParticleCount = cluster%ParticleCount + 1
         particle%next_particle => cluster%PARTICLE_LL
         cluster%PARTICLE_LL => particle
      ENDIF


      END SUBROUTINE ADD_PARTICLE_TO_CLUSTER


!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE getNextCluster(cluster)

      TYPE(CLUSTER_TYPE), INTENT(INOUT), POINTER :: cluster

      if(.NOT.associated(cluster)) then
! point cluster to first cluster in linked list of clusters
         cluster => CLUSTER_LL
      elseif(associated(cluster%next_cluster)) then
         cluster => cluster%next_cluster
      else
         print*,' You are looking for a cluster that does not exist'
         CALL MFIX_EXIT(myPE)
      endif
      END SUBROUTINE getNextCluster



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE GetNextParticle(cluster, particle)

      TYPE(CLUSTER_TYPE), INTENT(IN), POINTER :: cluster
      TYPE(PARTICLE_TYPE), INTENT(INOUT), POINTER :: particle

      if(.NOT.associated(particle)) then
         particle => cluster%PARTICLE_LL
      elseif(associated(particle%next_particle)) then
         particle => particle%next_particle
      else
         print*,' You are looking for a particle that does not exist'
         CALL MFIX_EXIT(myPE)
      endif
      END SUBROUTINE GetNextParticle



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE ADD_PARTICLE_TO_PSEARCHHISTORY(particle, pID)

      TYPE(PARTICLE_TYPE), INTENT(OUT), POINTER :: particle
      INTEGER, INTENT(IN) :: pID

      ALLOCATE(particle)

      NULLIFY(particle%next_particle)

!dbg
!      print *, '----------> addparticle_to_psearchhistory'

      if(pSearchHistoryCount == 0) then
! no particles in search list have been created/identified
         if(associated(PSEARCH_HISTORY_LL)) then
                 print*, ' Error - particle history pointer already &
                    & associated!'
            CALL MFIX_EXIT(myPE)
         else
            pSearchHistoryCount = 1
            particle%ID = pID
! create first particle of linked list of particles. with
! psearch_history_ll always pointing to the 'first' in the list
            PSEARCH_HISTORY_LL => particle
! dbg
!            write(*,"(/A,I7)") 'Create history w/ particle ', pID
         endif
      else
         if(.NOT.associated(PSEARCH_HISTORY_LL)) then
            print*, ' Error - particle history  pointer is not &
               & associated!'
            CALL MFIX_EXIT(myPE)
         else
            pSearchHistoryCount = pSearchHistoryCount + 1
            particle%ID = pID
! establish the link between the new particle and existing particle list
            particle%next_particle => PSEARCH_HISTORY_LL
! reassign/point psearch_history_ll to be 'first' in the linked list
            PSEARCH_HISTORY_LL => particle
! as a result the linked list of particles is created from 'bottom-up'
! with the new particle  always being inserted before any existing
! particles
! dbg
!            print *, 'Add particle in history ', pID

         endif
      ENDIF

      END SUBROUTINE ADD_PARTICLE_TO_PSEARCHHISTORY



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE GetTopParticle_in_PSearchHistory(particle)

      TYPE(PARTICLE_TYPE), INTENT(INOUT), POINTER :: particle

!dbg
!     print *, '----------> gettopparticle_in_psearchhistory'

      if(.NOT.associated(particle)) then
! point particle to first particle in linked list of particles
         particle => PSEARCH_HISTORY_LL
      elseif(associated(particle%next_particle)) then
! this will never happen if we nullify particle at start of this
! routine.
         particle => particle%next_particle
      else
         print*,' You are looking for a particle that does not exist'
         CALL MFIX_EXIT(myPE)
      endif
      END SUBROUTINE GetTopParticle_in_PSearchHistory



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE DeleteTopParticle_in_PSearchHistory()

      TYPE(PARTICLE_TYPE), POINTER :: particle

      NULLIFY(particle)
!dbg
!      print *, '----------> deletetopparticle_in_psearchhistory'

      if(.NOT.associated(PSEARCH_HISTORY_LL) .AND. &
         pSearchHistoryCount == 0) then
         write(*,"(//,3x,A,//)") 'No particles in history to delete!'
      elseif(.NOT.associated(PSEARCH_HISTORY_LL) .AND. &
         pSearchHistoryCount /= 0) then
         write(*,"(//,3x,A,//)") &
            ' Error with pSearchHistoryCount and pointer - Er:1'
      elseif(associated(PSEARCH_HISTORY_LL) .AND. &
         pSearchHistoryCount == 0) then
         write(*,"(//,3x,A,//)") &
            ' Error with pSearchHistoryCount and pointer - Er:2'
      else
         CALL GetTopParticle_In_PSearchHistory(particle)
! dbg
!         write(*,"(A,I7)")' Deleting particle: ',particle%ID

         if(associated(particle%next_particle)) then
            PSEARCH_HISTORY_LL => particle%next_particle
         else
            NULLIFY(PSEARCH_HISTORY_LL)
         endif
         pSearchHistoryCount = pSearchHistoryCount - 1
         NULLIFY(particle%next_particle)
         DEALLOCATE(particle)
         NULLIFY(particle)

         if(pSearchHistoryCount < 0) then
            write(*,"(4X,A)")'pSearchHistoryCount < 0'
            CALL MFIX_EXIT(myPE)
         endif

      endif

      END SUBROUTINE DeleteTopParticle_In_PSearchHistory



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE DELETE_PSEARCHHISTORY()

      TYPE(PARTICLE_TYPE), POINTER :: particle
      INTEGER pL, pDeleted

      NULLIFY(particle)
      pDeleted = 0
!dbg
!      print *, '----------> delete_psearch_history'

      if(.NOT.associated(PSEARCH_HISTORY_LL) .AND. &
         pSearchHistoryCount == 0) then
! dbg
!         write(*,"(//,3x,A,//)") 'No particles in history to delete!'
      elseif(.NOT.associated(PSEARCH_HISTORY_LL) .AND. &
         pSearchHistoryCount /= 0) then
         write(*,"(//,3x,A,//)") &
            ' Error with pSearchHistoryCount and pointer - Er:1'
      elseif(associated(PSEARCH_HISTORY_LL) .AND. &
         pSearchHistoryCount == 0) then
         write(*,"(//,3x,A,//)") &
            ' Error with pSearchHistoryCount and pointer - Er:2'
      else
         do pL =1, pSearchHistoryCount

            CALL GetTopParticle_In_PSearchHistory(particle)
! dbg
!            write(*,"(/A,I7)")' Deleting particle: ',particle%ID

            if(associated(particle%next_particle)) then
               PSEARCH_HISTORY_LL => particle%next_particle
            else
               NULLIFY(PSEARCH_HISTORY_LL)
            endif
            pDeleted = pDeleted + 1
            NULLIFY(particle%next_particle)
            DEALLOCATE(particle)
            NULLIFY(particle)
         enddo

         if(pDeleted == pSearchHistoryCount) then
! dbg
!            write(*,"(4X,A,I7)")'Number of particles deleted: ',pDeleted
            pSearchHistoryCount = 0
         else
            write(*,"()")' pDeleted /= PSearchHistoryCount'
            CALL MFIX_EXIT(myPE)
         endif

      endif

      END SUBROUTINE DELETE_PSEARCHHISTORY



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE PRINT_CLUSTERS

      use run

      implicit none

! Generic loop counters
      integer lc1, lc2, lc3, lc4
! Particle index (mapped to cluster send/recv)
      integer L

! Flags used for debugging.
!   0: No messages
!   1: Screen only messages
!   2: Detailed Files are written by each process
      INTEGER, parameter :: dbg_level = 0

! Particle properties (for cluster bound particles only)
      double precision, dimension(:),   allocatable :: lRad
      double precision, dimension(:,:), allocatable :: lPos
      double precision, dimension(:,:), allocatable :: lVel_s

! File properticles (for cluster bound particles only)
      double precision, dimension(:),   allocatable :: lEpg
      double precision, dimension(:),   allocatable :: lRog
      double precision, dimension(:),   allocatable :: lMug
      double precision, dimension(:,:), allocatable :: lVel_g

! An array for post process data storage.
      double precision, dimension(:), allocatable :: lPost

! Cluster analysis variables:
      integer lMin(1:3), lMax(1:3)

      double precision avgEpg, avgRe, avgSlip, lSlip
      double precision avgVel_s(1:3), avgVel_g(1:3)

      double precision posMin(1:3), posMax(1:3)
      double precision clSize(1:3), clDiameter

      double precision cSize

      Type(cType), pointer :: cThis
      Type(pType), pointer :: pThis

! Merge cluster data from all processes and set up send/recv datav.
      call init_print_clusters(dbg_level)

! Get particle data from all processes. Only collect what is needed.
      call getClusterParticleData(DES_RADIUS,  lRad)   ! Radius
      call getClusterParticleData(DES_POS_NEW, lPos)   ! Position
      call getClusterParticleData(DES_VEL_NEW, lVel_s) ! Velocity

! Get field variable data from all processes.
      call getClusterFieldData(Ep_g, lEpg)
      call getClusterFieldData(RO_G, lRog)
      call getClusterFieldData(Mu_g, lMug)
      call getClusterFieldData(U_g, V_g, W_g, lVel_g)

! Allocate output data variable. (Root process only.) This data must
! be sent back to the individual processes.
      allocate(lPost( recv_sum)); lPost = zero
! Initialize the global storage variable.
      postCluster = zero

! Calculate and report cluster stats.
      if(myPE == clusterPE) then
         open (unit=203, file='clusterInfo.dat', &
            status='unknown', position='append')
         if(clusterCount_all == 0) then
            write(203,"(/3X,'Time: ',F9.6,3x,'No clusters to print.')")&
               time
         else
            write(203,"(/3X,'Time: ',F9.6)") time

            nullify(cThis)
            nullify(pThis)
            lc2 = 0
            lp_lc10: do lc1 = 1, clusterCount_all
! Only clusters with more than 3 particles are processed. Keep track of
! the number reported vs. number processed.
               if(clusters(lc1)%size < 4) cycle lp_lc10
               lc2 = lc2 + 1

! Initialize variables.
               avgEpg = zero
               avgVel_g = zero
               avgVel_s = zero
               avgRe = zero

               posMin = large_number;  lMin = 0
               posMax = zero;          lMax = 0

! Loop over particles comprising this cluster.
               lc3 = 0
               cThis => clusters(lc1)
               if(associated(cThis%particle)) then
                  pThis => cThis%particle
               else
                  nullify(pThis)
               endif
! Although this loop could be based on the cluster size, it is safer
! to loop through the link-list. If there was an allocation error, this
! should prevent the code from crashing.
               do while(associated(pThis))
                  lc3 = lc3 + 1
! Get the particle map index. This index is used to get data from the
! send/recv variables.
                  L = pThis%map

! Report a problem if viscosity is zero. (sanity check)
                  if(lMug(L) == zero) then
                     write(*,"(3x,'Invalid Mu_g. ', &
                        &'Omitting cluster data: ',I4)") lc1
                     lc2 = lc2 - 1
                     cycle lp_lc10
                  endif

                  lSlip = zero
                  do lc4=1, merge(2, 3, NO_K)
! Calculate the cluster's dimensions:
                     if((lPos(L,lc4)-lRad(L)) < posMin(lc4)) then
                        posMin(lc4) = lPos(L,lc4)-lRad(L)
                        lMin(lc4) = L
                     endif
                     if((lPos(L,lc4)+lRad(L)) > posMax(lc4)) then
                        posMax(lc4) = lPos(L,lc4)+lRad(L)
                        lMax(lc4) = L
                     endif
! Calculate the slip velocity.
                     lSlip = lSlip + (lVel_g(L,lc4)-lVel_s(L,lc4))**2
                  enddo
                  lSlip = dsqrt(lSlip)

! This is a sanity check on the data. Although a zero slip velocity is
! possible, (hopefully) there is always some small variance.
                  if(lSlip == zero) then
                     write(*,"(3x,'Invalid lSlip. ', &
                        &'Omitting cluster data: ',I4)") lc1
                     lc2 = lc2 - 1
                     cycle lp_lc10
                  endif

                  avgVel_s = avgVel_s + lVel_s(L,:)
                  avgVel_g = avgVel_g + lVel_g(L,:)
                  avgEpg = avgEpg + lEpg(L)
                  avgRe = avgRe + lRog(L)*lEpg(L)*lSlip*2.0d0*lRad(L)/lMug(L)

! Store the cluster size with the particle. (Output data in vtk files.)
                  lPost(L) = float(cThis%size)

! Increment the loop.
                  if(associated(pThis%next)) then
                     pThis => pThis%next
                  else
                     nullify(pThis)
                  endif
               enddo
! This is where the above loop is checked against the previously
! calculated cluster size. They should match.
               if(lc3 /= cThis%size) then
                  write(*,"(3x,'Error processing particles. ', &
                     &'Omitting cluster data: ',I4)") lc1
                  lc2 = lc2 - 1
                  cycle lp_lc10
               endif

               cSize = dble(cThis%size)
               avgEpg   = avgEpg   / cSize
               avgVel_g = avgVel_g / cSize
               avgVel_s = avgVel_s / cSize
               avgRe    = avgRe    / cSize

               avgSlip = zero
               clSize = zero
               do lc4=1, merge(2, 3, NO_K)
                  avgSlip = avgSlip + (avgVel_g(lc4)-avgVel_s(lc4))**2
                  clSize(lc4) = posMax(lc4) - posMin(lc4)
               enddo
               avgSlip = dsqrt(avgSlip)

! Cluster diameter is calculated by using an equivalent circle diameter
! of an ellipse in each plane. Then weighting this area by the slip
! velocity normal to it. This in effect gives a hydrodynamic diameter
! of the cluster.- 1660
               if(DO_K) then
                  clDiameter = &
                     clSize(2)*clSize(3)*(avgVel_g(1)-avgVel_s(1))**2 &
                   + clSize(1)*clSize(3)*(avgVel_g(2)-avgVel_s(2))**2 &
                   + clSize(1)*clSize(2)*(avgVel_g(3)-avgVel_s(3))**2

                  clDiameter = dsqrt(clDiameter) / avgSlip
               else
                  clDiameter = undefined
               endif

               write(203,"(4X,I7,4X,I7,4(1X,G13.6))") lc1, &
                  cThis%size, avgEpg, avgRe, avgSlip, clDiameter

            enddo lp_lc10

            write(203,"(4X,'Clusters reported: ', I7, &
               &4x,'Clusters processed: ', I7)") clusterCount_all, lc2
         endif
         close(203)
      endif

      call sendClusterData(lPost, postCluster)

      nullify(cThis)
      nullify(pThis)

      if(allocated(lRad)) deallocate(lRad)
      if(allocated(lPos)) deallocate(lPos)
      if(allocated(lVel_s)) deallocate(lVel_s)

      if(allocated(lEpg)) deallocate(lEpg)
      if(allocated(lRog)) deallocate(lRog)
      if(allocated(lMug)) deallocate(lMug)
      if(allocated(lVel_g)) deallocate(lVel_g)

      if(allocated(lPost)) deallocate(lPost)

      call finl_print_clusters(dbg_level)

      END SUBROUTINE PRINT_CLUSTERS


!......................................................................!
!  Module name: INIT_PRINT_CLUSTERS                                    !
!                                                                      !
!  Purpose: Collect cluster information from all process on root.      !
!   > Number of clusters on each process.                              !
!   > Number of particles in each cluster.                             !
!                                                                      !
!  This information is used to construct a send/receive map for        !
!  cluster data.                                                       !
!                                                                      !
!  Author: J.Musser                                   Date:  Dec-12    !
!......................................................................!
      SUBROUTINE INIT_PRINT_CLUSTERS(dbg_level)

      use run

      implicit none

! Flags used for debugging.
!   0: No messages
!   1: Screen only messages
!   2: Detailed Files are written by each process
      INTEGER, intent(in) :: dbg_level


! Counter of processes
      INTEGER proc, ierr
! Generic loop counters
      INTEGER lc1, lc2, lc3, lc4

      character(LEN=255) :: filename

! Local process data:
!----------------------->>

      TYPE(CLUSTER_TYPE), POINTER :: cluster
      TYPE(PARTICLE_TYPE), POINTER :: particle

! Number of partcies in each cluster on local process.
      INTEGER, dimension(:), allocatable :: pCnt
! Number of clusters on local process.
      INTEGER, dimension(:) :: cCnt(0:numPEs-1)
! Global IDs of particles in clusters (array)
      INTEGER, dimension(:), allocatable :: gpIDs
! IS_GHOST entries for particles in clusters. (array)
      LOGICAL, dimension(:), allocatable :: gpGhost

! Root process data:
!----------------------->>

! The number of clusters detected on each process.
      INTEGER, dimension(:) :: cCnt_all(0:numPEs-1)
! Sum of cCnt_all entries
      INTEGER cCnt_sum

! The number of particles in each cluster.
      INTEGER, dimension(:), allocatable :: pCnt_all

! Displacements used for send/receive and cluster mapping.
      INTEGER, dimension(:) :: pCnt_dsp(0:numPEs-1)

      INTEGER, dimension(:), allocatable :: gp_dsp

! Global IDs of particles in clusters (array)
      INTEGER, dimension(:), allocatable :: gpIDs_all
! Ghost entries for particles in clusters. (array)
      LOGICAL, dimension(:), allocatable :: gpGhost_all

      INTEGER, dimension(:,:), allocatable :: mergeMap
      LOGICAL merging

      Type(cType), pointer :: this

      call des_mpi_barrier()

      if(myPE == clusterPE .and. dbg_level >=1) &
         write(*,"(//1x,'Start data dump ',50('-'),'>'/)")

      if(dbg_level >= 2) then
         filename=''
         write(filename,"('dbg_check_',I2.2,'.txt')")myPE
         open(convert='big_endian',unit=202, file=filename, &
            status='replace', position='append')
         write(202,"(3x,'check A')")
      endif

! Allocate receive buffer and receive count arrays.
      allocate( recv_cnt(0:numPEs-1) )
      allocate( recv_dsp(0:numPEs-1) )

! Set the gather counts array.
      cCnt = 0; cCnt_all = 0
      cCnt(mype) = clusterCount
      call global_sum(cCnt ,cCnt_all)
      call des_mpi_barrier()
      if(dbg_level >= 1) call dbg_print_clusters(1)
      if(dbg_level >= 2) write(202,"(5x,'check 1')")

! Create an array identifying the number of particles in each cluster.
      allocate(pCnt(cCnt(myPE)))
      NULLIFY(cluster)
      do lc1 = 1, cCnt(myPE)
         CALL getNextCluster(cluster)
         pCnt(lc1) = cluster%ParticleCount
      enddo
      if(dbg_level >= 2) then
         call dbg_print_clusters(2)
         write(202,"(5x,'check 2')")
      endif

! Calculate the total number of clusters in the entire domain.
      cCnt_sum = 1
      if(myPE == clusterPE) cCnt_sum = sum(cCnt_all)

! Allocate an array sized to the total number of clusters in the
! entire domain.
      allocate(pCnt_all(cCnt_sum)); pCnt_all(:) = 1
      if(myPE == clusterPE) then
         pCnt_dsp(0) = 0
! Calculate the array displacement for MPI_GATHER
         do proc = 1,numPEs-1
            pCnt_dsp(proc) = pCnt_dsp(proc-1) + cCnt_all(proc-1)
         enddo
      endif
      if(dbg_level >= 2) then
         call dbg_print_clusters(3)
         write(202,"(5x,'check 3')")
      endif

! Populate pCnt_all on clusterPE.
      call des_mpi_gatherv(pCnt, cCnt(myPE), pCnt_all, cCnt_all,       &
         pCnt_dsp, clusterPE, ierr)
      if(dbg_level >= 1) call dbg_print_clusters(4)
      if(dbg_level >= 2) write(202,"(5x,'check 4')")
      call des_mpi_barrier()


      send_cnt = sum(pCnt)
      if(myPE == clusterPE) then
         recv_cnt(:) = 0
         do proc = 0, numPEs-1
            if(cCnt_all(proc) > 0) then
               do lc1=1, cCnt_all(proc)
                  recv_cnt(proc) = recv_cnt(proc) + &
                     pCnt_all(lc1 + pCnt_dsp(proc))
               enddo
            else
               recv_cnt(proc) = 0
            endif
         enddo
         recv_dsp(:) = 0
         do proc = 1, numPEs-1
            recv_dsp(proc) = recv_dsp(proc-1) + recv_cnt(proc-1)
         enddo
         recv_sum = sum(recv_cnt)
      else
         recv_sum = 1
      endif
      if(dbg_level >= 1) call dbg_print_clusters(5)
      if(dbg_level >= 2) write(202,"(5x,'check 5')")


      allocate( gp_dsp(cCnt_sum) ); gp_dsp(:) = 0
      if(myPE == clusterPE) then
         do lc1 = 2, cCnt_sum
            gp_dsp(lc1) = gp_dsp(lc1-1) + pCnt_all(lc1-1)
         enddo
      endif
      if(dbg_level >= 1) call dbg_print_clusters(6)
      if(dbg_level >= 2) write(202,"(5x,'check 6')")


! Transfer global particle IDs:
      call getClusterParticleData(iglobal_id, gpIDs_all, gpIDs)
      if(dbg_level >= 2) then
         call dbg_print_clusters(7)
         write(202,"(5x,'check 7')")
      endif


! Transfer particle ghost array:
      call getClusterGhostData(gpGhost_all, gpGhost)
      if(dbg_level >= 2) then
         call dbg_print_clusters(7)
         write(202,"(5x,'check 7')")
         call dbg_print_clusters(9)
         write(202,"(5x,'check 9')")
         call dbg_print_clusters(10)
         write(202,"(5x,'check 10')")
         call dbg_print_clusters(11)
         write(202,"(5x,'check 11')")
      elseif(dbg_level >= 1) then
         call dbg_print_clusters(11)
      endif


      if(myPE == clusterPE) then
! If there are two or more clusters, check to see if there is any of
! the clusters intersect. (Merge clusters that span multiple processes.)
         if( cCnt_sum > 1) then
            allocate(mergeMap(cCnt_sum, cCnt_sum)); mergeMap = 0
            do lc1 = 1, cCnt_sum
               mergeMap(lc1,lc1) = 1
            enddo
! Loop over clusters, comparing the global IDs of the particles that
! they contain. Any clusters with overlapping particles are later
! merged as a single cluster.
            lp_lc1: do lc1 = 1, cCnt_sum-1
            lp_lc2: do lc2 = 1+gp_dsp(lc1), gp_dsp(lc1)+pCnt_all(lc1)
            lp_lc3: do lc3 = lc1+1, cCnt_sum
            lp_lc4: do lc4 = 1+gp_dsp(lc3), gp_dsp(lc3)+pCnt_all(lc3)
               if(gpIDs_all(lc2) == gpIDs_all(lc4)) then
                  mergeMap(lc1,lc3) = 1
                  mergeMap(lc3,lc1) = 1
                  if(dbg_level >= 1) &
                     write(*,"(3x,'Merge: ',I5,' and ',I5)") lc1, lc3
               endif
            enddo lp_lc4
            enddo lp_lc3
            enddo lp_lc2
            enddo lp_lc1
            if(dbg_level >= 1) then
               write(*,*)''
               call dbg_print_clusters(8, lmsg='before')
            endif

! Merge any clusters that share common particles. Common particles are
! identified by their global particle IDs. Note that if two processes
! contain a common particle, the particle should be 'real' on one
! process and a ghost on the other.
            merging = .true.
            do while(merging)
               merging = .false.
               lp_lc11: do lc1=1,cCnt_sum-1
                  if(sum(mergeMap(lc1,:)) == 0) cycle lp_lc11
                  lp_lc21: do lc2=lc1+1,cCnt_sum
                     if(mergeMap(lc1,lc2) /= 1) cycle lp_lc21
                     lc3=lc2
                     lp_lc41: do lc4 = 1, cCnt_sum
                        if(lc3 == lc4) cycle lp_lc41
                        if(mergeMap(lc3,lc4) == 1) then
                           mergeMap(lc1,lc4) = 1
                           mergeMap(lc3,lc4) = 0
                           mergeMap(lc3,lc3) = 0
                           merging = .true.
                        endif
                     enddo lp_lc41
                  enddo lp_lc21
               enddo lp_lc11
            enddo
            if(dbg_level >= 1) &
               call dbg_print_clusters(8, lmsg='after')

! Count the actual number of clusters in the domain.
            clusterCount_all = 0
            do lc1=1,cCnt_sum
               if(sum(mergeMap(lc1,:)) > 0) then
                  clusterCount_all = clusterCount_all + 1
               endif
            enddo
            if(dbg_level >= 1) then
               write(*,"(3x,'Number of clusters reported by ', &
                  &'all processes: ',I6)") cCnt_sum
               write(*,"(3x,'Actual number of clusters: ',I6)") &
                  clusterCount_all
               write(*,*)''
            endif


! Allocate the clusters pointer array. This array is populated with
! data used to map particles to send/recv data. Ghost particles are not
! included in the map so that cluster analysis only contains data from
! real particles.
            allocate(clusters(clusterCount_all))
            do lc1=1,clusterCount_all
               nullify(clusters(lc1)%particle)
               clusters(lc1)%size = 0
            enddo
            lc3 = 0
            lp_lc12: do lc1=1, cCnt_sum
               if(sum(mergeMap(lc1,:)) == 0) cycle lp_lc12
               lc3 = lc3 + 1
               this => clusters(lc3); this%size = 0
               lp_lc22: do lc2=1, cCnt_sum
                  if(mergeMap(lc1,lc2) == 0) cycle lp_lc22
                  lp_lc42: do lc4=1+gp_dsp(lc2),gp_dsp(lc2)+pCnt_all(lc2)
                     if(gpGhost_all(lc4)) cycle lp_lc42
                     call addParticle(this, lc4, gpIDs_all(lc4))
                  enddo lp_lc42
               enddo lp_lc22
            enddo lp_lc12
            deallocate(mergeMap)
            if(dbg_level >= 1) &
               call dbg_print_clusters(12, dbg=dbg_level)

! If there is only one cluster, there is no need to check for
! intersecting clusters.
         elseif( cCnt_sum == 1) then
            if(dbg_level >= 1) &
               write(*,"(3x,'Just one cluster: No merge.',/)")
            clusterCount_all = 1
            allocate(clusters(clusterCount_all))
            nullify(clusters(clusterCount_all)%particle)
            clusters(clusterCount_all)%size = 0
            this => clusters(clusterCount_all); this%size = 0
            do lc1=1, pCnt_all(1)
               if(gpGhost_all(lc1)) write(*,"(3x, &
                  &'Error: Ghost particle detected: ',I7)") lc1
               call addParticle(this, lc1, gpIDs_all(lc1))
            enddo
            if(dbg_level >= 1) &
               call dbg_print_clusters(12, dbg=dbg_level)
         else
            if(dbg_level >= 1) write(*,"(3x,'No clusters to merge.')")
            clusterCount_all = 0
            allocate(clusters(1))
            nullify(clusters(1)%particle)
            clusters(1)%size = 0
         endif
      else
         clusterCount_all = 0
         allocate( clusters(1) )
         nullify(clusters(1)%particle)
         clusters(1)%size = 0
      endif
      call des_mpi_barrier()

! clean up before returning.
      nullify(cluster)
      nullify(particle)

      if(allocated(pCnt))      deallocate(pCnt)
      if(allocated(gpIDs))     deallocate(gpIDs)
      if(allocated(gpGhost))     deallocate(gpGhost)

      if(allocated(pCnt_all))  deallocate(pCnt_all)
      if(allocated(gpIDs_all)) deallocate(gpIDs_all)
      if(allocated(gpGhost_all)) deallocate(gpGhost_all)

      if(allocated(gp_dsp))    deallocate(gp_dsp)

      return

      contains

!......................................................................!
!  Module name: dbg_PRINT_CLUSTERS                                     !
!                                                                      !
!  Purpose: A collection of various write statements for debugging.    !
!                                                                      !
!  Author: J.Musser                                   Date:  Dec-12    !
!......................................................................!
      SUBROUTINE dbg_print_clusters(lmsgID, lmsg, dbg)
      use run

      implicit none

! Index for printing a specific message.
      INTEGER, INTENT(IN) :: lmsgID

! A message to write with data.
      CHARACTER(len=*), intent(in), optional :: lmsg

      INTEGER, INTENT(IN), optional :: dbg

! Generic loop counters
      INTEGER proc, lc1, lc2, lc3, lc4
! Generic write buffers.
      CHARACTER(len=120) wbuff, wbuff2
! String for generating filenames.
      CHARACTER(LEN=255) :: filename

! Generic cluster pointer.
      Type(cType), pointer :: cThis
! Generic cluster pointer.
      Type(pType), pointer :: pThis

      SELECT CASE(lmsgID)

!``````````````````````````````````````````````````````````````````````!
! Root Process Only -                                                  !
!                                                                      !
! Purpose:: Message to write out the number of clusters detected on    !
! each process.                                                        !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      CASE (1); if(myPE /= clusterPE) return

         do proc=0,numPEs-1
            write(*,"(3x,'Process ',I2,' reporting ',I4,' clusters.')") &
               proc, cCnt_all(proc)
         enddo
         write(*,"(3x,'Total reported clusters: ',I6)") sum(cCnt_all)
         write(*,*)''


!``````````````````````````````````````````````````````````````````````!
! All Processes -                                                      !
!                                                                      !
! Purpose:: Each process creates a file and records the particles      !
! contained in clusters. Only the data from the current pass is kept   !
! (file status='replace').                                             !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      CASE (2) ! All processes

         filename = ''
         write(filename,"('dbg_pCnt_',I2.2,'.txt')") myPE
         open(convert='big_endian',unit=201, file=filename, status='replace')
         write(201,"(//3x,'Time:',F18.6)")Time
         write(201,"(3x,'Number of Clusters: ',I4)") cCnt(myPE)

         if(cCnt(myPE) > 0) then
            nullify(cluster)
            do lc2 = 1, clusterCount
               write(201,"(3x,'/Cluster ',I5,' has ',I8,' members.')")
               write(201,"(5x,'Membership:')")
               CALL getNextCluster(cluster)
               NULLIFY(particle)
               do lc1 = 1, cluster%ParticleCount
                  CALL GetNextParticle(cluster, particle)
                  write(201,"(5x,'Global ID: ',I8)") &
                     iglobal_id(particle%ID)
               enddo
            enddo
         endif
         close(201)


!``````````````````````````````````````````````````````````````````````!
! Root Process Only -                                                  !
!                                                                      !
! Purpose:: Root process constructs a file that lists all the particle !
! from all processes that are contained in clusters. Only data from    !
! the current pass is kept (file status='replace').                    !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      CASE (3); if(myPE /= clusterPE) return

         filename = 'dbg_pCnt_dsp.txt'
         open(convert='big_endian',unit=201, file=filename, status='replace')
         write(201,"(/3x,'Time:',F18.6)")Time
         write(201,"(3x,'cCnt_sum: ',I4)") cCnt_sum
         do lc1 =0, numPEs-1
            write(201,"(5x,'pCnt_dsp(',I4,'): ',I6)")lc1, pCnt_dsp(lc1)
         enddo
         close(201)


!``````````````````````````````````````````````````````````````````````!
! Root Process Only -                                                  !
!                                                                      !
! Purpose:: Root process reports the number of clusters on each        !
! process and the number of particles comprising each cluster.         !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      CASE (4); if(myPE /= clusterPE) return

         if(cCnt_sum > 0)then
            write(*,"(3x,'Total number of clusters: ',I6)") cCnt_sum
            do proc = 0, numPEs-1
               if(cCnt_all(proc) > 0) then
                  write(*,"(5x,'Process ',I2,': cluster count: ',I4)") &
                     proc, cCnt_all(proc)
                  do lc1 = 1, cCnt_all(proc)
                     write(*,"(7x,'Particles in cluster ',I4,': ',I6)")&
                        lc1, pCnt_all(lc1 + pCnt_dsp(proc))
                  enddo
               else
                  write(*,"(3x,'Process ',I2,' reports no clusters.')")&
                     proc
               endif
            enddo
            write(*,*)''
         endif


!``````````````````````````````````````````````````````````````````````!
! Root Process Only -                                                  !
!                                                                      !
! Purpose:: Root process reports the information on the send/recv      !
! process for cluster_gather routines.                                 !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      CASE (5); if(myPE /= clusterPE) return

         if(recv_sum > 0)then
            write(*,"(3x,'send_cnt: ',I8)") send_cnt
            write(*,"(3x,'recv_sum: ',I8)") recv_sum
            write(*,"(3x,'Generic send/recv setup:')")
            do proc = 0, numPEs-1
               if(recv_cnt(proc) > 0) then
                  write(*,"(5x,'Process ',I2,': count: ',I6,&
                     &' disp:',I6)") &
                     proc, recv_cnt(proc), recv_dsp(proc)
               else
                  write(*,"(3x,'Process ',I2,' is empty.')") proc
               endif
            enddo
         else
            write(*,"(3x,'No data to send/recv.')")
         endif
         write(*,*)''


!``````````````````````````````````````````````````````````````````````!
! Root Process Only -                                                  !
!                                                                      !
! Purpose:: Root process reports the displacements for particle data   !
! stored in the receive buffer for cluster_gather.                     !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      CASE (6); if(myPE /= clusterPE) return

         if(cCnt_sum > 0)then
            write(*,"(3x,'Global cluster-to-particle offset:')")
            do lc1 = 1, cCnt_sum
               write(*,"(7x,' Cluster ',I4,': ',I6)") &
                  lc1, gp_dsp(lc1)
            enddo
            write(*,*)''
         endif


!``````````````````````````````````````````````````````````````````````!
! Root Process Only -                                                  !
!                                                                      !
! Purpose:: Root process constructs a file containing the global IDs   !
! of all particles contained in clusters.                              !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      CASE (7); if(myPE /= clusterPE) return

         open(convert='big_endian',unit=201, file='dbg_gpIDs.txt', status='replace')
         do lc1=1,size(gpIDs_all)
            write(201,"(5x,' Global particle ID: ',I8)") gpIDs_all(lc1)
         enddo

         close(201)


!``````````````````````````````````````````````````````````````````````!
! Root Process Only -                                                  !
!                                                                      !
! Purpose:: Root process reports cluster merging data. Output is based !
! on the array data set in call to this routine.                       !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      CASE (8); if(myPE /= clusterPE) return

         if(.not.present(lmsg)) then
            write(*,"(/3x,'Invalid use of case debug report.')")
            return
         endif

         filename = ''
         write(filename,"('dbg_mergeMap_',A,I2.2,'.txt')") trim(lmsg)
         open(convert='big_endian',unit=201, file=filename, status='replace')

         write(201,"(/3x,'Time:',F18.6)")Time

         lc1 = size(mergeMap(:,1))         ! rows
         lc2 = min(size(mergeMap(1,:)),25) ! columns

         wbuff = ''
         write(wbuff,"(7x,'1')")
         do lc3=2,lc2
            write(wbuff,"(A,3x,I1)") trim(wbuff), lc3
         enddo

         wbuff2 = ''
         write(wbuff2,"(5x,'|')")
         do lc3=1,lc2
            write(wbuff2,"(2A)") trim(wbuff2),'---|'
         enddo

         write(201,*)trim(wbuff)
         write(201,*)trim(wbuff2)

         do lc3=1,lc1
            wbuff = ''
            write(wbuff,"(3X,I1,' |')") lc3
            do lc4 = 1, lc2
               write(wbuff,"(A,1x,I1,' |')") &
                  trim(wbuff), mergeMap(lc3,lc4)
            enddo
            write(201,*)trim(wbuff)
            write(201,*)trim(wbuff2)
         enddo
         close(201)


!``````````````````````````````````````````````````````````````````````!
! All Processes -                                                      !
!                                                                      !
! Purpose:: Each process constructs a file and records ghost particle  !
! information.                                                         !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      CASE (9) ! All processes

         filename = ''
         write(filename,"('dbg_isghost_',I2.2,'.txt')") myPE
         open(convert='big_endian',unit=201, file=filename, status='replace')
         write(201,"(/3x,'Time:',F18.6)")Time
         do lc1 =1, send_cnt
            write(201,"(5x,'IS_GHOST(',I8,'): ',L2)")lc1, gpGhost(lc1)
         enddo
         close(201)


!``````````````````````````````````````````````````````````````````````!
! Root Process Only -                                                  !
!                                                                      !
! Purpose:: Root process constructs a file and records ghost particle  !
! information for all particles contained in clusters.                 !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      CASE (10); if(myPE /= clusterPE) return

         open(convert='big_endian',unit=201, file='dbg_isghost_all.txt', status='replace')
         write(201,"(/3x,'Time:',F18.6)")Time
         do lc1 =1, recv_sum
            write(201,"(5x,'IS_GHOST(',I8,'): ',L2)")lc1, gpGhost_all(lc1)
         enddo
         close(201)


!``````````````````````````````````````````````````````````````````````!
! Root Process Only -                                                  !
!                                                                      !
! Purpose:: Root process reports particles in clusters that are        !
! identified as ghost particles on a process.                          !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      CASE (11); if(myPE /= clusterPE) return

         if(recv_sum <= 0) return

         lc2 = 0
         do proc = 0, numPEs-1
         do lc1 = 1, cCnt_all(proc)
         lc2 = lc2 + 1
         do lc3 = 1, pCnt_all(lc2)
            if(gpGhost_all(gp_dsp(lc2) + lc3)) then
               write(*,"(3x,'Particle ',I8,' in cluster ',I6, &
                  &' is a ghost on process ',I2,'.')")        &
                  gpIDs_all(gp_dsp(lc2) + lc3), lc2, proc
            endif
         enddo
         enddo
         enddo
         write(*,*)''


!``````````````````````````````````````````````````````````````````````!
! Root Process Only -                                                  !
!                                                                      !
! Purpose:: Root process reports cluster statistics.                   !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      CASE (12); if(myPE /= clusterPE) return

         if(.not.present(dbg)) return

         if(dbg >= 1) then
            do lc1=1, clusterCount_all
               write(*,"(3x,'Cluster ',I6,' reports ',I6,&
                  &' particles.')") lc1, clusters(lc1)%size
            enddo
         endif

         if(dbg >= 2) then
            do lc1=1, clusterCount_all
               filename = ''
               write(filename,"('dbg_cluster_',I3.3,'.txt')")lc1
               open(convert='big_endian',unit=201, file=filename, status='replace')
               write(201,"(3x,'Time:',F10.6)") Time

               cThis => clusters(lc1)

               write(201,"(3x,'Cluster ',I6,' reports ',I6,&
                  &' particles.')") lc1, cThis%size
               write(201,"(3x,' Particles in this cluster include:')")

               nullify(pThis)
               if(associated(cThis%particle)) pThis => cThis%particle
               do while(associated(pThis))
                  write(201,"(3x,I8)") pThis%id
                  if(associated(pThis%next)) then
                     pThis => pThis%next
                  else
                     nullify(pThis)
                  endif
               enddo
               close(201)
               nullify(cThis)
            enddo
         endif


!``````````````````````````````````````````````````````````````````````!
! Root Process Only -                                                  !
!                                                                      !
! Purpose:: Default message for invalid message id.                    !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      CASE DEFAULT; if(myPE /= clusterPE) return
         WRITE(*,"(3x,'No message exists for msgID: ',I4)")lmsgID
         write(*,*)''
      END SELECT

      return
      END SUBROUTINE dbg_print_clusters

      END SUBROUTINE INIT_PRINT_CLUSTERS


!......................................................................!
!  Module name: FINL_PRINT_CLUSTERS                                    !
!                                                                      !
!  Purpose: Deallocate storage arrays used for collecting cluster data.!
!                                                                      !
!  Author: J.Musser                                   Date:  Dec-12    !
!......................................................................!
      SUBROUTINE finl_print_clusters(dbg_level)

      implicit none

! Flags used for debugging.
!   0: No messages
!   1: Screen only messages
!   2: Detailed Files are written by each process
      INTEGER, intent(in) :: dbg_level

      integer lc1
      type(cType), pointer :: cThis
      type(pType), pointer :: pThis

      if(dbg_level >= 2)  write(202,"(3x,'check B')")
      call des_mpi_barrier()

      if(allocated(recv_cnt)) deallocate(recv_cnt)
      if(allocated(recv_dsp)) deallocate(recv_dsp)

      if(mype /= clusterPE ) then
         if(allocated(clusters)) deallocate(clusters)
      else
         if(allocated(clusters)) then
            do lc1=1,size(clusters)
               cThis => clusters(lc1)
               do while(associated(cThis%particle))
                  if(associated(cThis%particle%next)) then
                     pThis => cThis%particle
                     cThis%particle => pThis%next
                     nullify(pThis%next)
                     deallocate(pThis)
                     nullify(pThis)
                  else
                     deallocate(cThis%particle)
                     nullify(cThis%particle)
                  endif
                  cThis%size = cThis%size - 1
               enddo
               if(cThis%size /= 0) then
                  write(*,"(3x,'Error deallocating clusters: ', I6)")&
                  cThis%size
               endif
            enddo
            deallocate(clusters)
         endif
      endif

! Write debug message.
      if(dbg_level >= 2) close(202)
      if(myPE == clusterPE .and. dbg_level >=1) &
         write(*,"(/1x,'End data dump ',52('-'),'<'//)")

! Clean up local pointers.
      nullify(pThis)
      nullify(cThis)
      END SUBROUTINE finl_print_clusters


!......................................................................!
!  Module name: addParticle                                            !
!                                                                      !
!  Purpose: Given a specific cluster (this) a new particle object is   !
!  is created and added to the cluster link-list.                      !
!                                                                      !
!  Author: J.Musser                                   Date:  Dec-12    !
!......................................................................!
      SUBROUTINE addParticle(this, lMap, lID)

! Cluster object getting a new partilced added.
      Type(cType), pointer, intent(inout) :: this

! The particle map is used for locating a particle in send/recv data.
      integer, intent(in) :: lMap
! The global ID of the particle. Used mainly for debugging.
      integer, intent(in) :: lID

! The new particle being added to the cluster.
      Type(pType), pointer :: new

! New particle constructor: Allocate the memory and set its entries.
      allocate(new)
      nullify(new%next)
      new%map = lMap
      new%id = lID

! Update the cluster size (the number of particles in the cluster)
      this%size = this%size + 1

! Add the new particle to the cluster.
      if(this%size == 1) then
! If this is the first particle, then the cluster should have a null
! entry for the particle link-list. (Sanity check)
         if(associated(this%particle)) then
            write(*,"(3x,'Fatal Error (000) adding particle.')")
            return
         endif
      else
! If this isn't the first particle, then the cluster should have a
! valid particle in the link list.
         if(.not.associated(this%particle)) then
            write(*,"(3x,'Fatal Error (001) adding particle.')")
            return
         endif
         new%next => this%particle
      endif

      this%particle => new
      nullify(new)

      return
      END SUBROUTINE addParticle


!......................................................................!
!  Module name: getClusterParticleData_1i                              !
!                                                                      !
!  Purpose: Collect data from local processes and send to clusterPE.   !
!  The data from cluster-bound particles is placed in a temporary      !
!  array for send/recv.                                                !
!                                                                      !
!  Incoming data is a 1D Array of Integers                             !
!                                                                      !
!  Author: J.Musser                                   Date:  Dec-12    !
!......................................................................!
      SUBROUTINE getClusterParticleData_1i(lData, lrbuff, lOut)

      implicit none

! Data on local process being sent to clusterPE.
      integer, intent(in) :: lData(:)
! Data from local processes received by clusterPE.
      integer, allocatable, intent(inout) :: lrbuff(:)
! Data on local process in cluster-array format.
      integer, allocatable, intent(inout), optional :: lOut(:)

! Temporary storage array for local data. This data is in array format.
      integer, allocatable :: lsbuff(:)
! Clusters detected by local search
      type(cluster_type),  pointer :: cluster
! Particles classified as being in a cluster on the local process.
      type(particle_type), pointer :: particle

! Generic loop counters
      integer lc1, lc2, lc3
! Error flag.
      integer ierr

! The receive buffer should be unallocated. Deallocate it otherwise.
! This array is only important on clusterPE.
      if(allocated(lrbuff)) then
         if(myPE == clusterPE) write(*,"(3x,&
            &'Error in getClusterParticleData_1i: ',&
            &'Deallocating receive buffer.')")
         deallocate(lrbuff)
      endif

! Allocate the receive buffer (and returned data set)
      allocate(lrbuff(recv_sum))
! Allocate the temporary send buffer.
      allocate(lsbuff(send_cnt))
! Allocate the local return data set
      if(present(lOut)) allocate(lOut(send_cnt))

! Populate the send buffer. The send buffer only contains data about
! particles that belong to clusters.
      if(send_cnt > 0) then
         NULLIFY(cluster)
         lc3 = 0
         do lc2 = 1, clusterCount
            CALL getNextCluster(cluster)
            NULLIFY(particle)
            do lc1 = 1, cluster%ParticleCount
               lc3 = lc3 + 1
               CALL GetNextParticle(cluster, particle)
               lsbuff(lc3) = lData(particle%ID)
            enddo
         enddo
      endif

! Invoke MPI routines.
      call des_mpi_gatherv(lsbuff, send_cnt,             &
                           lrbuff, recv_cnt, recv_dsp,   &
                           clusterPE, ierr)

! Store the local return data set.
      if(present(lOut)) lOut = lsbuff

! Clean up.
      deallocate(lsbuff)
      nullify(cluster)
      nullify(particle)

      return
      END SUBROUTINE getClusterParticleData_1i


!......................................................................!
!  Module name: getClusterParticleData_2i                              !
!                                                                      !
!  Purpose: Collect data from local processes and send to clusterPE.   !
!  The data from cluster-bound particles is placed in a temporary      !
!  array for send/recv.                                                !
!                                                                      !
!  Incoming data is a 2D Array of Integers                             !
!                                                                      !
!  Author: J.Musser                                   Date:  Dec-12    !
!......................................................................!
      SUBROUTINE getClusterParticleData_2i(lData, lrbuff, lOut)

      implicit none

! Data on local process being sent to clusterPE.
      integer, intent(in) :: lData(:,:)
! Data from local processes received by clusterPE.
      integer, allocatable, intent(inout) :: lrbuff(:,:)
! Data on local process in cluster-array format.
      integer, allocatable, intent(inout), optional :: lOut(:,:)
! Temporary storage array for local data. This data is in array format.
      integer, allocatable :: lsbuff(:)

! Clusters detected by local search
      type(cluster_type),  pointer :: cluster
! Particles classified as being in a cluster on the local process.
      type(particle_type), pointer :: particle

! Generic loop counters
      integer lc1, lc2, lc3, lc4
! Upper/Lower array bounds of lData position 2
      integer lbnd, ubnd
! Error flag.
      integer ierr

! Get the upper and lower array bounds of the incoming data.
      lbnd = lbound(lData,1)
      ubnd = ubound(lData,1)

! The receive buffer should be unallocated. Deallocate it otherwise.
! This array is only important on clusterPE.
      if(allocated(lrbuff)) then
         if(myPE == 0) write(*,"(3x,&
            &'Error in getClusterParticleData_2i: ', &
            &'Deallocating receive buffer.')")
         deallocate(lrbuff)
      endif

! Allocate the receive buffer (and returned data set)
      allocate(lrbuff(lbnd:ubnd,recv_sum))
! Allocate the temporary send buffer.
      allocate(lsbuff(send_cnt))
! Allocate the local return data set
      if(present(lOut)) allocate(lOut(lbnd:ubnd,send_cnt))

! Populate the send buffer. The send buffer only contains data about
! particles that belong to clusters.
      do lc4=lbnd,ubnd
         if(send_cnt > 0) then
            lsbuff = -987654321
            NULLIFY(cluster)
            lc3 = 0
            do lc2 = 1, clusterCount
               CALL getNextCluster(cluster)
               NULLIFY(particle)
               do lc1 = 1, cluster%ParticleCount
                  lc3 = lc3 + 1
                  CALL GetNextParticle(cluster, particle)
                  lsbuff(lc3) = lData(lc4,particle%ID)
               enddo
            enddo
         endif

! Invoke MPI routines.
         call des_mpi_gatherv(lsbuff,        send_cnt,             &
                              lrbuff(lc4, :), recv_cnt, recv_dsp,   &
                              clusterPE, ierr)
! Store the local return data set.
         if(present(lOut)) lOut(lc4,:) = lsbuff

      enddo

! Clean up.
      deallocate(lsbuff)
      nullify(cluster)
      nullify(particle)

      return
      END SUBROUTINE getClusterParticleData_2i


!......................................................................!
!  Module name: getClusterParticleData_1d                              !
!                                                                      !
!  Purpose: Collect data from local processes and send to clusterPE.   !
!  The data from cluster-bound particles is placed in a temporary      !
!  array for send/recv.                                                !
!                                                                      !
!  Incoming data is a 1D Array of double precision values.             !
!                                                                      !
!  Author: J.Musser                                   Date:  Dec-12    !
!......................................................................!
      SUBROUTINE getClusterParticleData_1d(lData, lrbuff, lOut)

      implicit none

! Data on local process being sent to clusterPE.
      double precision, intent(in) :: lData(:)
! Data from local processes received by clusterPE.
      double precision, allocatable, intent(inout) :: lrbuff(:)
! Data on local process in cluster-array format.
      double precision, allocatable, intent(inout), optional :: lOut(:)

! Temporary storage array for local data. This data is in array format.
      double precision, allocatable :: lsbuff(:)

! Clusters detected by local search
      type(cluster_type),  pointer :: cluster
! Particles classified as being in a cluster on the local process.
      type(particle_type), pointer :: particle

! Generic loop counters
      INTEGER lc1, lc2, lc3
! Error flag.
      integer ierr

! The receive buffer should be unallocated. Deallocate it otherwise.
! This array is only important on clusterPE.
      if(allocated(lrbuff)) then
         if(myPE == 0) write(*,"(3x,&
            &'Error in getClusterParticleData_1d: ', &
            &'Deallocating receive buffer.')")
         deallocate(lrbuff)
      endif

! Allocate the receive buffer (and returned data set)
      allocate(lrbuff(recv_sum))
! Allocate the temporary send buffer.
      allocate(lsbuff(send_cnt))
! Allocate the local return data set
      if(present(lOut)) allocate(lOut(send_cnt))

! Populate the send buffer. The send buffer only contains data about
! particles that belong to clusters.
      if(send_cnt > 0) then
         NULLIFY(cluster)
         lc3 = 0
         do lc2 = 1, clusterCount
            CALL getNextCluster(cluster)
            NULLIFY(particle)
            do lc1 = 1, cluster%ParticleCount
               lc3 = lc3 + 1
               CALL GetNextParticle(cluster, particle)
               lsbuff(lc3) = lData(particle%ID)
            enddo
         enddo
      endif

! Invoke MPI routines.
      call des_mpi_gatherv(lsbuff, send_cnt,             &
                           lrbuff, recv_cnt, recv_dsp,   &
                           clusterPE, ierr)

! Store the local return data set.
      if(present(lOut)) lOut = lsbuff

! Clean up.
      deallocate(lsbuff)
      nullify(cluster)
      nullify(particle)

      return
      END SUBROUTINE getClusterParticleData_1d


!......................................................................!
!  Module name: getClusterParticleData_2d                              !
!                                                                      !
!  Purpose: Collect data from local processes and send to clusterPE.   !
!  The data from cluster-bound particles is placed in a temporary      !
!  array for send/recv.                                                !
!                                                                      !
!  Incoming data is a 2D Array of double precision values.             !
!                                                                      !
!  Author: J.Musser                                   Date:  Dec-12    !
!......................................................................!
      SUBROUTINE getClusterParticleData_2d(lData, lrbuff, lOut)

      implicit none

! Data on local process being sent to clusterPE.
      double precision, intent(in) :: lData(:,:)
! Data from local processes received by clusterPE.
      double precision, allocatable, intent(inout) :: lrbuff(:,:)
! Data on local process in cluster-array format.
      double precision, allocatable, intent(inout), optional :: lOut(:,:)

! Temporary storage array for local data. This data is in array format.
      double precision, allocatable :: lsbuff(:) ! send buffer

! Clusters detected by local search
      type(cluster_type),  pointer :: cluster
! Particles classified as being in a cluster on the local process.
      type(particle_type), pointer :: particle

! Generic loop counters
      integer lc1, lc2, lc3, lc4
! Upper/Lower array bounds of lData position 2
      integer lbnd, ubnd
! Error flag.
      integer ierr

! Get the upper and lower array bounds of the incoming data.
      lbnd = lbound(lData,1)
      ubnd = ubound(lData,1)

! The receive buffer should be unallocated. Deallocate it otherwise.
! This array is only important on clusterPE.
      if(allocated(lrbuff)) then
         if(myPE == 0) write(*,"(3x,&
            &'Error in getClusterParticleData_2d: ', &
            &'Deallocating receive buffer.')")
         deallocate(lrbuff)
      endif

! Allocate the receive buffer (and returned data set)
      allocate(lrbuff(lbnd:ubnd,recv_sum))
! Allocate the temporary send buffer.
      allocate(lsbuff(send_cnt))
! Allocate the local return data set
      if(present(lOut)) allocate(lOut(lbnd:ubnd,send_cnt))

! Populate the send buffer. The send buffer only contains data about
! particles that belong to clusters.
      do lc4=lbnd,ubnd
         if(send_cnt > 0) then
            lsbuff = -9.87654321
            NULLIFY(cluster)
            lc3 = 0
            do lc2 = 1, clusterCount
               CALL getNextCluster(cluster)
               NULLIFY(particle)
               do lc1 = 1, cluster%ParticleCount
                  lc3 = lc3 + 1
                  CALL GetNextParticle(cluster, particle)
                  lsbuff(lc3) = lData(lc4,particle%ID)
               enddo
            enddo
         endif

! Invoke MPI routines.
         call des_mpi_gatherv(lsbuff,        send_cnt,             &
                              lrbuff(lc4,:), recv_cnt, recv_dsp,   &
                              clusterPE, ierr)

! Store the local return data set.
         if(present(lOut)) lOut(lc4,:) = lsbuff

      enddo

! Clean up.
      deallocate(lsbuff)
      nullify(cluster)
      nullify(particle)

      return
      END SUBROUTINE getClusterParticleData_2d


!......................................................................!
!  Module name: getClusterParticleData_1l                              !
!                                                                      !
!  Purpose: Collect data from local processes and send to clusterPE.   !
!  The data from cluster-bound particles is placed in a temporary      !
!  array for send/recv.                                                !
!                                                                      !
!  Incoming data is a 1D Array of logicals.                            !
!                                                                      !
!  Note: Logicals are converted to integers for send/recv, then        !
!  converted back to logical values.                                   !
!                                                                      !
!  Author: J.Musser                                   Date:  Dec-12    !
!......................................................................!
      SUBROUTINE getClusterGhostData(lrbuff, lOut)

      implicit none

! Data from local processes received by clusterPE.
      logical, allocatable, intent(inout) :: lrbuff(:)
! Data on local process in cluster-array format.
      logical, allocatable, intent(inout), optional :: lOut(:)

! Data on local process in cluster-array format (converted to integer)
      integer, allocatable :: lsbuff_i(:)
! Data from local processes received by clusterPE (converted to integer)
      integer, allocatable :: lrbuff_i(:)

! Clusters detected by local search
      type(cluster_type),  pointer :: cluster
! Particles classified as being in a cluster on the local process.
      type(particle_type), pointer :: particle

! Generic loop counters
      INTEGER lc1, lc2, lc3
! Error flag.
      integer ierr

! The receive buffer should be unallocated. Deallocate it otherwise.
! This array is only important on clusterPE.
      if(allocated(lrbuff)) then
         if(myPE == 0) write(*,"(3x,&
            &'Error in getClusterParticleData_1l: ',&
            &'Deallocating receive buffer.')")
         deallocate(lrbuff)
      endif

! Allocate the receive buffer (and returned data set)
      allocate(lrbuff(recv_sum)); lrbuff = .false.
! Allocate the local return data set.
      if(present(lOut)) allocate(lOut(send_cnt))

! Allocate the temporary receiver buffer.
      allocate(lrbuff_i(recv_sum)); lrbuff_i = 0
! Allocate the temporary send buffer.
      allocate(lsbuff_i(send_cnt)); lsbuff_i = 0

! Populate the send buffer. The send buffer only contains data about
! particles that belong to clusters.
      if(send_cnt > 0) then
         NULLIFY(cluster)
         lc3 = 0
         do lc2 = 1, clusterCount
            CALL getNextCluster(cluster)
            NULLIFY(particle)
            do lc1 = 1, cluster%ParticleCount
               lc3 = lc3 + 1
               CALL GetNextParticle(cluster, particle)
! Convert logical to integer.
               if(IS_GHOST(particle%ID)) lsbuff_i(lc3) = 1
            enddo
         enddo
      endif
! Invoke MPI routines.
      call des_mpi_gatherv(lsbuff_i, send_cnt,             &
                           lrbuff_i, recv_cnt, recv_dsp,   &
                           clusterPE, ierr)

! Convert integer back to logical.
      do lc1=1,recv_sum
         if(lrbuff_i(lc1) == 1) lrbuff(lc1) = .true.
      enddo

! Convert and store the local return data set.
      if(present(lOut)) then
         lOut = .false.
         do lc1 = 1,  send_cnt
            if(lsbuff_i(lc1) == 1) lOut(lc1) = .true.
         enddo
      endif

! Clean up.
      deallocate(lrbuff_i)
      deallocate(lsbuff_i)
      nullify(cluster)
      nullify(particle)

      return
    END SUBROUTINE getClusterGhostData

!......................................................................!
!  Module name: getClusterFieldData_1d                                 !
!                                                                      !
!  Purpose: Collect field (continuum) data from local processes and    !
!  send to clusterPE.  The data from the field associated with         !
!  cluster-bound particles is accessed by the particke map and not     !
!  the global IJK value.                                               !
!                                                                      !
!  Incoming data is a 1D Array of double precision values.             !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date:  Dec-12    !
!......................................................................!
      SUBROUTINE getClusterFieldData_1d(lData, lrbuff, lOut)

      implicit none

! Field data on local process being sent to clusterPE.
      double precision, intent(in) :: lData(:)
! Field data from local processes received by clusterPE.
      double precision, allocatable, intent(inout) :: lrbuff(:)
! Field ata on local process in cluster-array format.
      double precision, allocatable, intent(inout), optional :: lOut(:)

      double precision, allocatable :: lsbuff(:) ! send buffer

! Clusters detected by local search
      type(cluster_type),  pointer :: cluster
! Particles classified as being in a cluster on the local process.
      type(particle_type), pointer :: particle

! Generic loop counters
      integer lc1, lc2, lc3
! IJK index of fluid cell containing particle
      integer ijk
! Error flag.
      integer ierr

! The receive buffer should be unallocated. Deallocate it otherwise.
! This array is only important on clusterPE.
      if(allocated(lrbuff)) then
         if(myPE == 0) write(*,"(3x,&
            &'Error in getClusterFieldData_1d: ', &
            &'Deallocating receive buffer.')")
         deallocate(lrbuff)
      endif

! Allocate the receive buffer (and returned data set)
      allocate(lrbuff(recv_sum))
! Allocate the temporary send buffer.
      allocate(lsbuff(send_cnt))
! Allocate the local return data set.
      if(present(lOut)) allocate(lOut(send_cnt))

! Populate the send buffer. The send buffer only contains data about
! particles that belong to clusters.
      if(send_cnt > 0) then
         NULLIFY(cluster)
         lc3 = 0
         do lc2 = 1, clusterCount
            CALL getNextCluster(cluster)
            NULLIFY(particle)
            do lc1 = 1, cluster%ParticleCount
               lc3 = lc3 + 1
               CALL GetNextParticle(cluster, particle)
! Use the particke fluid cell IJK to access the field data.
               ijk = PIJK(particle%ID,4)
! Ghost particles have an IJK of 0. Since data from ghost particles is
! not needed, set the value as undefined.
               if(ijk == 0) then
                  lsbuff(lc3) = undefined
               else
                  lsbuff(lc3) = lData(ijk)
               endif
            enddo
         enddo
      endif

! Invoke MPI routines.
      call des_mpi_gatherv(lsbuff, send_cnt,             &
                           lrbuff, recv_cnt, recv_dsp,   &
                           clusterPE, ierr)

! Store the local return data set.
      if(present(lOut)) lOut = lsbuff

! Clean up.
      deallocate(lsbuff)
      nullify(cluster)
      nullify(particle)

      return
      END SUBROUTINE getClusterFieldData_1d


!......................................................................!
!  Module name: getClusterFieldData_3d                                 !
!                                                                      !
!  Purpose: Collect field (continuum) data from local processes and    !
!  send to clusterPE.  The data from the field associated with         !
!  cluster-bound particles is placed in a temporary and is no longer   !
!  IJK bound.                                                          !
!                                                                      !
!  Incoming data are three 1D Array of double precision values.        !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date:  Dec-12    !
!......................................................................!
      SUBROUTINE getClusterFieldData_3d(lData_1, lData_2, lData_3, &
         lrbuff, lOut)

      implicit none

! Data on local process being sent to clusterPE.
      double precision, intent(in) :: lData_1(:) ! x-axis
      double precision, intent(in) :: lData_2(:) ! y-axis
      double precision, intent(in) :: lData_3(:) ! z-axis

! Data from local processes received by clusterPE.
      double precision, allocatable, intent(inout) :: lrbuff(:,:)
! Data on local process in cluster-array format.
      double precision, allocatable, intent(inout), optional :: lOut(:,:)
! Temporary storage array for local data. This data is in array format.
      double precision, allocatable :: lsbuff(:) ! send buffer

! Clusters detected by local search
      type(cluster_type),  pointer :: cluster
! Particles classified as being in a cluster on the local process.
      type(particle_type), pointer :: particle

! Generic loop counters
      integer lc1, lc2, lc3
! IJK index of fluid cell containing particle
      integer ijk
! Error flag.
      integer ierr

! The receive buffer should be unallocated. Deallocate it otherwise.
! This array is only important on clusterPE.
      if(allocated(lrbuff)) then
         if(myPE == 0) write(*,"(3x,&
            &'Error in getClusterFieldData_3d: ', &
            &'Deallocating receive buffer.')")
         deallocate(lrbuff)
      endif

! Allocate the receive buffer (and returned data set)
      allocate(lrbuff(recv_sum,1:3))
! Allocate the temporary send buffer.
      allocate(lsbuff(send_cnt))
! Allocate the local return data set
      if(present(lOut)) allocate(lOut(send_cnt,1:3))

! Populate the send buffer. The send buffer only contains data about
! particles that belong to clusters.
      if(send_cnt > 0) then
         lsbuff = -9.87654321
         NULLIFY(cluster)
         lc3 = 0
         do lc2 = 1, clusterCount
            CALL getNextCluster(cluster)
            NULLIFY(particle)
            do lc1 = 1, cluster%ParticleCount
               lc3 = lc3 + 1
               CALL GetNextParticle(cluster, particle)
! Use the particke fluid cell IJK to access the field data.
               ijk = PIJK(particle%ID,4)
! Ghost particles have an IJK of 0. Since data from ghost particles is
! not needed, set the value as undefined.
               if(ijk == 0) then
                  lsbuff(lc3) = undefined
               else
                  lsbuff(lc3) = lData_1(ijk)
               endif
            enddo
         enddo
      endif
! Invoke MPI routines.
      call des_mpi_gatherv(lsbuff,      send_cnt,             &
                           lrbuff(:,1), recv_cnt, recv_dsp,   &
                           clusterPE, ierr)
! Store the local return data set.
      if(present(lOut)) lOut(:,1) = lsbuff

      if(send_cnt > 0) then
         lsbuff = -9.87654321
         NULLIFY(cluster)
         lc3 = 0
         do lc2 = 1, clusterCount
            CALL getNextCluster(cluster)
            NULLIFY(particle)
            do lc1 = 1, cluster%ParticleCount
               lc3 = lc3 + 1
               CALL GetNextParticle(cluster, particle)
! Use the particke fluid cell IJK to access the field data.
               ijk = PIJK(particle%ID,4)
! Ghost particles have an IJK of 0. Since data from ghost particles is
! not needed, set the value as undefined.
               if(ijk == 0) then
                  lsbuff(lc3) = undefined
               else
                  lsbuff(lc3) = lData_2(ijk)
               endif
            enddo
         enddo
      endif
! Invoke MPI routines.
      call des_mpi_gatherv(lsbuff,      send_cnt,             &
                           lrbuff(:,2), recv_cnt, recv_dsp,   &
                           clusterPE, ierr)
! Store the local return data set.
      if(present(lOut)) lOut(:,2) = lsbuff

      if(DO_K) then
         if(send_cnt > 0) then
            lsbuff = -9.87654321
            NULLIFY(cluster)
            lc3 = 0
            do lc2 = 1, clusterCount
               CALL getNextCluster(cluster)
               NULLIFY(particle)
               do lc1 = 1, cluster%ParticleCount
                  lc3 = lc3 + 1
                  CALL GetNextParticle(cluster, particle)
! Use the particke fluid cell IJK to access the field data.
                   ijk = PIJK(particle%ID,4)
! Ghost particles have an IJK of 0. Since data from ghost particles is
! not needed, set the value as undefined.
                  if(ijk == 0) then
                     lsbuff(lc3) = undefined
                  else
                     lsbuff(lc3) = lData_3(ijk)
                  endif
               enddo
            enddo
         endif
! Invoke MPI routines.
         call des_mpi_gatherv(lsbuff,      send_cnt,             &
                              lrbuff(:,3), recv_cnt, recv_dsp,   &
                              clusterPE, ierr)
! Store the local return data set.
         if(present(lOut)) lOut(:,3) = lsbuff

      endif

! Clean up.
      deallocate(lsbuff)
      nullify(cluster)
      nullify(particle)

      return
      END SUBROUTINE getClusterFieldData_3d


!......................................................................!
!  Module name: sendClusterData_1d                                     !
!                                                                      !
!  Purpose: Take cluster data from clusterPE and send it to each       !
!  process. (Storage of output data for visualization)                 !
!                                                                      !
!  Note: The local storage array is initialized to zero by default.    !
!  A different value may be used by specifying initValue.              !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date:  Dec-12    !
!......................................................................!
      SUBROUTINE sendClusterData_1d(lsbuff, lData)

      implicit none

! Data on clusterPE being distributed to the other processes.
      double precision, intent(in) :: lsbuff(:)
! Storage array on local process.
      double precision, intent(inout) :: lData(:)

! Data from clusterPE received by local process.
      double precision, allocatable :: lrbuff(:)

! Clusters detected by local search
      type(cluster_type),  pointer :: cluster
! Particles classified as being in a cluster on the local process.
      type(particle_type), pointer :: particle

! Generic loop counters
      integer lc1, lc2, lc3
! Error flag.
      integer ierr

! Allocate the temporary receive buffer on local process.
      allocate(lrbuff(send_cnt))

! Invoke MPI routines.
      call des_mpi_scatterv(lsbuff, recv_cnt, recv_dsp, lrbuff, &
         send_cnt, clusterPE, ierr)

! Map the received data back to individual particles using the local
! cluster data.
      if(send_cnt > 0) then
         NULLIFY(cluster)
         lc3 = 0
         do lc2 = 1, clusterCount
            CALL getNextCluster(cluster)
            NULLIFY(particle)
            do lc1 = 1, cluster%ParticleCount
               lc3 = lc3 + 1
               CALL GetNextParticle(cluster, particle)
               lData(particle%ID) = lrbuff(lc3)
            enddo
         enddo
      endif

! Clean up.
      deallocate(lrbuff)
      nullify(cluster)
      nullify(particle)

      return
      END SUBROUTINE sendClusterData_1d



      END MODULE DES_CLUSTER
