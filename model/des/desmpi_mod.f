!----------------------------------------------------------------------!
!  Module: DESMPI                                                      !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
!  Purpose: Contains routines for packing real and ghost particles     !
!     into the MPI send buffers.                                       !
!----------------------------------------------------------------------!
      module desmpi

! Ghost particle packet size.
      INTEGER :: iGhostPacketSize
      INTEGER :: iParticlePacketSize
      INTEGER :: iPairPacketSize

! Flags and constants for interfaces
      integer, dimension(:), allocatable :: ineighproc
      logical, dimension(:), allocatable :: iexchflag

! offset for periodic boundaries
      double precision, dimension(:,:), allocatable :: dcycl_offset

      type array
         double precision, dimension(:), allocatable :: facebuf
      end type array

! following variables used for sendrecv ghost particles and particle exchange
      type(array), dimension(:), allocatable :: dsendbuf
      type(array), dimension(:), allocatable :: drecvbuf

      integer,dimension(:),allocatable:: isendcnt
      integer,dimension(:),allocatable:: isendreq
      integer,dimension(:),allocatable:: irecvreq

      integer,parameter :: ibufoffset = 2

! The maximum size of the receive buffer.
      integer :: imaxbuf
      integer :: ispot

! following variables are used for gather and scatter
      double precision, dimension(:), allocatable :: drootbuf
      double precision, dimension(:), allocatable :: dprocbuf
      integer, dimension(:), allocatable :: irootbuf
      integer, dimension(:), allocatable :: iprocbuf

      integer,dimension(:), allocatable:: idispls
      integer,dimension(:), allocatable:: iscattercnts
      integer,dimension(:), allocatable:: igathercnts

      integer :: iscr_recvcnt
      integer :: igath_sendcnt

! following variables are used to identify the cell number for ghost cells
      integer,dimension(:,:),allocatable :: isendindices
      integer,dimension(:,:),allocatable :: irecvindices

! variables used to read initial particle properties
      double precision, dimension(:,:), allocatable:: dpar_pos
      double precision, dimension(:,:), allocatable:: dpar_vel
      double precision, dimension(:), allocatable:: dpar_den
      double precision, dimension(:), allocatable:: dpar_rad

      contains

!------------------------------------------------------------------------
! subroutine       : des_dbgmpi
! Purpose          : For printing the flags and values set for interface
!                    communication
! Parameters       : ptype - based on this following info is printed to
!                    the file
!                    1 - interface flags
!                    2 - send buffer for ghost particles
!                    3 - recv buffer for ghost particles
!                    4 - particle information
!                    5 - send buffer for particles exchanging processor
!                    6 - particles info
!                    7 - neighinfo
!------------------------------------------------------------------------
      subroutine des_dbgmpi(ptype)

      use discretelement, only: DES_POS_NEW
      use discretelement, only: iGLOBAL_ID

      use discretelement, only: S_TIME
      use discretelement, only: DIMN
      use discretelement, only: DO_NSEARCH
      use discretelement, only: iGHOST_CNT
      use discretelement, only: MAX_PIP, PIP
      use functions, only: is_ghost, is_nonexistent, is_normal, is_entering_ghost, is_exiting_ghost

      use geometry, only: NO_K
      use compar, only: myPE

      use desgrid, only: dg_funijk, iofpos, jofpos

!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer ptype
!-----------------------------------------------
! local varaiables
!-----------------------------------------------
      character (255) filename
      integer lcurpar,lpacketsize,lface,lparcnt,lbuf,lindx,ltordimn
      integer lneighcnt,lneighindx
      integer lsize
      double precision xpos,ypos
      integer li,lj,lparcount
!-----------------------------------------------

      write(filename,'("dbg_desmpi",I4.4,".dat")') mype
      open(44,file=filename,convert='big_endian')
      select case(ptype)
      case (1)
         write(44,*)&
            "------------------------------------------------------"
         write(44,*) "Flag Information"
         do lface =1,dimn*2
            write(44,*) "details for face =" , lface
            write(44,*) "Exchflag, cyclfac, neighproc" ,iexchflag(lface),ineighproc(lface)
         end do
         write(44,*) &
            "------------------------------------------------------"
      case (2)
         ltordimn = merge(1,3,NO_K)
         lpacketsize = 2*dimn + ltordimn+ 5
         do lface =1,dimn*2
            if (.not.iexchflag(lface))cycle
            lparcnt = dsendbuf(1+mod(lface,2))%facebuf(1)
            if (lparcnt .gt. 0) then
               write(44,*) &
                "------------------------------------------------------"
               write(44,*) "ghost send buffer for face", lface
               write(44,*) "Number of particles in sendbuf",lparcnt
               write(44,*) "particle number global_id ijk prvijk ",&
                  "radius material new_pos new_vel omega_new"
               write(44,*) &
                "-----------------------------------------------------"
               do lcurpar = 1,lparcnt
                  lbuf = (lcurpar-1) * lpacketsize + ibufoffset
                  write(44,*) lcurpar,(dsendbuf(1+mod(lface,2))%facebuf(lindx),lindx=lbuf,lbuf+lpacketsize-1)
               end do
            end if
         end do
      case (3)
         ltordimn = merge(1,3,NO_K)
         lpacketsize = 2*dimn + ltordimn+ 5
         do lface =1,dimn*2
            if (.not.iexchflag(lface))cycle
            lparcnt = drecvbuf(1+mod(lface,2))%facebuf(1)
            if (lparcnt .gt. 0) then
               write(44,*) &
                "------------------------------------------------------"
               write(44,*) "ghost recv buffer for face", lface
               write(44,*) "Number of particles in recvbuf",lparcnt
               write(44,*) "particle number global_id ijk prvijk ",&
                  "radius material new_pos new_vel omega_new"
               write(44,*) &
                 "-----------------------------------------------------"
               do lcurpar = 1,lparcnt
                  lbuf = (lcurpar-1) * lpacketsize + ibufoffset
                  write(44,*) lcurpar,(drecvbuf(1+mod(lface,2))%facebuf(lindx),lindx=lbuf,lbuf+lpacketsize-1)
               end do
            end if
         end do
      case (4)
          write(44,*) &
             "---------------------------------------------------------"
          write(44,*) "Particle info"
          write(44,*) "max_pip,pip =" , max_pip,pip
          write(44,*) "ghost position                        ",&
             "i       j     k    ijk"
          write(44,*) &
             "---------------------------------------------------------"
          lparcount = 1
          do lcurpar=1,max_pip
             if (lparcount.gt.pip) exit
             if (is_nonexistent(lcurpar))cycle
             lparcount=lparcount + 1
             xpos = des_pos_new(lcurpar,1)
             ypos = des_pos_new(lcurpar,2)
             li=iofpos(xpos);lj=jofpos(ypos)
             write(44,*)(is_ghost(lcurpar).or.is_entering_ghost(lcurpar).or.is_exiting_ghost(lcurpar)),xpos,ypos,li,lj,dg_funijk(li,lj,1)
          end do
      case (5)
         ltordimn = merge(1,3,NO_K)
         lpacketsize = 9*dimn + ltordimn*4 + 13
         do lface =1,dimn*2
            if (.not.iexchflag(lface))cycle
            lparcnt = dsendbuf(1+mod(lface,2))%facebuf(1)
            if (lparcnt .gt. 0) then
               write(44,*) &
                "------------------------------------------------------"
               write(44,*) "particle crossing info send buffer", lface
               write(44,*) "Number of particles in sendbuf",lparcnt
               do lcurpar = 1,lparcnt
                  lbuf = (lcurpar-1) * lpacketsize + ibufoffset
                  write(44,*) "global_id  ijk prvijk radius  i,j,k, ijk"
                  write(44,*) &
                 "-----------------------------------------------------"
                  lsize = 8
                  write(44,'(5(2x,f8.4))') (dsendbuf(1+mod(lface,2))%facebuf(lindx),lindx=lbuf,lbuf+lsize-1)
                  lbuf = lbuf + lsize

                  write(44,*) "phase density vol mass omoi pos_old"
                  write(44,*) &
                 "-----------------------------------------------------"
                  lsize = 5+dimn
                  write(44,'(5(2x,f8.4))') (dsendbuf(1+mod(lface,2))%facebuf(lindx),lindx=lbuf,lbuf+lsize-1)
                  lbuf = lbuf + lsize

                  write(44,*) "pos_new     vel_old   vel_new"
                  write(44,*) &
                 "-----------------------------------------------------"
                  lsize = 3*dimn
                  write(44,'(5(2x,f8.4))') (dsendbuf(1+mod(lface,2))%facebuf(lindx),lindx=lbuf,lbuf+lsize-1)
                  lbuf = lbuf + lsize

                  write(44,*) "omega_old     omega_new"
                  write(44,*) &
                 "-----------------------------------------------------"
                  lsize = ltordimn*2
                  write(44,'(5(2x,f8.4))') (dsendbuf(1+mod(lface,2))%facebuf(lindx),lindx=lbuf,lbuf+lsize-1)
                  lbuf = lbuf + lsize

                  write(44,*) "acc_old     rot_acc_old   fc "
                  write(44,*) &
                 "-----------------------------------------------------"
                  lsize = 2*dimn + ltordimn
                  write(44,'(5(2x,f8.4))') (dsendbuf(1+mod(lface,2))%facebuf(lindx),lindx=lbuf,lbuf+lsize-1)
                  lbuf = lbuf + lsize

                  write(44,*) "fn ft tow"
                  write(44,*) &
                 "-----------------------------------------------------"
                  lsize = 2*dimn + ltordimn
                  write(44,'(5(2x,f8.4))') (dsendbuf(1+mod(lface,2))%facebuf(lindx),lindx=lbuf,lbuf+lsize-1)
                  lbuf = lbuf + lsize

! print neighbour information
                  lneighcnt =dsendbuf(1+mod(lface,2))%facebuf(lbuf);lbuf = lbuf + 1
                  write(44,*) "total neighbour=",lneighcnt
                  write(44,*) "neighbou",lneighcnt
                  do lneighindx = 1, lneighcnt
                     lsize = 3
                     write(44,'(5(2x,f8.4))') (dsendbuf(1+mod(lface,2))%facebuf(lindx),lindx=lbuf,lbuf+lsize-1)
                     lbuf = lbuf + lsize
                  enddo
               enddo
            endif
         enddo
      case (6)
         write(44,*) "-----------------------------------------------"
         write(44,*) "at Time =",s_time
         write(44,*) "Total paticles =",pip
         write(44,*) "Total ghost paticles =",ighost_cnt
         write(44,*) "do_nsearch =",do_nsearch
         lparcnt = 1
         do lcurpar = 1,max_pip
            if(lparcnt.gt.pip) exit
            lparcnt = lparcnt + 1
            write(44,*) "particle position =",des_pos_new(lcurpar,1:dimn)
         end do
         write(44,*) "-----------------------------------------------"
      case (7)
         write(44,*) "-----------------------------------------------"
         write(44,*) "pip and max_pip" , pip, max_pip,.not.is_nonexistent(1)
         write(44,*) s_time
         lparcnt = 1
         do lcurpar =1,max_pip
            if(lparcnt.gt.pip) exit
            if(is_nonexistent(lcurpar)) cycle
            lparcnt = lparcnt+1
            if(is_ghost(lcurpar).or.is_entering_ghost(lcurpar).or.is_exiting_ghost(lcurpar)) cycle
            write(44,*) "Info for particle", iglobal_id(lcurpar)
            write(44,*) "position new ", des_pos_new(lcurpar,:)
         end do
         write(44,*) "-----------------------------------------------"
      end select
      close(44)
      end subroutine des_dbgmpi


      end module
