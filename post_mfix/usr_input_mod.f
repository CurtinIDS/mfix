!
! *********************************************************************
! ************************** module : usr_input ***********************
! *********************************************************************
!
!
      module usr_input

         use param1

         real :: usr_t1              ! start_time
         real :: usr_t2              ! end time
         logical :: usr_tavg         ! time average flag
         logical :: usr_allow_tavg   ! flag to allow/not time averaging

         logical :: usr_done         ! true if no more inputs desired

         character(len=8) :: usr_var      ! variable name
         integer   :: usr_var_num    ! index corresponding to usr_var

         integer   :: usr_m          ! solid phase index
         integer   :: usr_n          ! species index

         integer   :: usr_direction

         logical   :: usr_sum

         integer   :: usr_i1
         integer   :: usr_i2
         logical   :: usr_i_avg
         logical   :: usr_must_i_avg

         integer   :: usr_j1
         integer   :: usr_j2
         logical   :: usr_j_avg
         logical   :: usr_must_j_avg

         integer   :: usr_k1
         integer   :: usr_k2
         logical   :: usr_k_avg
         logical   :: usr_must_k_avg

         integer   :: usr_status

         character(len=120) :: usr_fname
         logical   :: usr_ask_for_file


         real :: usr_t1_a              ! start_time
         real :: usr_t2_a              ! end time
         logical :: usr_tavg_a         ! time average flag
         logical :: usr_allow_tavg_a   ! flag to allow/not time averaging

         character(len=8) :: usr_var_a      ! variable name
         integer   :: usr_var_num_a    ! index corresponding to usr_var

         integer   :: usr_m_a          ! solid phase index
         integer   :: usr_n_a          ! species index

         integer   :: usr_direction_a

         logical   :: usr_sum_a

         integer   :: usr_i1_a
         integer   :: usr_i2_a
         logical   :: usr_i_avg_a
         logical   :: usr_must_i_avg_a

         integer   :: usr_j1_a
         integer   :: usr_j2_a
         logical   :: usr_j_avg_a
         logical   :: usr_must_j_avg_a

         integer   :: usr_k1_a
         integer   :: usr_k2_a
         logical   :: usr_k_avg_a
         logical   :: usr_must_k_avg_a



         INTEGER   :: REC_POINTER(N_SPX)
         LOGICAL   :: READ_SPX(N_SPX) , AT_EOF(N_SPX)
         real      :: time_real(n_spx)


      end module usr_input
