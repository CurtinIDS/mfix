      MODULE time_cpu

! cpu time/second
      DOUBLE PRECISION :: CPUos
! old cpu time and time for calculating CPUos
      DOUBLE PRECISION :: CPU_NLOG, TIME_NLOG
! Initial value of CPU time.
      DOUBLE PRECISION :: CPU0
! Time for IO
      DOUBLE PRECISION :: CPU_IO = 0.0d0

! Initial value of CPU time at the begin of MFIX, prior any I/O
      DOUBLE PRECISION :: CPU00
      DOUBLE PRECISION :: WALL0

! Time at start of simulation
      DOUBLE PRECISION :: TIME_START
! Wall time at the beginning
      DOUBLE PRECISION :: WALL_START

      END MODULE time_cpu
