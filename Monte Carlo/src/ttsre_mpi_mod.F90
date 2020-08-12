module ttsre_mpi_mod

!
! this routine only needs to do something if MPI defined, but it is compiled
! in both cases.
!

  implicit none

#ifdef MPI
  include 'mpif.h'

! ... variables added for communication and log control in MPI
  character :: cfltmp*30,ext*4
  integer :: ibuf(4),ierr,iready,myid,numprocs,irc,&
       sender,status(MPI_STATUS_SIZE)
  integer :: COMM_REDUCE,myid_reduce,numprocs_reduce

! ... constants which control message tags
  integer :: tag_ready,tag_start

! ... some work variables
  integer :: icountbuf
  real,allocatable,dimension(:,:,:,:) :: Rscratch4


CONTAINS
   subroutine ttsre_mpi_init
     implicit none
    
! define a few useful things
     iready = 1
     tag_ready = 0
     tag_start = 1
   end subroutine ttsre_mpi_init
#endif
 
end module ttsre_mpi_mod
