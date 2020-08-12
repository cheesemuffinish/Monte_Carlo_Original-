module wrimsg_mod

  implicit none
  integer :: IDVOUT,IDVERR

contains

  subroutine initialize_wrimsg_serial

    implicit none

    IDVOUT = 0
    IDVERR = 0
  end subroutine initialize_wrimsg_serial

  subroutine initialize_wrimsg_parallel

    implicit none

    IDVOUT = 12
    IDVERR = 12
  end subroutine initialize_wrimsg_parallel

end module wrimsg_mod
