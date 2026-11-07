program main
  use sdl
  use iso_c_binding
  implicit none
  
  integer :: screen_width = 1600
  integer :: screen_height = 800
  logical :: ok = .true.
  logical :: quit = .false.

  integer(c_int32_t), allocatable :: pixels(:)
  allocate(pixels(screen_width*screen_height))

  Call sdl_open("Hitchhiker's Tracer",screen_width,screen_height,ok)
  if (ok .neqv. .true.) then
    print *,  "There was an error during the subroutine sdl_open"
    Call sdl_close
  end if

  do while(.not. quit)
    Call sdl_poll_quit(quit)
    Call sdl_present(pixels)
  end do


  Call sdl_close

end program main
