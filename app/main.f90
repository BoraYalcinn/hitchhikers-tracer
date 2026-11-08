program main
  use sdl
  use iso_c_binding
  implicit none
  
  integer :: screen_width = 1600
  integer :: screen_height = 800
  logical :: ok = .true.
  logical :: quit = .false.

  integer :: i, j
  integer :: x, y, r, g
  integer :: b = 64
  integer :: a = 255

  integer(c_int32_t), allocatable :: pixels(:,:)
  allocate(pixels(screen_width, screen_height))

  outer_loop: do y = 0, screen_height - 1
    inner_loop: do x = 0, screen_width - 1
      r = x * 255 / (screen_width - 1)
      g = y * 255 / (screen_height - 1)
      pixels(x + 1, y + 1) = ior(shiftl(a,24), ior(shiftl(b,16), ior(shiftl(g,8), r)))
    end do inner_loop
  end do outer_loop



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
