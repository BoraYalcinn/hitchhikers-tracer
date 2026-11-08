program main
  use cudafor
  use kernel
  use sdl
  use iso_c_binding
  implicit none
  
  integer :: screen_width = 1600
  integer :: screen_height = 800
  logical :: ok = .true.
  logical :: quit = .false.

  integer :: x, y, r, g
  integer :: b = 64
  integer :: a = 255

  ! allocate host and device pixel buffers 
  integer(c_int32_t), allocatable, device :: pixels_d(:,:)
  type(dim3) :: blocks, threads
  integer(c_int32_t), allocatable :: pixels(:,:)
  allocate(pixels(screen_width, screen_height))
  allocate(pixels_d(screen_width,screen_height))

  Call sdl_open("Hitchhiker's Tracer",screen_width,screen_height,ok)
  if (ok .neqv. .true.) then
    print *,  "There was an error during the subroutine sdl_open"
    Call sdl_close
  end if

  threads = dim3(16, 16, 1)
  blocks  = dim3((screen_width + 15) / 16, (screen_height + 15) / 16, 1)
  call fill_gradient<<<blocks, threads>>>(screen_width, screen_height, pixels_d)
  pixels = pixels_d  ! GPU -> CPU copy (from host -> device)

  do while(.not. quit)
    Call sdl_poll_quit(quit)
    Call sdl_present(pixels)
  end do


  Call sdl_close

end program main
