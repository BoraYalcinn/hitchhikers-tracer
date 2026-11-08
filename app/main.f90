program main
  use cudafor
  use iso_c_binding, only: c_int32_t
  use kernel
  use sdl
  use frame_output
  implicit none

  integer :: screen_width = 1600
  integer :: screen_height = 800
  logical :: ok = .true.
  logical :: quit = .false.

  integer(c_int32_t), allocatable :: pixels(:,:)
  integer(c_int32_t), allocatable, device :: pixels_d(:,:)
  type(dim3) :: blocks, threads

  allocate(pixels(screen_width, screen_height))
  allocate(pixels_d(screen_width, screen_height))

  call sdl_open("Hitchhiker's Tracer", screen_width, screen_height, ok)
  if (.not. ok) then
    print *, "There was an error during the subroutine sdl_open"
    stop 1
  end if

  threads = dim3(16, 16, 1)
  blocks = dim3((screen_width + 15) / 16, (screen_height + 15) / 16, 1)
  call fill_gradient<<<blocks, threads>>>(screen_width, screen_height, pixels_d)

  do while (.not. quit)
    call sdl_poll_quit(quit)
    call present_frame(pixels_d, pixels)
  end do

  call sdl_close()
end program main