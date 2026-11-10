program main
  use cudafor
  use iso_c_binding, only: c_int32_t
  use kernel
  use sdlü
  use frame_output
  use geodesic_kernel
  use camera_rays, only: degrees_to_radians
  implicit none

  integer :: screen_width = 1600
  integer :: screen_height = 800
  logical :: ok = .true.
  logical :: quit = .false.

  integer(8) :: t_start, t_now, rate 
  integer :: frames

  integer(c_int32_t), allocatable :: pixels(:,:)
  integer(c_int32_t), allocatable, device :: pixels_d(:,:)
  double precision, parameter :: r_obs = 20.0d0
  double precision, parameter :: r_capture = 2.05d0
  double precision, parameter :: r_escape = 100.0d0
  double precision, parameter :: step_scale = 0.005d0
  integer, parameter :: max_steps = 200000
  double precision :: fov_y

  allocate(pixels(screen_width, screen_height))
  allocate(pixels_d(screen_width, screen_height))

  call sdl_open("Hitchhiker's Tracer", screen_width, screen_height, ok)
  if (.not. ok) then
    print *, "There was an error during the subroutine sdl_open"
    stop 1
  end if

  fov_y = degrees_to_radians(60.0d0)
  call render_frame(pixels_d, screen_width, screen_height, fov_y, r_obs, r_capture, r_escape, step_scale, max_steps)
    
  call system_clock(count_rate=rate)
  call system_clock(t_start)
  frames = 0

  do while (.not. quit)
    frames = frames + 1
    call system_clock(t_now)
    if(t_now - t_start >= rate) then
      print *, "FPS: ", real(frames) / (real(t_now - t_start) / real(rate))
    end if
    call sdl_poll_quit(quit)
    call present_frame(pixels_d, pixels)
  end do

  call sdl_close()
end program main