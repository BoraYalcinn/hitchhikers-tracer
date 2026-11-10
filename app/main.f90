program main
  use cudafor
  use iso_c_binding, only: c_int32_t
  use kernel
  use sdl
  use frame_output
  use geodesic_kernel
  use camera_rays, only: degrees_to_radians
  implicit none

  integer :: screen_width = 1600
  integer :: screen_height = 800
  logical :: ok = .true.
  logical :: quit = .false.

  integer(8) :: t_start, t_now, rate 
  integer(8) :: t0,t1
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

  ! The tick rate is needed by every timing below, so read it once up front.
  call system_clock(count_rate=rate)

  fov_y = degrees_to_radians(60.0d0)
  call system_clock(t0)
  call render_frame(pixels_d, screen_width, screen_height, fov_y, r_obs, r_capture, r_escape, step_scale, max_steps)
  call system_clock(t1)
  print *, "render ms: ", 1000.0d0 * real(t1 - t0, 8) / real(rate, 8)

  call system_clock(t_start)
  frames = 0

  do while (.not. quit)
    frames = frames + 1
    call system_clock(t_now)
    if(t_now - t_start >= rate) then
      print *, "FPS: ", real(frames) / (real(t_now - t_start) / real(rate))
      t_start = t_now
      frames = 0
    end if
    call sdl_poll_quit(quit)
    call present_frame(pixels_d, pixels)
  end do

  call sdl_close()
end program main
