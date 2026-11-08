! sdl.f90 -- minimal SDL2 bindings 

! Fortran is case-insensitive, so constants must not share a name with an SDL
! function (SDL_QUIT would clash with SDL_Quit). Constants here use short
! names without the SDL_ prefix.

module sdl
  use, intrinsic :: iso_c_binding
  implicit none
  private

  public :: sdl_open, sdl_close, sdl_poll_quit, sdl_present, sdl_delay_ms

  ! SDL constants
  integer(c_int32_t), parameter :: INIT_VIDEO         = int(z'00000020', c_int32_t)
  integer(c_int),     parameter :: WINDOWPOS_CENTERED = int(z'2FFF0000', c_int)
  integer(c_int32_t), parameter :: WINDOW_SHOWN       = int(z'00000004', c_int32_t)
  integer(c_int32_t), parameter :: RENDERER_ACCEL     = int(z'00000002', c_int32_t)
  integer(c_int32_t), parameter :: RENDERER_VSYNC     = int(z'00000004', c_int32_t)
  ! SDL_PIXELFORMAT_ABGR8888 == bytes R,G,B,A in memory on little-endian
  integer(c_int32_t), parameter :: PIXFMT_RGBA        = int(z'16762004', c_int32_t)
  integer(c_int),     parameter :: TEXACCESS_STREAM   = 1_c_int
  integer(c_int32_t), parameter :: EV_QUIT            = int(z'00000100', c_int32_t)

  interface
    function SDL_Init(flags) bind(C, name="SDL_Init")
      import :: c_int, c_int32_t
      integer(c_int32_t), value :: flags
      integer(c_int) :: SDL_Init
    end function

    subroutine SDL_Quit() bind(C, name="SDL_Quit")
    end subroutine

    function SDL_GetError() bind(C, name="SDL_GetError")
      import :: c_ptr
      type(c_ptr) :: SDL_GetError
    end function

    function SDL_CreateWindow(title, x, y, w, h, flags) bind(C, name="SDL_CreateWindow")
      import :: c_ptr, c_char, c_int, c_int32_t
      character(kind=c_char), dimension(*), intent(in) :: title
      integer(c_int), value :: x, y, w, h
      integer(c_int32_t), value :: flags
      type(c_ptr) :: SDL_CreateWindow
    end function

    subroutine SDL_DestroyWindow(window) bind(C, name="SDL_DestroyWindow")
      import :: c_ptr
      type(c_ptr), value :: window
    end subroutine

    function SDL_CreateRenderer(window, index, flags) bind(C, name="SDL_CreateRenderer")
      import :: c_ptr, c_int, c_int32_t
      type(c_ptr), value :: window
      integer(c_int), value :: index
      integer(c_int32_t), value :: flags
      type(c_ptr) :: SDL_CreateRenderer
    end function

    subroutine SDL_DestroyRenderer(renderer) bind(C, name="SDL_DestroyRenderer")
      import :: c_ptr
      type(c_ptr), value :: renderer
    end subroutine

    function SDL_CreateTexture(renderer, format, access, w, h) bind(C, name="SDL_CreateTexture")
      import :: c_ptr, c_int, c_int32_t
      type(c_ptr), value :: renderer
      integer(c_int32_t), value :: format
      integer(c_int), value :: access, w, h
      type(c_ptr) :: SDL_CreateTexture
    end function

    subroutine SDL_DestroyTexture(texture) bind(C, name="SDL_DestroyTexture")
      import :: c_ptr
      type(c_ptr), value :: texture
    end subroutine

    function SDL_UpdateTexture(texture, rect, pixels, pitch) bind(C, name="SDL_UpdateTexture")
      import :: c_ptr, c_int
      type(c_ptr), value :: texture, rect, pixels
      integer(c_int), value :: pitch
      integer(c_int) :: SDL_UpdateTexture
    end function

    function SDL_RenderClear(renderer) bind(C, name="SDL_RenderClear")
      import :: c_ptr, c_int
      type(c_ptr), value :: renderer
      integer(c_int) :: SDL_RenderClear
    end function

    function SDL_RenderCopy(renderer, texture, srcrect, dstrect) bind(C, name="SDL_RenderCopy")
      import :: c_ptr, c_int
      type(c_ptr), value :: renderer, texture, srcrect, dstrect
      integer(c_int) :: SDL_RenderCopy
    end function

    subroutine SDL_RenderPresent(renderer) bind(C, name="SDL_RenderPresent")
      import :: c_ptr
      type(c_ptr), value :: renderer
    end subroutine

    function SDL_PollEvent(event) bind(C, name="SDL_PollEvent")
      import :: c_int, c_int8_t
      ! SDL_Event is a 56-byte union; we only read the first 4 bytes (which is the type).
      integer(c_int8_t), intent(out) :: event(*)
      integer(c_int) :: SDL_PollEvent
    end function

    subroutine SDL_Delay(ms) bind(C, name="SDL_Delay")
      import :: c_int32_t
      integer(c_int32_t), value :: ms
    end subroutine
  end interface

  type(c_ptr), save :: window   = c_null_ptr
  type(c_ptr), save :: renderer = c_null_ptr
  type(c_ptr), save :: texture  = c_null_ptr
  integer,     save :: tex_w = 0

contains

  ! Open a window of w x h pixels with a vsynced renderer and one streaming
  ! texture of the same size. Returns ok = .false. on any failure.
  subroutine sdl_open(title, w, h, ok)
    character(len=*), intent(in)  :: title
    integer,          intent(in)  :: w, h
    logical,          intent(out) :: ok

    ok = .false.
    if (SDL_Init(INIT_VIDEO) /= 0) then
      call report('SDL_Init')
      return
    end if

    window = SDL_CreateWindow(trim(title)//c_null_char, WINDOWPOS_CENTERED, &
                              WINDOWPOS_CENTERED, int(w, c_int), int(h, c_int), WINDOW_SHOWN)
    if (.not. c_associated(window)) then
      call report('SDL_CreateWindow')
      call SDL_Quit()
      return
    end if

    renderer = SDL_CreateRenderer(window, -1_c_int, ior(RENDERER_ACCEL, RENDERER_VSYNC))
    if (.not. c_associated(renderer)) then
      call report('SDL_CreateRenderer')
      call sdl_close()
      return
    end if

    texture = SDL_CreateTexture(renderer, PIXFMT_RGBA, TEXACCESS_STREAM, &
                                int(w, c_int), int(h, c_int))
    if (.not. c_associated(texture)) then
      call report('SDL_CreateTexture')
      call sdl_close()
      return
    end if

    tex_w = w
    ok = .true.
  end subroutine sdl_open

  ! Drain the event queue. Sets quit = .true. if the window was closed.
  subroutine sdl_poll_quit(quit)
    logical, intent(out) :: quit
    integer(c_int8_t) :: buf(64)
    integer(c_int32_t) :: etype

    quit = .false.
    do while (SDL_PollEvent(buf) /= 0)
      etype = transfer(buf(1:4), etype)
      if (etype == EV_QUIT) quit = .true.
    end do
  end subroutine sdl_poll_quit

  ! Show one finished frame. pixels holds w*h RGBA pixels, one 32-bit integer
  ! each (byte order R,G,B,A in memory), row-major, top row first.
  subroutine sdl_present(pixels)
    integer(c_int32_t), intent(in), target, contiguous :: pixels(:,:)
    integer(c_int) :: rc

    rc = SDL_UpdateTexture(texture, c_null_ptr, c_loc(pixels), int(4 * tex_w, c_int))
    rc = SDL_RenderClear(renderer)
    rc = SDL_RenderCopy(renderer, texture, c_null_ptr, c_null_ptr)
    call SDL_RenderPresent(renderer)   ! blocks until vsync
  end subroutine sdl_present

  subroutine sdl_delay_ms(ms)
    integer, intent(in) :: ms
    call SDL_Delay(int(ms, c_int32_t))
  end subroutine sdl_delay_ms

  subroutine sdl_close()
    if (c_associated(texture))  call SDL_DestroyTexture(texture)
    if (c_associated(renderer)) call SDL_DestroyRenderer(renderer)
    if (c_associated(window))   call SDL_DestroyWindow(window)
    texture  = c_null_ptr
    renderer = c_null_ptr
    window   = c_null_ptr
    call SDL_Quit()
  end subroutine sdl_close

  subroutine report(what)
    character(len=*), intent(in) :: what
    character(kind=c_char), pointer :: msg(:)
    type(c_ptr) :: p
    integer :: n

    p = SDL_GetError()
    print '(a,a)', what, ' failed: '
    if (c_associated(p)) then
      call c_f_pointer(p, msg, [256])
      n = 0
      do while (n < 256)
        if (msg(n + 1) == c_null_char) exit
        n = n + 1
      end do
      print '(256a)', msg(1:n)
    end if
  end subroutine report

end module sdl