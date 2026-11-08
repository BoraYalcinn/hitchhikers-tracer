module frame_output
    use cudafor
    use sdl
    use iso_c_binding, only: c_int32_t
    implicit none
    
    contains
    subroutine present_frame(pixels_d,pixels)
        integer(c_int32_t), device,intent(in) :: pixels_d(:,:)
        integer(c_int32_t), intent(inout) :: pixels(:,:)
        pixels = pixels_d        ! GPU -> CPU
        call sdl_present(pixels) ! CPU -> Screen
    end subroutine present_frame
end module frame_output
