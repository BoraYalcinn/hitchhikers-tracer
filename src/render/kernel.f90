module kernel
    use cudafor
    use iso_c_binding, only: c_int32_t
    implicit none
    contains
        attributes(global) subroutine fill_gradient(w,h,pixels)
            integer, value :: w,h
            integer(c_int32_t), device :: pixels(w,h)
            integer :: x, y, r, g, b, a

            x = (blockIdx%x - 1) * blockDim%x + threadIdx%x - 1
            y = (blockIdx%y - 1) * blockDim%y + threadIdx%y - 1

            if(x >= w .or. y >= h) then
                return
            end if

            r = x * 255 / (w - 1)
            g = y * 255 / (h - 1)
            b = 64
            a = 255
            pixels(x + 1, y + 1) = ior(shiftl(a,24), ior(shiftl(b,16), ior(shiftl(g,8), r)))
        end subroutine fill_gradient

end module kernel