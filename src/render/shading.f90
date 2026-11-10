module shading
    use cudafor
    use camera_rays, only: PI
    implicit none

    double precision, parameter :: cell_size = PI/12.0d0

    contains

    attributes(host,device) subroutine sky_color(d,rgb)

        double precision, intent(in) :: d(3)
        double precision, intent(out) :: rgb(3)

        double precision :: lon 
        double precision :: lat
        double precision :: z
        integer :: ix,iy

        lon = atan2(d(2),d(1))
        z = max(-1.0d0, min(1.0d0, d(3)))
        lat = asin(z)

        ix = floor(lon/cell_size)
        iy = floor(lat/cell_size)

        if (modulo(ix + iy, 2) == 0) then
            rgb(1) = 0.02d0
            rgb(2) = 0.02d0
            rgb(3) = 0.05d0
        else
            rgb(1) = 0.10d0
            rgb(2) = 0.12d0
            rgb(3) = 0.20d0
        end if

    end subroutine sky_color

end module shading  