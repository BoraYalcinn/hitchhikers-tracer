module camera_rays
    use cudafor
    implicit none

    double precision, parameter :: PI = 4.0d0*atan(1.0d0)
    
    contains

    attributes(host,device) function degrees_to_radians(degrees) result(radians)
        double precision, intent(in) :: degrees
        double precision :: radians
 
        radians = degrees * PI / 180.0d0
    end function degrees_to_radians

    attributes(host,device) subroutine pixel_direction(i,j,ox,oy,width,height,fov_y,n)
        integer, intent(in) :: i 
        integer, intent(in) :: j
        double precision, intent(in) :: ox
        double precision, intent(in) :: oy 
        integer, intent(in) :: width
        integer, intent(in) :: height
        double precision, intent(in) :: fov_y
        double precision, intent(out) :: n(3)  
        
        double precision :: tan_half_fov
        double precision :: aspect_ratio
        double precision :: x_ndc, y_ndc
        double precision :: u,v
        double precision :: norm_n

        tan_half_fov = tan(0.5d0*fov_y)
        
        aspect_ratio = dble(width)/dble(height)
 
        ! Pixel position -> [-1,1] screen coordinates (-1 left/bottom, +1 right/top).
        ! Rows grow downward, hence the flipped sign for y.
        x_ndc = 2.0d0*(dble(i) + ox)/dble(width) - 1.0d0
        y_ndc = 1.0d0 - 2.0d0*(dble(j) + oy)/dble(height)
 
        ! Screen coordinates -> position on the image plane.
        u = x_ndc * aspect_ratio * tan_half_fov
        v = y_ndc * tan_half_fov
 
        ! Camera looks along -r; up is -theta, right is +phi.
        n(1) = -1.0d0
        n(2) = -v
        n(3) = u
 
        norm_n = sqrt(n(1)**2 + n(2)**2 + n(3)**2)
        n = n/norm_n

    end subroutine pixel_direction


end module camera_rays