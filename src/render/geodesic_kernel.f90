module geodesic_kernel
    use cudafor
    use iso_c_binding, only: c_int32_t
    use geodesic, only: rk4_step
    use photon, only: init_ray
    use camera_rays, only: pi, pixel_direction
    use ray_plane, only: to_equatorial, from_equatorial
    use shading, only: sky_color
    implicit none

    contains

    attributes(host,device) subroutine trace_camera_ray(n_plane,r_obs,r_capture, r_escape, step_scale, max_steps, status, phi_end)

        double precision, intent(in) :: n_plane(3)
        double precision, intent(in) :: r_obs, r_capture,r_escape 
        double precision, intent(in) :: step_scale
        integer, intent(in)  :: max_steps
        integer, intent(out) :: status
        double precision, intent(out) :: phi_end

        double precision :: y(8)
        double precision :: h 
        integer :: step 

        call init_ray(r_obs, 0.5d0*pi, 0.0d0, n_plane, y)
        status = 0
        phi_end = 0.0d0

        do step = 1,max_steps
            h = step_scale * y(2)
            call rk4_step(y,h)
            if(y(2) < r_capture) then 
                status = 1
                exit
            else if (y(2) > r_escape .and. y(6) > 0.0d0) then 
                status = 2
                phi_end = y(4)
                exit
            end if
        end do

    end subroutine trace_camera_ray

    attributes(global) subroutine render_kernel(pixels,width,height,fov_y,r_obs, r_capture, r_escape,step_scale,max_steps)

    integer(c_int32_t), device :: pixels(:,:)
        integer, value :: width, height
        integer, value :: max_steps
        double precision, value :: fov_y
        double precision, value :: r_obs, r_capture, r_escape
        double precision, value :: step_scale

        integer :: i, j
        integer :: status
        integer :: red, green, blue
        double precision :: n(3), n_plane(3) 
        double precision :: tangent(2)
        double precision :: d(3)
        double precision :: rgb(3)
        double precision :: phi_end

        i = (blockIdx % x - 1) * blockDim % x + threadIdx % x
        j = (blockIdx % y - 1) * blockDim % y + threadIdx % y
        if(i > width .or. j > height) then
            return
        end if

        call pixel_direction(i-1,j-1,0.5d0,0.5d0,width,height,fov_y,n)
        call to_equatorial(n,n_plane,tangent)
        call trace_camera_ray(n_plane,r_obs,r_capture,r_escape,step_scale, max_steps, status, phi_end)
        
        if(status == 2) then
            call from_equatorial(phi_end,tangent,d)
            call sky_color(d,rgb)
        else 
            rgb = 0.0d0
        end if

        red   = int(255.0d0*rgb(1) + 0.5d0)
        green = int(255.0d0*rgb(2) + 0.5d0)
        blue  = int(255.0d0*rgb(3) + 0.5d0)

        pixels(i,j) = red + 256*green + 65536*blue - 16777216

    end subroutine render_kernel

    subroutine render_frame(pixels, width, height, fov_y, r_obs, r_capture, r_escape, step_scale, max_steps)
        integer(c_int32_t), device :: pixels(:,:)
        integer, intent(in) :: width, height 
        integer, intent(in) :: max_steps
        double precision, intent(in) :: fov_y 
        double precision, intent(in) :: r_obs, r_capture, r_escape 
        double precision, intent(in) :: step_scale

        type(dim3) :: block, grid
        integer :: ierr

        block = dim3(16,16,1)
        grid  = dim3((width + 15)/16,(height + 15)/16,1)
        call render_kernel<<<grid,block>>>(pixels, width, height, fov_y, r_obs,r_capture, r_escape, step_scale, max_steps)
        ierr = cudaDeviceSynchronize()

    end subroutine render_frame

end module geodesic_kernel