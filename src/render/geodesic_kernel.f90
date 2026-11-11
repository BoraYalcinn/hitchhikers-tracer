module geodesic_kernel
    use cudafor
    use iso_c_binding, only: c_int32_t
    use geodesic, only: rk4_step, rk4_step_sp
    use photon, only: init_ray, init_ray_sp
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
                phi_end = y(4) + asin( (y(8) / (-y(5))) / y(2) )
                exit
            end if
        end do

    end subroutine trace_camera_ray

    ! Single precision copy of trace_camera_ray. This is what render_kernel
    ! runs: consumer GPUs do double precision at 1/64 of the single
    ! precision rate, and the picture does not need more than float
    ! accuracy. trace_camera_ray above stays the double precision
    ! reference; test_geodesic_kernel checks the two agree.
    attributes(host,device) subroutine trace_camera_ray_sp(n_plane,r_obs,r_capture, r_escape, step_scale, max_steps, status, phi_end)

        real(4), intent(in) :: n_plane(3)
        real(4), intent(in) :: r_obs, r_capture,r_escape
        real(4), intent(in) :: step_scale
        integer, intent(in)  :: max_steps
        integer, intent(out) :: status
        real(4), intent(out) :: phi_end

        real(4) :: y(8)
        real(4) :: h
        integer :: step

        call init_ray_sp(r_obs, 0.5e0*real(pi,4), 0.0e0, n_plane, y)
        status = 0
        phi_end = 0.0e0

        do step = 1,max_steps
            h = step_scale * y(2)
            call rk4_step_sp(y,h)
            if(y(2) < r_capture) then
                status = 1
                exit
            else if (y(2) > r_escape .and. y(6) > 0.0e0) then
                status = 2
                phi_end = y(4) + asin( (y(8) / (-y(5))) / y(2) )
                exit
            end if
        end do

    end subroutine trace_camera_ray_sp

    attributes(global) subroutine render_kernel(pixels,width,height,fov_y,r_obs, r_capture, r_escape,step_scale,max_steps,samples)

        integer(c_int32_t), device :: pixels(:,:)
        integer, value :: width, height
        integer, value :: max_steps
        double precision, value :: fov_y
        double precision, value :: r_obs, r_capture, r_escape
        double precision, value :: step_scale
        integer, value :: samples
        double precision :: rgb_sum(3)
        double precision :: ox, oy
        double precision :: sx, sy

        integer :: i, j
        integer :: status
        integer :: red, green, blue
        double precision :: n(3), n_plane(3) 
        double precision :: tangent(2)
        double precision :: d(3)
        double precision :: rgb(3)
        real(4) :: n_plane_sp(3)
        real(4) :: phi_end

        i = (blockIdx % x - 1) * blockDim % x + threadIdx % x
        j = (blockIdx % y - 1) * blockDim % y + threadIdx % y
        if(i > width .or. j > height) then
            return
        end if

        rgb_sum = 0.0d0
        do sy = 1, samples
            do sx = 1, samples
                ox = (dble(sx) - 0.5d0) / dble(samples)
                oy = (dble(sy) - 0.5d0) / dble(samples)

                call pixel_direction(i-1,j-1,ox,oy,width,height,fov_y,n)
                call to_equatorial(n,n_plane,tangent)
                n_plane_sp = real(n_plane, 4)
                call trace_camera_ray_sp(n_plane_sp, real(r_obs,4), real(r_capture,4), real(r_escape,4), real(step_scale,4), max_steps, status, phi_end)

                if(status == 2) then
                    call from_equatorial(dble(phi_end),tangent,d)
                    call sky_color(d,rgb)
                else
                    rgb = 0.0d0
                end if

                rgb_sum = rgb_sum + rgb
            end do
        end do

        rgb = rgb_sum / dble(samples*samples)
        red   = int(255.0d0*rgb(1) + 0.5d0)
        green = int(255.0d0*rgb(2) + 0.5d0)
        blue  = int(255.0d0*rgb(3) + 0.5d0)

        pixels(i,j) = red + 256*green + 65536*blue - 16777216

    end subroutine render_kernel

    subroutine render_frame(pixels, width, height, fov_y, r_obs, r_capture, r_escape, step_scale, max_steps, samples)
        integer(c_int32_t), device :: pixels(:,:)
        integer, intent(in) :: width, height
        integer, intent(in) :: max_steps
        integer, intent(in) :: samples
        double precision, intent(in) :: fov_y
        double precision, intent(in) :: r_obs, r_capture, r_escape
        double precision, intent(in) :: step_scale

        type(dim3) :: block, grid
        integer :: ierr

        block = dim3(16,16,1)
        grid  = dim3((width + 15)/16,(height + 15)/16,1)
        call render_kernel<<<grid,block>>>(pixels, width, height, fov_y, r_obs,r_capture, r_escape, step_scale, max_steps, samples)
        ierr = cudaDeviceSynchronize()

    end subroutine render_frame

end module geodesic_kernel