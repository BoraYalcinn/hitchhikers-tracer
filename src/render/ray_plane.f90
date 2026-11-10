module ray_plane
    use cudafor
    implicit none

    ! Below this the tangential part of the direction is treated as zero
    ! (the ray points straight at the hole, so its plane is undefined).
    double precision, parameter :: tangent_epsilon = 1.0d-14

contains

    ! Why: init_ray needs a direction whose theta component is zero, while the
    ! camera gives a general direction n = (n_r, n_theta, n_phi).
    attributes(host,device) subroutine to_equatorial(n, n_plane, tangent)
        double precision, intent(in) :: n(3)
        double precision, intent(out) :: n_plane(3)
        double precision, intent(out) :: tangent(2)

        double precision :: n_t

        n_t = sqrt(n(2)**2 + n(3)**2)

        n_plane(1) = n(1)
        n_plane(2) = 0.0d0
        n_plane(3) = n_t

        if (n_t > tangent_epsilon) then
            tangent(1) = n(2)/n_t
            tangent(2) = n(3)/n_t
        else
            ! Straight radial ray: any tangential direction gives the same
            ! result, so pick (0, 1) to avoid dividing by zero.
            tangent(1) = 0.0d0
            tangent(2) = 1.0d0
        end if
    end subroutine to_equatorial

    ! Why: after integrating in the equatorial plane we only know how far the
    ! ray swung around, the angle phi_end. The shader needs a direction in the
    ! real 3D sky to look up the background, so the rotation is undone here.
    attributes(host,device) subroutine from_equatorial(phi_end, tangent, d)
        double precision, intent(in) :: phi_end
        double precision, intent(in) :: tangent(2)
        double precision, intent(out) :: d(3)

        d(1) = cos(phi_end)
        d(2) = sin(phi_end)*tangent(2)
        d(3) = -sin(phi_end)*tangent(1)
    end subroutine from_equatorial

end module ray_plane