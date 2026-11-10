! ray_plane: rotates a camera ray into the equatorial plane and back.
!
! Why this module exists: the geodesic integrator only follows rays in the
! equatorial plane (theta = pi/2, p_theta = 0), but a pixel in the upper half of
! the screen looks out of that plane. A Schwarzschild black hole is spherically
! symmetric, so every ray stays forever in one plane through the hole's centre,
! the plane spanned by the observer, the centre and the ray direction. We rotate
! that plane onto the equator, integrate in 2D, and rotate the result back.
! This module is only that geometry; it contains no physics.
module ray_plane
    use cudafor
    implicit none

    ! Below this the tangential part of the direction is treated as zero
    ! (the ray points straight at the hole, so its plane is undefined).
    double precision, parameter :: tangent_epsilon = 1.0d-14

contains

    ! Why: init_ray needs a direction whose theta component is zero, while the
    ! camera gives a general direction n = (n_r, n_theta, n_phi).
    !
    ! The radial part n_r is unchanged by the rotation. The two tangential parts
    ! are merged into one length n_t and placed entirely on the phi axis:
    !   n_t = sqrt(n_theta^2 + n_phi^2),   n_plane = (n_r, 0, n_t)
    ! The unit direction of the tangential part, tangent = (n_theta, n_phi)/n_t,
    ! is saved because from_equatorial needs it to undo the rotation.
    !
    ! Example: n = (-0.6, 0.48, 0.64) gives n_t = 0.8, n_plane = (-0.6, 0, 0.8),
    ! tangent = (0.6, 0.8).
    !
    !   n        in : local camera direction (n_r, n_theta, n_phi), unit length
    !   n_plane  out: direction in the equatorial plane (n_r, 0, n_t)
    !   tangent  out: unit tangential direction (theta part, phi part)
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
    !
    ! Axes of the ray's plane in Cartesian coordinates (observer at theta = pi/2,
    ! phi = 0, so r-hat = x, theta-hat = -z, phi-hat = y):
    !   e1 = (1, 0, 0)                       toward the observer
    !   e2 = (0, tangent(2), -tangent(1))    along the saved tangential direction
    ! The final direction is  d = cos(phi_end) e1 + sin(phi_end) e2.
    !
    ! Example (no bending): tangent = (0.6, 0.8), phi_end = atan2(0.8, -0.6)
    ! gives d = (-0.6, 0.64, -0.48), which is the original camera direction
    ! written in Cartesian axes.
    !
    !   phi_end  in : final phi of the ray, y(4) after integration
    !   tangent  in : saved by to_equatorial
    !   d        out: sky direction (x, y, z), unit length
    attributes(host,device) subroutine from_equatorial(phi_end, tangent, d)
        double precision, intent(in) :: phi_end
        double precision, intent(in) :: tangent(2)
        double precision, intent(out) :: d(3)

        d(1) = cos(phi_end)
        d(2) = sin(phi_end)*tangent(2)
        d(3) = -sin(phi_end)*tangent(1)
    end subroutine from_equatorial

end module ray_plane