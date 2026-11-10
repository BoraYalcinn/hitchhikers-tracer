module geodesic
    use cudafor
    use metric, only: metric_rhs, metric_rhs_sp
    implicit none

    contains
    attributes(host,device) pure subroutine rhs(y,dy) 
        ! right hand side of the equation
        ! tells us the speed at that time  
        ! y(1)=t, y(2)=r, y(3)=θ, y(4)=φ, y(5)=p_t, y(6)=p_r, y(7)=p_θ, y(8)=p_φ
        !
        ! dx^mu/dlambda = sum_nu g^{mu nu} p_nu
        ! dp_mu/dlambda = -1/2 sum_{alpha,beta} (d_mu g^{alpha beta}) p_alpha p_beta
        !
        ! The sums are expanded by tools/gen_metric.py, so metric_rhs only
        ! evaluates the terms that are not identically zero. Looping over the
        ! full ginv/dginv arrays here cost ~80 loads and stores per call for
        ! a handful of real terms, and that array spilled to local memory on
        ! the GPU.
        double precision, intent(in)  :: y(8)
        double precision, intent(out) :: dy(8)

        double precision :: dx(4)
        double precision :: dp(4)

        call metric_rhs(y(1:4), y(5:8), dx, dp)
        dy(1:4) = dx
        dy(5:8) = dp

    end subroutine rhs
 
    attributes(host,device) pure subroutine rk4_step(y,h)
        
        double precision, intent(in)  :: h 
        double precision, intent(inout) :: y(8)

        double precision :: k1(8)
        double precision :: k2(8)
        double precision :: k3(8)
        double precision :: k4(8)
        double precision :: ytemp(8)

        call rhs(y,k1)
        ytemp = y + 0.5d0*h*k1
        call rhs(ytemp,k2)
        ytemp = y + 0.5d0*h*k2
        call rhs(ytemp,k3)
        ytemp = y + h*k3
        call rhs(ytemp,k4)

        y = y + h/6.0d0 * (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4)

    end subroutine rk4_step


    ! Single precision copies of rhs and rk4_step for the GPU renderer.
    ! Consumer GPUs run double precision at 1/64 of the single precision
    ! rate, so the kernel integrates in real(4). The double precision
    ! routines above stay the reference the tests check against.
    attributes(host,device) pure subroutine rhs_sp(y,dy)
        real(4), intent(in)  :: y(8)
        real(4), intent(out) :: dy(8)

        real(4) :: dx(4)
        real(4) :: dp(4)

        call metric_rhs_sp(y(1:4), y(5:8), dx, dp)
        dy(1:4) = dx
        dy(5:8) = dp

    end subroutine rhs_sp

    attributes(host,device) pure subroutine rk4_step_sp(y,h)

        real(4), intent(in)  :: h
        real(4), intent(inout) :: y(8)

        real(4) :: k1(8)
        real(4) :: k2(8)
        real(4) :: k3(8)
        real(4) :: k4(8)
        real(4) :: ytemp(8)

        call rhs_sp(y,k1)
        ytemp = y + 0.5e0*h*k1
        call rhs_sp(ytemp,k2)
        ytemp = y + 0.5e0*h*k2
        call rhs_sp(ytemp,k3)
        ytemp = y + h*k3
        call rhs_sp(ytemp,k4)

        y = y + h/6.0e0 * (k1 + 2.0e0 * k2 + 2.0e0 * k3 + k4)

    end subroutine rk4_step_sp


    ! Why: every physics test (photon sphere, capture, critical b, deflection)
    ! and, later, the game's photon need the same loop: integrate the ray and
    ! stop when its fate is known. Writing it once keeps the stopping rules in
    ! one place. This is the CPU double-precision reference; the GPU renderer
    ! has its own float kernel and is compared against this one.
    !
    attributes(host,device) pure subroutine trace_ray(y,h, max_steps, r_capture,r_escape, status, nsteps)

        double precision, intent(inout) :: y(8)
        double precision, intent(in) :: h 
        double precision, intent(in) :: r_capture
        double precision, intent(in) :: r_escape
        integer, intent(in) :: max_steps
        integer, intent(out)  :: status
        integer, intent(out) :: nsteps

        integer :: i
        status = 0
        nsteps = 0
        do i = 1,max_steps
            call rk4_step(y,h)
            nsteps = i
            if(y(2) < r_capture) then
                status = 1
                exit
            endif
            if(y(2) > r_escape .and. y(6) > 0.0d0) then
                status = 2
                exit
            endif
        end do

    end subroutine trace_ray


end module geodesic