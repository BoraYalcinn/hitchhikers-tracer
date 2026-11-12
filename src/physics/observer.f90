module observer
    use metric, only: metric_cov
    implicit none

    contains

    attributes(host,device) pure subroutine static_tetrad(x,e)
        double precision, intent(in) :: x(4)
        double precision, intent(out):: e(4,4)

        double precision :: gcov(4,4)
        integer :: a

        call metric_cov(x,gcov)

        e = 0.0d0
        do a = 1,4
            e(a,a) = 1.0d0 / sqrt(abs(gcov(a,a)))
        end do
    end subroutine static_tetrad

endmodule observer
