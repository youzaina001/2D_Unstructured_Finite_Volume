module cartesian_mesh

    use parameters

    implicit none
    
contains

    subroutine load_mesh(x,y,x_m,y_m)

        real(PR), dimension(0:imax), intent(out) :: x
        real(PR), dimension(0:jmax), intent(out) :: y
        real(PR), dimension(1:imax), intent(out) :: x_m
        real(PR), dimension(1:imax), intent(out) :: y_m

        integer :: i, j

        x = (/(i*dx, i = 0, imax)/)
        y = (/(j*dy, j = 0, jmax)/)
        x_m = (/((i - 0.5_PR)*dx, i = 1, imax)/)
        y_m = (/((j - 0.5_PR)*dy, j = 1, jmax)/)

    end subroutine load_mesh
    
end module cartesian_mesh