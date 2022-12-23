module cartesian_mesh

    use parameters

    implicit none
    
contains

    subroutine load_mesh(x,y,x_m,y_m)

        real(PR), dimension(1:imax+1), intent(out) :: x, x_m
        real(PR), dimension(1:jmax+1), intent(out) :: y, y_m

        integer :: i, j

        x = (/(-Lx/2._PR + i*dx, i = 1, imax+1)/)
        y = (/(-Ly/2._PR + j*dy, j = 1, jmax+1)/)
        x_m = (/(-Lx/2._PR + (i - 0.5_PR)*dx, i = 1, imax+1)/)
        y_m = (/(-Ly/2._PR + (j - 0.5_PR)*dy, j = 1, jmax+1)/)

    end subroutine load_mesh
    
end module cartesian_mesh