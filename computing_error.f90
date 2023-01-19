module computing_error

    use parameters
    use flux

    implicit none
    
contains

    function L1_Error(Rho_ex,Rho_num) result(l1error)

        real(PR), intent(in) :: Rho_ex(0:imax+1,0:jmax+1), Rho_num(0:imax+1,0:jmax+1)
        real(PR) :: l1error
        real(PR) :: SumL1
        integer :: i, j

        ! Initialisation
        SumL1 = 0

        ! Computing L1 error on density
        do j = 1, jmax
            do i = 1, imax

                SumL1 = SumL1 + abs(Rho_ex(i,j) - Rho_num(i,j))
                
            end do
        end do

        l1error = SumL1/(imax*jmax)

    end function L1_Error

    function L2_Error(Rho_ex,Rho_num) result(l2error)

        real(PR), intent(in) :: Rho_ex(0:imax+1,0:jmax+1), Rho_num(0:imax+1,0:jmax+1)
        real(PR) :: l2error
        real(PR) :: SumL2
        integer :: i, j

        ! Initialisation
        SumL2 = 0

        ! Computing L2 error on density
        do j = 1, jmax
            do i = 1, imax

                SumL2 = SumL2 + (Rho_ex(i,j) - Rho_num(i,j))**2
                
            end do 
        end do

        l2error = sqrt(SumL2/(imax*jmax))

    end function L2_Error

    function Linf_Error(Rho_ex,Rho_num) result(linferror)

        real(PR), intent(in) :: Rho_ex(0:imax+1,0:jmax+1), Rho_num(0:imax+1,0:jmax+1)
        real(PR) :: linferror
        real(PR) :: SumLinf
        integer :: i, j

        ! Initialisation
        SumLinf = 0

        ! Computing Linf error on density
        do j = 1, jmax
            do i = 1, imax

                SumLinf = max(SumLinf,abs(Rho_ex(i,j) - Rho_num(i,j)))
                
            end do 
        end do

        linferror = SumLinf

    end function Linf_Error
    
end module computing_error
