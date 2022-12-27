module computing_error

    use parameters

    implicit none
    
contains

    function L1_Error(Uex,Unum) result(l1error)

        real(PR), intent(in) :: Uex(1:imax+1,1:jmax+1,1:3), Unum(1:imax+1,1:jmax+1,1:3)
        real(PR) :: Rho_ex(1:imax+1,1:jmax+1), Rho_num(1:imax+1,1:jmax+1), u_ex(1:imax+1,1:jmax+1)
        real(PR) :: u_num(1:imax+1,1:jmax+1), v_ex(1:imax+1,1:jmax+1), v_num(1:imax+1,1:jmax+1)
        real(PR) :: l1error
        real(PR) :: SumL1
        integer :: i, j

        ! Initialisation
        Sum = 0

        ! Extracting non conservative variables
        do j = 1, jmax + 1

            do i = 1, imax + 1

                call conservative_to_non_conservative(Uex(i,j,1:3),Rho_ex(i,j),u_ex(i,j),v_ex(i,j))
                call conservative_to_non_conservative(Unum(i,j,1:3),Rho_num(i,j),u_num(i,j),v_num(i,j))
                
            end do
            
        end do

        ! Computing L1 error on density
        do j = 1, jmax + 1

            do i = 1, imax + 1

                SumL1 = SumL1 + abs(Rho_ex(i,j) - Rho_num(i,j))
                
            end do
            
        end do

        l1error = SumL1

    end function L1_Error

    function L2_Error(Uex,Unum) result(l2error)

        real(PR), intent(in) :: Uex(1:imax+1,1:jmax+1,1:3), Unum(1:imax+1,1:jmax+1,1:3)
        real(PR) :: Rho_ex(1:imax+1,1:jmax+1), Rho_num(1:imax+1,1:jmax+1), u_ex(1:imax+1,1:jmax+1)
        real(PR) :: u_num(1:imax+1,1:jmax+1), v_ex(1:imax+1,1:jmax+1), v_num(1:imax+1,1:jmax+1)
        real(PR) :: l2error
        real(PR) :: SumL2
        integer :: i, j

        ! Initialisation
        Sum = 0

        ! Extracting non conservative variables
        do j = 1, jmax + 1

            do i = 1, imax + 1

                call conservative_to_non_conservative(Uex(i,j,1:3),Rho_ex(i,j),u_ex(i,j),v_ex(i,j))
                call conservative_to_non_conservative(Unum(i,j,1:3),Rho_num(i,j),u_num(i,j),v_num(i,j))
                
            end do
            
        end do

        ! Computing L1 error on density
        do j = 1, jmax + 1

            do i = 1, imax + 1

                SumL2 = SumL2 + (Rho_ex(i,j) - Rho_num(i,j))**2
                
            end do
            
        end do

        l2error = SumL2

    end function L2_Error
    
end module computing_error
