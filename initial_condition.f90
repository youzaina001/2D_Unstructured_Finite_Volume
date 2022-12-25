module initial_condition

    use parameters
    use flux

    implicit none
    
contains

    subroutine init(Rho,u,v,xm,ym,UU,case)

        character(3), intent(in) :: case
        real(PR), intent(in) :: xm(1:imax+1), ym(1:jmax+1)
        real(PR), intent(inout) :: Rho(1:imax+1,1:jmax+1), u(1:imax+1,1:jmax+1), v(1:imax+1,1:jmax+1)
        real(PR), intent(out) ::  UU(1:imax+1,1:jmax+1,1:3)
        integer :: i, j

        if (case == 'one') then

            do j = 1, jmax + 1

                do i = 1, imax + 1

                    if ((xm(i)-0.5_PR)**2 + (ym(j)-0.4_PR)**2 <= 0.3_PR**2) then

                        Rho(i,j) = 5._PR
                        u(i,j) = 0._PR
                        v(i,j) = 0._PR

                    else

                        Rho(i,j) = 1._PR
                        u(i,j) = 0._PR
                        v(i,j) = 0._PR
                        
                    end if
                    
                    call non_conservative_to_conservative(Rho(i,j),u(i,j),v(i,j),UU(i,j,1:3))
                   
                end do
          
            end do
            
        end if
    
        
    end subroutine init

end module initial_condition
