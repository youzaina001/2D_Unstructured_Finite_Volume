module time_scheme

    use parameters
    use flux
    use solver
    use initial_condition

    implicit none
    
contains

    subroutine explicit_euler(Un,t,dt,kms,xm,ym,Unp1)

        integer, intent(in) :: kms 
        real(PR), intent(in) :: Un(1:imax+1,1:jmax+1,1:3), xm(1:imax+1), ym(1:jmax+1), t, dt
        real(PR), intent(out) ::  Unp1(1:imax+1,1:jmax+1,1:3)
        real(PR) :: Sn(1:imax+1,1:jmax+1,1:3), phi(1:imax+1,1:jmax+1,1:3), sigma(2:imax,2:jmax)
        integer :: i, j

        do j = 2, jmax

            do i = 2, imax
   
                call spatial_discretization(Un(i,j,1:3),Un(i-1,j,1:3),Un(i+1,j,1:3),Un(i,j+1,1:3),Un(i,j-1,1:3),phi(i,j,1:3))
                call source_term(Un(i,j,1:3),Sn(i,j,1:3))
                sigma(i,j) = compute_sigma(xm(i),ym(j),t,kms)
                Unp1(i,j,1:3) = Un(i,j,1:3) + dt * phi(i,j,1:3) - dt * sigma(i,j) * Sn(i,j,1:3)
               
            end do
            
        end do
        
    end subroutine explicit_euler
    
end module time_scheme
