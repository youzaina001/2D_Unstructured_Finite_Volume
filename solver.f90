module solver

    use parameters
    use flux

    implicit none
    
contains

    subroutine spatial_discretization(Uk,Ul,Ur,Ut,Ub,phi)

        real(PR), intent(in) :: Uk(1:3), Ul(1:3), Ur(1:3), Ut(1:3), Ub(1:3)
        real(PR), intent(out) :: phi(1:3)
        real(PR) :: n_i_to_ip1(2), n_i_to_im1(2)
        real(PR) :: n_j_to_jp1(2), n_j_to_jm1(2)
        integer :: i, j

        n_i_to_im1 = (/-1,0/)
        n_i_to_ip1 = (/1,0/)
        n_j_to_jm1 = (/0,-1/)
        n_j_to_jp1 = (/0,1/)

        phi = -(1._PR/dx) * (Rusanov(Uk,Ur,n_i_to_ip1) - Rusanov(Ul,Uk,n_i_to_im1)) &
                   & -(1._PR/dy) * (Rusanov(Uk,Ut,n_j_to_jp1) - Rusanov(Ub,Uk,n_j_to_jm1))

    end subroutine spatial_discretization

    subroutine source_term(Uk,Sk)

        real(PR), intent(in) :: Uk(1:3)
        real(PR), intent(out) :: Sk(1:3)

        Sk(1) = 0._PR
        Sk(2) = -1._PR*Uk(2)
        Sk(3) = -1._PR*Uk(3)
        
    end subroutine source_term
    
end module solver
