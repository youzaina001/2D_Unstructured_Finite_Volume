module solver

    use parameters
    use flux

    implicit none
    
contains

    subroutine spatial_discretization(Uk,Ul,Ur,Ut,Ub,Pk,Pl,Prr,Pt,Pb,PPk,PPl,PPr,PPt,PPb,phi)

        real(PR), intent(in) :: Uk(1:3), Ul(1:3), Ur(1:3), Ut(1:3), Ub(1:3)
        real(PR), intent(in) :: Pk, Pl, Prr, Pt, Pb, PPk, PPl, PPr, PPt, PPb
        real(PR), intent(out) :: phi(1:3)
        real(PR) :: n_i_to_ip1(1:2), n_i_to_im1(1:2)
        real(PR) :: n_j_to_jp1(1:2), n_j_to_jm1(1:2)
        integer :: i, j

        n_i_to_im1(1) = -1._PR; n_i_to_im1(2) = 0._PR
        n_i_to_ip1(1) = 1._PR; n_i_to_ip1(2) = 0._PR
        n_j_to_jm1(1) = 0._PR; n_j_to_jm1(2) = -1._PR
        n_j_to_jp1(1) = 0._PR; n_j_to_jp1(2) = 1._PR

        phi = -(1._PR/dx) * (Rusanov(Uk,Ur,Pk,Prr,PPk,PPr,n_i_to_ip1) + Rusanov(Ul,Uk,Pl,Pk,PPl,PPk,n_i_to_im1)) &
                   & -(1._PR/dy) * (Rusanov(Uk,Ut,Pk,Pt,PPk,PPt,n_j_to_jp1) + Rusanov(Ub,Uk,Pb,Pk,PPb,PPk,n_j_to_jm1))

    end subroutine spatial_discretization

    subroutine source_term(Uk,Sk)

        real(PR), intent(in) :: Uk(1:3)
        real(PR), intent(out) :: Sk(1:3)

        Sk(1) = 0._PR
        Sk(2) = -1._PR*Uk(2)
        Sk(3) = -1._PR*Uk(3)
        
    end subroutine source_term
    
end module solver
