module high_order

  use parameters
  !use lapack

  implicit none
  
contains

    subroutine inverse(A, Am1)

      real(PR), dimension(:,:), intent(in) :: A
      real(PR), dimension(size(A,1),size(A,2)), intent(out) :: Am1

      integer :: i, j, pivot, n
      real(PR) :: temp, pivot_val

      n = size(A,1)

      Am1 = A

      do pivot = 1, n-1

        pivot_val = Am1(pivot, pivot)
        i = pivot

        do j = pivot+1, n

          if (abs(Am1(j, pivot)) .gt. abs(pivot_val)) then

            pivot_val = Am1(j, pivot)
            i = j

          end if

        end do

        if (i .ne. pivot) then

          do j = 1, n

            temp = Am1(pivot, j)
            Am1(pivot, j) = Am1(i, j)
            Am1(i, j) = temp

          end do

        end if

        pivot_val = Am1(pivot, pivot)

        do i = pivot+1, n

          temp = Am1(i, pivot) / pivot_val

          do j = pivot+1, n

            Am1(i, j) = Am1(i, j) - temp * Am1(pivot, j)
          
          end do
          
          Am1(i, pivot) = temp
        
        end do
      
      end do

      Am1(n, n) = 1._PR / Am1(n, n)

      do i = n-1, 1, -1

        temp = 0.0

        do j = i+1, n

          temp = temp + Am1(i, j) * Am1(j, n)
        
        end do
        
        Am1(i, n) = (1._PR / Am1(i, i)) * (Am1(i, n) - temp)
      
      end do

      do i = n-1, 1, -1

        do j = i-1, 1, -1

          Am1(j, n) = Am1(j, n) - Am1(j, i) * Am1(i, n)
        
        end do
      
      end do

    end subroutine inverse


    subroutine polyomial_reconstruction_cell_IJ(dmax,indx,indy,xm,ym,value,coeffs)

      integer, intent(in) :: dmax, indx, indy
      real(PR), intent(out) ::  xm(1:imax+1), ym(1:jmax+1), value(1:imax+1,1:jmax+1)
      real(PR), dimension(:,:), allocatable, intent(out) :: coeffs
      integer :: info, ipiv(2)
      real(PR), dimension(:,:), allocatable :: RHS
      real(PR), dimension(:,:), allocatable :: P_IJ, MTM, MTB, MTM_inv

      if (dmax == 1) then

        ! Allocation des éléments du système à résoudre
        allocate(coeffs(1:2,1:1))
        allocate(P_IJ(1:9,1:2),RHS(1:9,1:1))
        allocate(MTM(1:2,1:2), MTB(1:2,1:1), MTM_inv(1:2,1:2))

        ! Allocation du système à résoudre

        P_IJ(1,1) = xm(indx-1) - xm(indx); RHS(1,1) = value(indx-1,indy-1) - value(indx,indy)
        P_IJ(1,2) = ym(indy-1) - ym(indy)

        P_IJ(2,1) = 0._PR
        P_IJ(2,2) = ym(indy-1) - ym(indy); RHS(2,1) = value(indx,indy-1) - value(indx,indy)

        P_IJ(3,1) = xm(indx+1) - xm(indx); RHS(3,1) = value(indx,indy+1) - value(indx,indy)
        P_IJ(3,2) = ym(indy-1) - ym(indy)

        P_IJ(4,1) = xm(indx-1) - xm(indx); RHS(4,1) = value(indx-1,indy) - value(indx,indy)
        P_IJ(4,2) = 0._PR

        P_IJ(5,1) = 0._PR                ; RHS(5,1) = 0._PR 
        P_IJ(5,2) = 0._PR

        P_IJ(6,1) = xm(indx+1) - xm(indx); RHS(6,1) = value(indx+1,indy) - value(indx,indy)
        P_IJ(6,2) = 0._PR

        P_IJ(7,1) = xm(indx-1) - xm(indx); RHS(7,1) = value(indx+1,indy+1) - value(indx,indy)
        P_IJ(7,2) = ym(indy+1) - ym(indy)

        P_IJ(8,1) = 0._PR
        P_IJ(8,2) = ym(indy+1) - ym(indy); RHS(8,1) = value(indx+1,indy) - value(indx,indy)

        P_IJ(9,1) = xm(indx+1) - xm(indx); RHS(9,1) = value(indx+1,indy+1) - value(indx,indy)
        P_IJ(9,2) = ym(indy+1) - ym(indy)

        ! Résolution du système ci-dessus
        ! D'abord, on calcule les matrices M^t*M et M^t*b 
        MTM = matmul(transpose(P_IJ),P_IJ)
        MTB = matmul(transpose(P_IJ),RHS)

        ! On calcule l'inverse de M^t*M
        call inverse(MTM, MTM_inv)
        !call dgetrf(2,2,MTM,2,ipiv,info)

        ! Check for errors
        !if (info /= 0) then
        !  write(*,*) "Error in LU factorization"
        !stop
        !endif

        ! Compute inverse of A
        !call dgetri(2,MTM,2,ipiv,MTM_inv,2,info)

        ! Check for errors
        !if (info /= 0) then
        !  write(*,*) "Error in inverse computation"
        !stop
        !endif

        ! Finalement, on peut obtenir les poids
        coeffs = matmul(MTM_inv,MTB)

      end if
      
    end subroutine polyomial_reconstruction_cell_IJ

    !--------------------------------------------------------------------!
    !> @brief
    !> Check if the value of the solution on a cell is physically admissible.
    !>
    !> @param[in] U_K conservative variables
    !>
    !> @author Florian Blachère
    !--------------------------------------------------------------------!
    !>function IsPhysicallyAdmissible(U_K)
    !>  real(PR), dimension(:), intent(in) :: U_K

    !>  logical :: IsPhysicallyAdmissible

    !>  if (U_K(1) < 0._PR) then
    !>    IsPhysicallyAdmissible = .FALSE.
    !>  else if (compute_pressure(U_K) < 0._PR) then
    !>    IsPhysicallyAdmissible = .FALSE.
    !>  else
    !>    IsPhysicallyAdmissible = .TRUE.
    !>  end if

    !>end function IsPhysicallyAdmissible
  
end module high_order