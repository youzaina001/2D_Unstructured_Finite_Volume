program reconstruct_polynomial
    implicit none
    
    ! Declare variables
    integer :: i, j, n, m, p
    real, dimension(:), allocatable :: x, y
    real, dimension(:,:), allocatable :: A, b, c
    
    ! Read in the number of nodes and the degree of the polynomial
    read(*,*) n, p
    
    ! Allocate arrays
    allocate(x(n), y(n), A(n,p+1), b(p+1), c(p+1))
    
    ! Read in the x and y values of the nodes
    read(*,*) x, y
    
    ! Construct the matrix A and right-hand side vector b
    do i = 1, n
      A(i,1) = 1.0
      b(1) = b(1) + y(i)
      do j = 2, p+1
        A(i,j) = A(i,j-1)*x(i)
        b(j) = b(j) + A(i,j)*y(i)
      end do
    end do
    
    ! Solve the linear system using a Gaussian elimination method
    call gauss(A, b, c, n, p+1)
    
    ! Print the coefficients of the reconstructed polynomial
    write(*,*) c
  
  end program reconstruct_polynomial
  
  ! Subroutine to perform Gaussian elimination
  subroutine gauss(A, b, x, n, m)
    implicit none
    
    ! Declare variables
    integer :: i, j, k
    real, dimension(:,:), intent(in) :: A
    real, dimension(:), intent(in) :: b
    real, dimension(:), intent(out) :: x
    integer, intent(in) :: n, m
    real :: temp
    
    ! Perform Gaussian elimination
    do k = 1, n-1
      do i = k+1, n
        temp = A(i,k)/A(k,k)
        do j = k+1, m
          A(i,j) = A(i,j) - temp*A(k,j)
        end do
        b(i) = b(i) - temp*b(k)
      end do
    end do
    
    ! Perform back substitution to solve for the unknowns
    x(n) = b(n)/A(n,n)
    do i = n-1, 1, -1
      temp = 0.0
      do j = i+1, n
        temp = temp + A(i,j)*x(j)
      end do
      x(i) = (b(i) - temp)/A(i,i)
    end do
  
  end subroutine gauss
  