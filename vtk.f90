module vtk

    use parameters

    implicit none
    
contains

   subroutine sortie_vtk(numero,imax,jmax,x,y,Rho,u,v)

      character(100) :: numero  
      integer, intent(in) :: imax, jmax
      real(PR), dimension(0:imax+1), intent(in) :: x
      real(PR), dimension(0:jmax+1), intent(in) :: y
      real(PR), dimension(:,:), intent(in) :: Rho, u, v

      integer :: i, j

      open(9, file = 'Output/results.'//trim(adjustl(numero))//'.vtk')
        
      write(9,'(1A26)') '# vtk DataFile Version 2.0'
      write(9,*) 'Euler 2D'
      write(9,*) 'ASCII'
      write(9,*) 'DATASET STRUCTURED_GRID'
      write(9,*) 'DIMENSIONS ',imax+1,jmax+1,1
      
    
      write(9,*) 'X_COORDINATES',imax+1,' double'
      do i=0,imax
         write(9,*) x(i)
      end do

      write(9,*) 'Y_COORDINATES',jmax+1,' double'
      do j=0,jmax
         write(9,*) y(j)
      end do

      write(9,*) 'Z_COORDINATES',1,' double'
      write(9,*) 0.0_PR
    
      write(9,*) 'CELL_DATA',imax*jmax
      write(9,*) 'SCALARS density double'
      write(9,*) 'LOOKUP_TABLE default'
      do j=1,jmax
         do i=1,imax
            write(9,*) Rho(i,j)
         end do
      end do

      write(1,*) 'VECTORS U double'
      do j=1,jmax
         do i=1,imax
            write(9,*) u(i,j),v(i,j),0.0_PR
         end do
      end do
    
      close(9)

   end subroutine sortie_vtk
    
end module vtk