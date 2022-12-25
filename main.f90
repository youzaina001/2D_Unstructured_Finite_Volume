program main

   use parameters
   use cartesian_mesh
   use flux
   use solver
   use vtk

   implicit none

   ! Déclaration des variables
   integer :: i, j, k, n
   real(PR), dimension(:), allocatable :: x, y, xm, ym
   real(PR), dimension(:,:), allocatable :: Rho, u, v
   real(PR), dimension(:,:,:), allocatable :: Un, Unp1, Sn, phi
   real(PR) :: dt
   character(100) :: numero

   ! Initialisation
   k = 0


   ! Allocation des vecteurs 1D pour le maillage
   allocate(x(1:imax+1))
   allocate(y(1:jmax+1))
   allocate(xm(1:imax+1))
   allocate(ym(1:jmax+1))

   ! Allocation des vecteurs pour la solution
   allocate(Rho(1:imax+1,1:jmax+1))
   allocate(u(1:imax+1,1:jmax+1))
   allocate(v(1:imax+1,1:jmax+1))
   allocate(Un(1:imax+1,1:jmax+1,1:3))
   allocate(Sn(1:imax+1,1:jmax+1,1:3))
   allocate(Unp1(1:imax+1,1:jmax+1,1:3))
   allocate(phi(1:imax+1,1:jmax+1,1:3))

   ! Chargement du maillage 2D cartésien
   call load_mesh(x,y,xm,ym)

   ! Initialisation
   Rho = 0._PR; u = 0._PR; v = 0._PR
   Un = 0._PR; Unp1 = 0._PR; phi = 0._PR

   ! Conditions initiales
   do j = 1, jmax + 1

      do i = 1, imax + 1

         Rho(i,j) = 1._PR
         u(i,j) = 1._PR
         v(i,j) = 1._PR
         call non_conservative_to_conservative(Rho(i,j),u(i,j),v(i,j),Un(i,j,1:3))
         
      end do

   end do

   ! Boucle en temps
   do n = 1, nmax

      dt = 1._PR*10**(-5)

      do j = 2, jmax

         do i = 2, imax

            call spatial_discretization(Un(i,j,1:3),Un(i-1,j,1:3),Un(i+1,j,1:3),Un(i,j+1,1:3),Un(i,j-1,1:3),phi(i,j,1:3))
            call source_term(Un(i,j,1:3),Sn(i,j,1:3))
            Unp1(i,j,1:3) = Un(i,j,1:3) + dt * phi(i,j,1:3) + dt * Sn(i,j,1:3)
            call conservative_to_non_conservative(Unp1(i,j,1:3),Rho(i,j),u(i,j),v(i,j))
            
         end do
         
      end do

      if (mod(n,pas_affichage) == 0) then ! On stocke la solution par fréquence de 
                                          ! 10 pas de temps  
         k=k+1

         write(numero,*) k

         if (k <= 9) then

            numero = '000'//trim(adjustl(numero))

         elseif (9 < k .and. k <= 99) then

            numero = '00'//trim(adjustl(numero))

         elseif (99 < k .and. k <= 999) then

            numero = '0'//trim(adjustl(numero))

         end if
         
         call sortie_vtk(numero,imax,jmax,x,y,Rho,u,v)

      end if

   end do
  
end program main