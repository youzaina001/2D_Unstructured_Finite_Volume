program main

   use parameters
   use cartesian_mesh
   use flux
   use solver
   use high_order
   use initial_condition
   use boundary_conditions
   use time_scheme
   use computing_error
   use vtk

   implicit none

   ! Déclaration des variables
   integer :: i, j, k, n, iteration
   real(PR), dimension(:), allocatable :: x, y, xm, ym
   real(PR), dimension(:,:), allocatable :: Rho, u, v, p
   real(PR), dimension(:,:), allocatable :: Rhop1, up1, vp1, sigma, sigmab
   real(PR), dimension(:,:,:), allocatable :: Un, Unp1, Sn, phi
   real(PR) :: t0, t, dt
   real(PR) :: PD, SPEED
   real(PR) :: l1error, l2error
   character(100) :: numero, num

   ! Creating the file in which we are gonna store the L1 and L2 errors for different meshes
   open(25, file='Output/Convergence_errors.dat', access = 'append')


   ! Initialisation
   k = 0
   iteration = 0
   t0 = 0._PR
   t = 0._PR
   PD = 0._PR
   SPEED = 0._PR
   l1error = 0._PR
   l2error = 0._PR

   ! Affectation des parametres
   call read_params()

   ! Allocation des vecteurs 1D pour le maillage
   allocate(x(0:imax))
   allocate(y(0:jmax))
   allocate(xm(1:imax))
   allocate(ym(1:jmax))

   ! Allocation des vecteurs pour la solution
   allocate(Rho(0:imax+1,0:jmax+1))
   allocate(u(0:imax+1,0:jmax+1))
   allocate(v(0:imax+1,0:jmax+1))
   allocate(p(0:imax+1,0:jmax+1))
   allocate(sigma(0:imax+1,0:jmax+1))
   allocate(sigmab(1:imax,1:jmax))
   allocate(Rhop1(0:imax+1,0:jmax+1))
   allocate(up1(0:imax+1,0:jmax+1))
   allocate(vp1(0:imax+1,0:jmax+1))
   allocate(Un(1:imax,1:jmax,1:3))
   allocate(Sn(1:imax,1:jmax,1:3))
   allocate(Unp1(1:imax,1:jmax,1:3))
   allocate(phi(1:imax,1:jmax,1:3))

   ! Chargement du maillage 2D cartésien
   call load_mesh(x,y,xm,ym)

   ! Initialisation
   Rho = 0._PR; u = 0._PR; v = 0._PR; p = 0._PR
   Rhop1 = 0._PR; up1 = 0._PR; vp1 = 0._PR; sigma = 0._PR; sigmab = 0._PR
   Un = 0._PR; Unp1 = 0._PR; phi = 0._PR

   ! Conditions initiales
   call init(Rho,u,v,p,xm,ym,kms,t0,Un)

   ! Boundary conditions
   if (case == 'one') then
      call Dirichlet(Rho,u,v,x,y,t,kms)
   else if (case == 'two') then
      call Neumann(Rho,u,v)
   else if (case == 'thr') then
      call Wall(Rho,u,v)
   end if

   ! Remplissange des vecteurs Un et Unp1
   do j = 1, jmax
      do i = 1, imax

         Un(i,j,1) = Rho(i,j)
         Un(i,j,2) = Rho(i,j)*u(i,j)
         Un(i,j,3) = Rho(i,j)*v(i,j)
         Unp1(i,j,1) = Un(i,j,1)
         Unp1(i,j,2) = Un(i,j,2)
         Unp1(i,j,3) = Un(i,j,3)

      end do
   end do

   ! Stockage de la solution exacte
   write(num,*) 112233445
   call sortie_vtk(num,imax,jmax,x,y,Rho(1:imax,1:jmax),u(1:imax,1:jmax),v(1:imax,1:jmax))

   ! Boucle en temps
   do while (t <= tf)

      ! Computing time step
      call time_step(Rho,u,v,xm,ym,t,kms,dt,sigmab)

      print*, "The time step is equal to:", dt, "."

      ! Paramètre thete_e permettant la limitation
      print*, "La valeur de sigma est:", maxval(sigmab)
      theta_e = 1._PR!/(1._PR + maxval(sigmab)*t)

      ! Solution using an explicit Euler solver or RK3 scheme
      call explicit_euler(Un,t,dt,kms,xm,ym,Unp1)

      ! Mise a jour des solutions
      Un = Unp1

      ! Boundary conditions
      if (case == 'one') then
         call Dirichlet(Rho,u,v,x,y,t,kms)
      else if (case == 'two') then
         call Neumann(Rho,u,v)
      else if (case == 'thr') then
         call Wall(Rho,u,v)
      end if

      do j = 1, jmax

         do i = 1, imax

            ! Extracting non-conservative variables
            call conservative_to_non_conservative(Un(i,j,1:3),Rho(i,j),u(i,j),v(i,j))
            
         end do
         
      end do

      if (mod(iteration,pas_affichage) == 0) then  

         k=k+1
         write(numero,*) k
         if (k <= 9) then
            numero = '000'//trim(adjustl(numero))
         elseif (9 < k .and. k <= 99) then
            numero = '00'//trim(adjustl(numero))
         elseif (99 < k .and. k <= 999) then
            numero = '0'//trim(adjustl(numero))
         end if
         call sortie_vtk(numero,imax,jmax,x,y,Rho(1:imax,1:jmax),u(1:imax,1:jmax),v(1:imax,1:jmax))

      end if

      iteration = iteration + 1
      t = t + dt
      print*, "On est à la ", iteration, "itération."
      print*, "Le temps est égal à :", t, "."

   end do

   ! Calcul des erreurs L1 et L2
   !l1error = L1_Error(Rho_ex,Rho)
   !l2error = L2_Error(Rho_ex,Rho)
   !write(25,*) l1error, l2error

   deallocate(x,y,xm,ym)
   deallocate(Rho,u,v,p)
   deallocate(sigma,Rhop1,up1,vp1)
   deallocate(Un,Sn,Unp1,phi)
  
end program main